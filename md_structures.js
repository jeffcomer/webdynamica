//////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Jeffrey Robert Comer
//
// This file is part of WebDynamica, a browser-based interactive
// molecular dynamics program using WebGL
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <https://www.gnu.org/licenses/>.
//////////////////////////////////////////////////////////////////////

'use strict';

// Define the different kinds of terms
class TermInfo {
    static bond = 0;
    static angle = 1;
    static dihedral = 2;
    static exclude = 3;
    static lj = 4;
    
    static num = 5;
    static name = ['bond', 'angle', 'dihedral', 'exclusion', 'lj'];
    static parFields = [2, 2, 3, 2, 2];
    static bodies = [2, 3, 4, 2, 1];
    // Bonds and excludes have just one pixel per term (TermInfo.pixels[termKind] == 1)
    // Angles and dihedrals have two pixels per term (TermInfo.pixels[termKind] == 2)
    static pixels = [1, 2, 2, 1, 1];
}

class MaterialEnum {
    static normal = 0;
    static water = 1;
    static background = 2;
}

class TermTexInfo {
    constructor(termKind, atomWidth, atomHeight, termsPerAtom) {
	this.width = termsPerAtom*TermInfo.pixels[termKind]*atomWidth;
	this.height = atomHeight;
	this.size = this.width*this.height;
	this.fields = 4*this.size; // RGBA
	this.termsPerAtom = termsPerAtom;
    }
}


// Store atom data to be inserted into WebGL textures
class AtomTextureData {
    // baseMolecule: a molecular structure of the form produced by charmm_to_javascript.js
    // 
    // termsPerAtom: a 4-member list with the maximum number of
    // bonds, angles, dihedals, and exclusions per atom
    //
    // additionalCapacity: minimum extra atom spaces in texture map beyond
    // the number of atoms in the baseMolecule
    constructor(baseMolecule, termsPerAtom, additionalCapacity = 0) {
	// Set the target number of atoms represented in the texture
	const targetNum = baseMolecule.num + additionalCapacity;
	
	// Because the width of the texture will be expanded in the bond, angle, dihedral textures,
	// we choose a smaller width for the atom texture so that none of the texture dimensions get too big
	const heightScale = Math.max(2,...termsPerAtom);
	console.log(`Texture height scale factor: ${heightScale}`);
	
	// Calculate the dimensions for the atom texture
	// The goal is to make the largest term texture nearly square
	this.width = Math.ceil(Math.sqrt(targetNum/heightScale));
	this.height = Math.ceil(targetNum/this.width);
	this.size = this.width*this.height;
	this.fields = 4*this.size; // RGBA
	//console.log(`${baseMolecule.name} has ${baseMolecule.num} atoms`);
	console.log(`Texture provides space for ${this.size} atoms`);
	console.log(`Initially active atoms: ${baseMolecule.num}`);
	console.log(`Initially inactive atoms: ${this.size - baseMolecule.num}`);
	console.log(`atom texture size: ${this.width} ${this.height}`);

	// Maximum number of terms per atom
	this.termsPerAtom = [termsPerAtom[TermInfo.bond], termsPerAtom[TermInfo.angle],
			     termsPerAtom[TermInfo.dihedral], termsPerAtom[TermInfo.exclude]];

	// Special residue names
	this.waterResNames = ['TIP3', 'TIP4', 'SPCE', 'SPCF', 'WAT', 'HOH'];
	this.backgroundResNames = ['C', 'GRPH', 'CELU', 'QTZ'];
	
	// Create the texture data 
	// Define the dimensions of the textures
	// The size of these textures is similar to the atom texture, but
	// scaled in x by the number of terms per atom (maximum necessary)
	// and the number of pixels per term.
	// Bonds and excludes have just one pixel per term (TermInfo.pixels[termKind] == 1)
	// Angles and dihedrals have two pixels per term (TermInfo.pixels[termKind] == 2)
	this.termDim = [];
	for (let k = 0; k <= TermInfo.exclude; k++) {
	    this.termDim[k] = new TermTexInfo(k, this.width, this.height, this.termsPerAtom[k]);
	    console.log(`${TermInfo.name[k]} texture size: ${this.termDim[k].width} ${this.termDim[k].height}`);
	}
	// Some aliases for the dimensions
	this.atomDim = {width: this.width, height: this.height, size: this.size, fields: this.fields}; 
	this.bondDim = this.termDim[TermInfo.bond];
	this.angleDim = this.termDim[TermInfo.angle];
	this.dihedralDim = this.termDim[TermInfo.dihedral];
	this.excludeDim = this.termDim[TermInfo.exclude];

	// Create the texture data arrays
	this.pos = new Float32Array(this.fields); // [x y z active]
	this.bond = new Float32Array(this.bondDim.fields); // [i1 0 K bondLen]
	this.angle = new Float32Array(this.angleDim.fields); // [i1 i2 0 0] [role K theta 0]
	this.dihedral = new Float32Array(this.dihedralDim.fields); // [i1 i2 i3 0] [role K mult delta]
	this.exclude = new Float32Array(this.excludeDim.fields);  // [i1 0 eps RMin]
	this.nonbond = new Float32Array(this.fields); // [mass charge eps RMin]
	this.color = new Float32Array(this.fields); // [red green blue radius]
	this.select = new Float32Array(this.fields); // [fragment monomer material radius]
	this.restrain = new Float32Array(this.fields); // [x0 y0 z0 spring]

	// Do not allow the default mass to be zero
	for (let ai = 0; ai < this.size; ai++) {
	    this.nonbond[4*ai] = 12.0;
	}
	
	// Create additional arrays
	this.typeList = new Array(this.size).fill('UNK');
	this.resNameList = new Array(this.size).fill('unknown');
	this.nameList = new Array(this.size).fill('UNK');

	// Insert the base molecule into the texture data arrays
	this.hashDim = {};
	this.insertMolecule(baseMolecule, this.pos);
	// Hash table is initialized in insertMolecule
    }

    static countActive(posData) {
	const size = Math.floor(posData.length/4);
	let activeNum = 0;

	// Active atoms have a nonzero w-component
	for (let ai = 0; ai < size; ai++) {
	    if (posData[4*ai+3] != 0.0) activeNum++;
	}
	return activeNum;
    }

    static getCenter(mol) {
	if (mol.index_mass != null) {
	    // Geometric centroid
	    let sum = [0.0, 0.0, 0.0];
	    for (let i = 0; i < mol.atom_coord.length; i++) {
		const j = i % 3;
		sum[j] += mol.atom_coord[i];
	    }
	    return [sum[0]/mol.num, sum[1]/mol.num, sum[2]/mol.num];
	} else {
	    // Center of mass
	    let sum = [0.0, 0.0, 0.0];
	    let weightSum = 0.0;
	    for (let ai = 0; ai < mol.num; ai++) {
		const mass = mol.par_mass[mol.index_mass[ai]];
		sum[0] += mass*mol.atom_coord[3*ai];
		sum[1] += mass*mol.atom_coord[3*ai+1];
		sum[2] += mass*mol.atom_coord[3*ai+2];
		weigthSum += mass;
	    }
	    return [sum[0]/mass, sum[1]/mass, sum[2]/mass];
	}
    }

    setPositions(posData) {
	// Validate the input
	if (posData.length != this.pos.length) {
	    console.log(`Warning! setPositions(): Input position data for texture is not equal to the number of fields in the AtomTextureInfo object (${posData.length} != ${this.pos.length})`);
	    return false;
	}

	// Update the position data
	for (let i = 0; i < this.pos.length; i++)
	    this.pos[i] = posData[i];
	return true;
    }

    // Delete a molecule from the textures
    deleteMolecule(lastPosData, selectId, selectColumn) {
	if (selectColumn >= 4) {
	    console.log(`Warning! deleteMolecule() invalid column: ${selectColumn} (must be 0–3 inclusive)`);
	    return 0;
	}
	
	// Update the positions
	this.setPositions(lastPosData);

	// Make a list of atoms to delete
	const deleteList = [];
	for (let ai = 0; ai < this.size; ai++) {
	    const ourId= this.select[4*ai + selectColumn];
	    if (ourId == selectId) {
		deleteList.push(ai);
	    }
	}
	//console.log(`Deleting ${deleteList.length} atoms`);

	// Zero out the atom texture data for these selections
	for (const ai of deleteList) {
	    for (let c = 0; c < 4; c++) {
		this.pos[4*ai+c] = 0.0;
		this.color[4*ai+c] = 0.0;
		this.select[4*ai+c] = 0.0;
		this.restrain[4*ai+c] = 0.0;
	    }
	}
	// Same for nonbond, but don't set the mass to zero
	for (const ai of deleteList) {
	    this.nonbond[4*ai] = 12.0;
	    for (let c = 1; c < 3; c++) {
		this.nonbond[4*ai+c] = 0.0;
	    }
	}

	// Delete all terms associated with this atom
	for (let k = 0; k <= TermInfo.exclude; k++) {
	    const termCount = this.deleteTerm(deleteList, k);
	    //console.log(`Deleted ${this.termsPerAtom[k]*deleteList.length} ${TermInfo.name[k]} terms for atoms being deleted`);
	    //console.log(`Deleted ${termCount} ${TermInfo.name[k]} terms altogether (including for other atoms)`);
	}

	// Clear any selected items in the hash table
	// This needs to be done since subsequently inserted atoms might
	// have the same indices.
	// However, it means that when the hash table is searched, we
	// can't stop at an empty node and need to go all the way to
	// maxCollisions. This is fine since the SIMD nature of the
	// shader means there is no advantage to shortcutting the search
	// anyway.
	let hashDelete = 0;
	for (let k = 0; k <= this.hashDim.size; k++) {
	    const atomIndex0 = this.hash[4*k];
	    const atomIndex1 = this.hash[4*k+1];
	    const RMin = this.hash[4*k+3]; 

	    // Delete active exclusions containing this atom
	    if (RMin != 0.0 && (deleteList.includes(atomIndex0) || deleteList.includes(atomIndex1))) {
		this.hash[4*k] = 0.0;
		this.hash[4*k+1] = 0.0;
		this.hash[4*k+2] = 0.0;
		this.hash[4*k+3] = 0.0;
		hashDelete++;
	    }
	}
	//console.log(`Deleted ${hashDelete} exclusions from the hash table.`);

	return deleteList.length;
    }

    // Delete a bonded term including all terms for the selected atoms
    // and all terms involving the selected atoms for other atoms
    deleteTerm(deleteList, termKind) {
	const termsPerAtom = this.termsPerAtom[termKind];
	const termFields = 4*TermInfo.pixels[termKind];
	const fieldsPerAtom = termsPerAtom*termFields;
	let termCount = 0;

	let texData = null;
	if (termKind == TermInfo.bond) {
	    texData = this.bond;
	} else if (termKind == TermInfo.angle) {
	    texData = this.angle;
	} else if (termKind == TermInfo.dihedral) {
	    texData = this.dihedral;
	} else if (termKind == TermInfo.exclude) {
	    texData = this.exclude;
	} else {
	    console.log(`Warning! deleteBondedTerm(): Unrecognized termKind ${termKind}`);
	    return 0;
	}

	// Delete all terms pertaining to the selected atoms
	for (const ai of deleteList) {
	    const start = ai*fieldsPerAtom;
	    for (let i = 0; i < fieldsPerAtom; i++) {
		texData[start + i] = 0.0;
	    }
	    termCount += termsPerAtom;
	}

	// Delete all terms that contain this atom
	const totalTerms = this.size * termsPerAtom;
	for (let t = 0; t < totalTerms; t++) {
	    const termStart = t*termFields;
	    
	    // Does the index list contain an atom being deleted?
	    // Atom indices are in the first pixel from zero to bodies-1
	    const neighNum = TermInfo.bodies[termKind] - 1;
	    let deleting = false;
	    for (let c = 0; c < neighNum; c++) {
		const neighIndex = texData[termStart + c];
		if (deleteList.includes(neighIndex)) {
		    deleting = true;
		    break;
		}
	    } // loop over neighbor atom indices

	    // Delete this term if appropriate
	    if (deleting) {
		termCount++;
		for (let c = 0; c < termFields; c++) {
		    texData[termStart + c] = 0.0;
		} // Loop over term fields
	    }
	} // Loop over all terms

	return termCount;
    }
    
    insertMolecule(mol, lastPosData, shiftCenter, consecutive = false) {
	const centering = (shiftCenter != null);
	if (mol.num == null || mol.name == null) {
	    console.log(`Warning! Could not insert molecule. mol object is invalid.`);
	    return false;
	}
	if (mol.index_type == null || mol.par_type == null) {
	    console.log(`Warning! Could not insert molecule. Molecule ${mol.name} has no atom types.`);
	    return false;
	}
	if (mol.index_type.length != mol.num) {
	    console.log(`Warning! Could not insert molecule ${mol.name} because types are not provided for all atoms (${mol.index_type.length} != ${mol.num})`);
	    return false;
	}
	
	// Check if there is enough space in the textures to insert this molecule
	const currNum = AtomTextureData.countActive(lastPosData);
	const avail = this.size - currNum;
	console.log('currNum: ', currNum);
	console.log('avail: ', avail);
	if (avail < mol.num) {
	    console.log(`Warning! Could not insert molecule ${mol.name} because there isn't enough space in the textures (space avail: ${avail}, space needed: ${mol.num}). Molecule not added.`);
	    return false;
	}

	// Make a table for mapping from atom indices in mol to empty texture indices
	const map = new Int32Array(mol.num);
	if (consecutive) {
	    // Find a sufficiently large consecutive space in the texture
	    let start = -1;
	    let size = 0;
	    for (let ai = 0; ai < this.size; ai++) {
		if (lastPosData[4*ai+3] != 0.0) {
		    // Non-empty
		    start = -1;
		} else {
		    // Set a new start position if this is the first empty slot since the last occupied slot
		    if (start < 0) start = ai;

		    // Some molecules have length one, so it's okay if we just set start above.
		    const size = ai - start + 1;
		    // We've found a sufficiently large consecutive region
		    if (size >= mol.num) break;
		}
	    }
	    // Not found.
	    if (start < 0) {
		console.log(`Warning! Could not insert molecule ${mol.name} because there isn't a large enough consecutive inactive region in the texture. Molecule not added.`);
		return false;
	    }

	    // Create the atom index map by shifting from start
	    // Map from mol indices to texture indices
	    for (let i = 0; i < mol.num; i++) map[i] = i + start;
	} else {
	    // Fill in any empty spaces.
	    let currAtomIndex = 0;
	    for (let i = 0; i < mol.num; i++) {
		// Skip forward in the texture until we find the next empty space
		while (lastPosData[4*currAtomIndex+3] != 0.0) {
		    currAtomIndex++;
		    if (4*currAtomIndex >= lastPosData.size) {
			// This shouldn't happen since we counted above.
			console.log(`Warning! Failed to insert molecule because we ran out of space in the texture`);
			return false;
		    }
		}

		// Map from mol indices to texture indices
		map[i] = currAtomIndex;
		currAtomIndex++;
	    }
	}

	// We have decided to insert the molecule
	// First, we should update the position data
	this.setPositions(lastPosData);
	
	// Insert the new position data from the molecule
	// We shift the center of mass to "shiftCenter"
	const shift = [0, 0, 0];
	// Are we centering the structure about it's centroid?
	if (centering) {
	    const center = AtomTextureData.getCenter(mol);
	    for (let c = 0; c < 3; c++) shift[c] = shiftCenter[c] - center[c];
	}
	for (let ai = 0; ai < mol.num; ai++) {
	    // Set the coordinates
	    for (let c = 0; c < 3; c++) 
		this.pos[4*map[ai]+c] = mol.atom_coord[3*ai+c] + shift[c];

	    // The atom is active
	    this.pos[4*map[ai] + 3] = 1.0;
	}

	// Update the type array
	for (let ai = 0; ai < mol.num; ai++) 
	    this.typeList[map[ai]] = mol.par_type[mol.index_type[ai]];
	
	// Update the resName array
	if (mol.index_resname != null) {
	    for (let ai = 0; ai < mol.num; ai++)
		this.resNameList[map[ai]] = mol.par_resname[mol.index_resname[ai]];
	} else {
	    for (let ai = 0; ai < mol.num; ai++)
		this.resNameList[map[ai]] = 'unknown';
	}

	// Update the name array
	if (mol.index_name != null) {
	    for (let ai = 0; ai < mol.num; ai++)
		this.nameList[map[ai]] = mol.par_name[mol.index_name[ai]];
	} else {
	    // If there are no atom names, use the types
	    for (let ai = 0; ai < mol.num; ai++)
		this.nameList[map[ai]] = this.typeList[map[ai]];
	}

	// Insert data into the nonbonded texture
	this.insertNonbond(mol, map);
	// Set the color texture using the nonbond texture, type, and resName
	this.setColors(mol, map);
			
	// Insert data into the bond, angle, dihedral, and exclusion textures
	this.insertTerms(mol, map, TermInfo.bond);
	this.insertTerms(mol, map, TermInfo.angle);
	this.insertTerms(mol, map, TermInfo.dihedral);
	this.insertTerms(mol, map, TermInfo.exclude);

	// Insert data into the selection texture (fragment number, monomer number, material)
	this.insertSelect(mol, map);
	// Insert the initial positions into the restraint texture
	this.insertRestrain(mol, map);

	// Insert into the hash table or generate it from scratch
	if (this.hash != null) {
	    // If the hash table exists, insert the molecule
	    // Note that we may increase maxCollisions!
	    // This may require recompiling the shaders
	    this.insertMoleculeIntoHashTable(map);
	} else {
	    // Optimize the hash table from scratch
	    this.makeHashTable();
	}
	
	return true;
    }

    // Stores [red, green, blue, radius] in the texture data array
    setColors() {
	// Reset the colors for all atoms
	for (let ai = 0; ai < this.typeList.length; ai++) {
	    const type = this.typeList[ai];
	    const resName = this.resNameList[ai];
	    
	    let color = null;
	    if (this.backgroundResNames.some((rn) => rn == resName ) && type[0] == 'C') {
		// Gray for graphene carbon
		color = [0.75,0.75,0.75];
	    } else {
		// Determine the color based on the type
		color = AtomTextureData.getColor(type);
	    }
	    this.color[4*ai] = color[0];
	    this.color[4*ai+1] = color[1];
	    this.color[4*ai+2] = color[2];
	    // This field will store half the Lennard-Jones Rmin (not the alpha)
	    let radius = this.nonbond[4*ai+3];
	    // Don't show super-small radii for H-bonding hydrogen
	    if (radius < 2.4) radius = 2.4;

	    // Inactive atoms get a zero radius
	    const active = this.pos[4*ai+3];
	    this.color[4*ai+3] = active ? 0.5*radius : 0.0;
	}
    }

    // Restrain atoms marked as the background in "select" with spring constant
    restrainBackground(spring) {
	let n = 0;
	for (let ai = 0; ai < this.size; ai++) {
	    const material = this.select[4*ai+2];
	    const mass = this.nonbond[4*ai]; // Don't restrain hydrogen

	    if (material == MaterialEnum.background && mass > 4.0) {
		this.restrain[4*ai+3] = spring;
		n++;
	    }
	}
	return n;
    }

    // Add newly inserted molecule to the restraint texture (with no restraints)
    insertRestrain(mol, map) {
	// Assume that we've already inserted into this.pos
	for (let ai = 0; ai < map.length; ai++) {
	    const texAtomIndex = map[ai];
	    for (let c = 0; c < 3; c++) {
		// Get the original position
		const orig = mol.atom_coord[3*ai+c];
		this.restrain[4*texAtomIndex+c] = this.pos[4*texAtomIndex+c];
	    }
	    // Default no restraints (spring constant of zero)
	    this.restrain[4*texAtomIndex + 3] = 0.0;
	}
    }

    // Add newly inserted molecule to the selection texture
    insertSelect(mol, map) {
	// Zero out the terms for this atom
	for (let ai = 0; ai < mol.num; ai++) {
	    const texAtomIndex = map[ai];
	    for (let c = 0; c < 4; c++) {
		this.select[4*texAtomIndex + c] = 0.0;
	    }
	}
	
	// Get the current maximum fragment number
	let maxFragment = 0;
	for (let ai = 0; ai < this.select.length; ai++) {
	    let fragment = this.select[4*ai]; // fragment is red channel
	    if (fragment > maxFragment) maxFragment = fragment;
	}
	//console.log('max fragment', maxFragment);
	
	// Get the current maximum monomer number
	let maxMonomer = 0;
	for (let ai = 0; ai < this.select.length; ai++) {
	    let monomer = this.select[4*ai+1]; // monomer is green channel
	    if (monomer > maxMonomer) maxMonomer = monomer;
	}

	// Fragment
	if (mol.atom_fragment != null) {
	    for (let ai = 0; ai < mol.num; ai++)
		this.select[4*map[ai]] = mol.atom_fragment[ai] + maxFragment + 1;
	} else {
	    // If there is no fragment array, the whole inserted structure will be a single fragment
	    for (let ai = 0; ai < mol.num; ai++)
		this.select[4*map[ai]] = maxFragment + 1;
	}

	// Monomer (residue number)
	if (mol.atom_monomer != null) {
	    for (let ai = 0; ai < mol.num; ai++)
		this.select[4*map[ai]+1] = mol.atom_monomer[ai] + maxMonomer + 1;
	} else {
	    // If there is no monomer array, the whole inserted structure will be a single monomer
	    for (let ai = 0; ai < mol.num; ai++)
		this.select[4*map[ai]+1] = maxMonomer + 1;
	}

	// Material
	if (mol.index_resname != null) {
	    for (let ai = 0; ai < mol.num; ai++) {
		const resName = mol.par_resname[mol.index_resname[ai]];
		let material = MaterialEnum.normal;
		if (this.waterResNames.includes(resName)) {
		    material = MaterialEnum.water;
		} else if (this.backgroundResNames.includes(resName)) {
		    material = MaterialEnum.background;
		}
		this.select[4*map[ai]+2] = material;
	    }
	}

	// Radius
	for (let ai = 0; ai < mol.num; ai++) {
	    // Just copy from the colorData (which must be filled first)
	    this.select[4*map[ai]+3] = this.color[4*map[ai]+3];
	}
    }

    // Stores [mass, charge, LJ_epsilon, LJ_Rmin] in the texture data array
    // The texture has the same size as other atom textures (one RGBA per atom)
    insertNonbond(mol, map) {
	// Get the masses
	const molMass = new Float32Array(mol.num);
	if (mol.index_mass == null) {
	    console.log('Warning! Inserting molecule ${mol.name}: Masses are not present. They will be guessed');
	    for (let ai = 0; ai < mol.num; ai++)
		molMass[ai] = this.guessMass(mol.par_type[mol.index_type[ai]]);
	} else {
	    for (let ai = 0; ai < mol.num; ai++)
		molMass[ai] = mol.par_mass[mol.index_mass[ai]];
	}

	// Get the charges
	const molCharge = new Float32Array(mol.num);
	if (mol.index_charge == null) {
	    console.log('Warning! Inserting molecule ${mol.name}: Charges are not present. They will be zero.');
	} else {
	    for (let ai = 0; ai < mol.num; ai++)
		molCharge[ai] = mol.par_charge[mol.index_charge[ai]];
	}

	// Set the mass, charge, LJ epsilon, and LJ RMin
	for (let i = 0; i < map.length; i++) {
	    const texAtomIndex = map[i];
	    const ljIndex = mol.index_lj[i];
	    const ljArrayIndex = TermInfo.parFields[TermInfo.lj]*ljIndex;
	    
	    this.nonbond[4*texAtomIndex] = molMass[i]; 
	    this.nonbond[4*texAtomIndex+1] = molCharge[i];
	    this.nonbond[4*texAtomIndex+2] = mol.par_lj[ljArrayIndex];
	    this.nonbond[4*texAtomIndex+3] = mol.par_lj[ljArrayIndex+1];
	}
    }

    // Insert bond, angle, dihedral, or exclusion terms into the appropriate texture
    insertTerms(mol, map, termKind) {
	const degToRad = Math.PI/180.0;
	
	// What arrays in the molecular structure do we need
	let molTerm = null;
	let molPar = null;
	let molIndex = null;
	let texData = null;
	if (termKind == TermInfo.bond) {
	    molTerm = mol.term_bond;
	    molPar = mol.par_bond;
	    molIndex = mol.index_bond;
	    texData = this.bond;
	} else if (termKind == TermInfo.angle) {
	    molTerm = mol.term_angle;
	    molPar = mol.par_angle;
	    molIndex = mol.index_angle;
	    texData = this.angle;
	} else if (termKind == TermInfo.dihedral) {
	    molTerm = mol.term_dihedral;
	    molPar = mol.par_dihedral;
	    molIndex = mol.index_dihedral;
	    texData = this.dihedral;
	} else if (termKind == TermInfo.exclude) {
	    molTerm = mol.term_exclusion;
	    molPar = mol.par_exclusion;
	    molIndex = mol.index_exclusion;
	    texData = this.exclude;
	} else {
	    console.log(`Warning! insertTerms() for molecule ${mol.name}: Unrecognized termKind ${termKind}. L-J terms have their own function insertNonbonded()`);
	    return;
	}

	// Set the size of the texture region
	const termFields = 4*TermInfo.pixels[termKind];
	const termsPerAtom = this.termsPerAtom[termKind];

	// Zero out the terms for this atom
	for (let ai = 0; ai < mol.num; ai++) {
	    const texAtomIndex = map[ai];
	    const termTexStart = texAtomIndex * termsPerAtom * termFields;
	    for (let j = 0; j < termsPerAtom * termFields; j++) {
		texData[termTexStart + j] = 0.0;
	    }
	}

	// If there are no terms of this kind defined, we are done
	if (molTerm == null || molIndex == null || molPar == null) return;

	// Consistency check for number of parameter indices and number of parameter values
	const termNum = molTerm.length/TermInfo.bodies[termKind];
	if (termNum != molIndex.length) {
	    const termName = TermInfo.name[termKind];
	    console.log(`ERROR! Number of terms (${termNum}) in molecule ${mol.name} term_${termName} is not equal to number of term parameter indices (${molIndex.length}) in mol.index_${termName}`);
	    return;
	}

	// Consistency check for maximum parameter index and number of parameter sets
	let parIndexMax = molIndex[0];
	const parSetNum = molPar.length/TermInfo.parFields[termKind];
	for (const parIndex of molIndex) {
	    if (parIndex > parIndexMax) parIndexMax = parIndex;
	    if (parIndex < 0) {
		console.log(`ERROR! molecule ${mol.name} index_${termName} has a negative index!`);
		return;
	    }
	}
	if (parIndexMax >= parSetNum) {
	    console.log(`ERROR! Maximum parameter index (${parIndexMax}) in molecule ${mol.name} index_${termName} is too large for the number of parameter sets (${parSetNum}) in par_${termName}`);
	    return;
	}

	// For each atom, there are entries in each texture specifying the index
	// of the other atoms (atomNeighList) and the FF parameters (atomParList)
	const atomNeighList = [];
	const atomParList = [];
	for (let ai = 0; ai < mol.num; ai++) {
	    atomNeighList[ai] = []; // neighbor list (atom indices) for each term
	    atomParList[ai] = []; // parameters for each term

	    // Find all terms involving the current atom
	    for (let ti = 0; ti < molTerm.length; ti++) {	    
		// This term contains our current atom
		if (molTerm[ti] == ai) {
		    const termIndex = Math.floor(ti/TermInfo.bodies[termKind]);
		    const termStart = termIndex*TermInfo.bodies[termKind];
		    const termEnd = termStart + TermInfo.bodies[termKind];

		    // Get the neighbors, including the current atom
		    let neighList = molTerm.slice(termStart, termEnd);
		    // Is our atom toward the beginning of the list or toward the end?
		    if (ti - termStart >= TermInfo.bodies[termKind]/2) {
			// If it is toward the end, reverse it so it's at the beginning
			neighList.reverse();
		    }
		    // Check if our atom is at the edge or internal
		    const role = (neighList[0] == ai)?0:1;
		    // Filter our current atom out of the neighList
		    let ourIndex = neighList.indexOf(ai);
		    neighList = neighList.slice(0,ourIndex).concat(neighList.slice(ourIndex+1));

		    // Get the parameters (the first parameter is the role)
		    const parStart = molIndex[termIndex]*TermInfo.parFields[termKind];
		    const parEnd = parStart + TermInfo.parFields[termKind];
		    const parList = [role, ...molPar.slice(parStart, parEnd)];

		    // Add the term to this atom
		    atomNeighList[ai].push(neighList);
		    atomParList[ai].push(parList);
		} // end if: term matches our current atom
	    } // end loop over terms
	} // end loop over atoms

	// Consistency check:
	// What is the maximum number of terms per atom?
	let maxAtomTerms = atomParList[0].length;
	for (let ai = 1; ai < mol.num; ai++) {
	    if (atomParList[ai].length > maxAtomTerms) maxAtomTerms = atomParList[ai].length;
	}
	if (maxAtomTerms > termsPerAtom) {
	    console.log(`ERROR! insertTerms() for molecule ${mol.name}: there are too many ${TermInfo.name[termKind]} terms (${maxAtomTerms}) and the texture can only hold ${termsPerAtom}`);
	    return;
	}

	// Remap the neighbor list to the indices used in the texture
	for (let ai = 0; ai < mol.num; ai++) {
	    for (let ti = 0; ti < atomNeighList[ai].length; ti++) {
		for (let c = 0; c < atomNeighList[ai][ti].length; c++) {
		    atomNeighList[ai][ti][c] = map[atomNeighList[ai][ti][c]];
		}
	    }
	}

	// Set the texture data for each atom
	for (let ai = 0; ai < mol.num; ai++) {
	    const texAtomIndex = map[ai];
	    const termTexStart = texAtomIndex * termsPerAtom * termFields;
	    for (let ti = 0; ti < atomParList[ai].length; ti++) {
		const j = termTexStart + ti * termFields;
		if (TermInfo.pixels[termKind] == 1) {
		    // This is for bonds and exclusions
		    texData[j] = atomNeighList[ai][ti][0]; // Index of other atom
		    texData[j+1] = atomParList[ai][ti][0]; // Role (always zero)
		    texData[j+2] = atomParList[ai][ti][1]; // bond constant or epsilon
		    texData[j+3] = atomParList[ai][ti][2]; // bond dist. or RMin
		} else {
		    // Angles and dihedrals
		    // Atom indices are in the first RGBA pixel and parameters are in the second
		    // The first parameter is the "role" (edge or internal atom)
		    for (let c = 0; c < atomNeighList[ai][ti].length && c<4; c++)
			texData[j+c] = atomNeighList[ai][ti][c];
		    
		    for (let c = 0; c < atomParList[ai][ti].length && c<4; c++) {
			let scale = 1.0;
			// Convert angle and dihedral terms from degrees to radians
			if ((termKind == TermInfo.angle || termKind == TermInfo.dihedral)
			    && c == atomParList[ai][ti].length - 1) scale = degToRad;

			// Set the second pixel
			texData[j+c+4] = scale*atomParList[ai][ti][c];
		    } // loop over components of the second pixel
		} // conditional for one pixel or two
		// We have completed setting one term
	    } // loop over terms applied to this atom
	} // loop over atoms in the molecule to be inserted

	return;
    }
    
    // mol is a molecular structure of the kind produced by charmm_to_javascript.js
    // termKind is one of TermInfo.bond, TermInfo.angle, TermInfo.dihedral, TermInfo.exclude
    static maxTermsPerAtom(mol, termKind) {
	// Get the term array in the molecular structure
	let molTerm = null;
	if (termKind == TermInfo.bond) {
	    molTerm = mol.term_bond;
	} else if (termKind == TermInfo.angle) {
	    molTerm = mol.term_angle;
	} else if (termKind == TermInfo.dihedral) {
	    molTerm = mol.term_dihedral;
	} else if (termKind == TermInfo.exclude) {
	    molTerm = mol.term_exclusion;
	} else {
	    console.log(`Warning! maxTermsPerAtom(): Unrecognized termKind ${termKind}.`);
	    return 1;
	}

	// There are no terms of this kind
	// The texture must have at least one term, so return 1
	if (molTerm == null) return 1;

	// The size of the term texture must be at least one
	let maxTerms = 1;
	// Loop over the atoms
	for (let ai = 0; ai < mol.num; ai++) {
	    // Count all terms involving the current atom
	    let termCount = 0;
	    for (let ti = 0; ti < molTerm.length; ti++) {	    
		// This term contains our current atom
		if (molTerm[ti] == ai) termCount++;
	    }
	    // Update the maximum
	    if (termCount > maxTerms) maxTerms = termCount;
	}

	return maxTerms;
    }

    sizeHashTable(hashTableFillFactor = 23) {
	// We use a large default fill factor so we only have 1 or 2 collisions at worst
	//hashTableFillFactor = 201; // try to get zero collisions (not worth it)
	
	// Count the number of active exclusions (or special L-J parameters) in the exclusion data
	const excludePerAtom = this.termsPerAtom[TermInfo.exclude];
	const excludeList = [];
	for (let ai = 0; ai < this.size; ai++) {
	    // Loop through the exclusion terms
	    const atomStart = 4*excludePerAtom*ai; // starting index in exclusion data
	    for (let ti = 0; ti < excludePerAtom; ti++) {
		const termStart = atomStart + 4*ti;
		// Active exclusions have non-zero RMin
		if (this.exclude[termStart+3] != 0.0 && ai < this.exclude[termStart]) {
		    excludeList.push([this.exclude[termStart], this.exclude[termStart+1]]);
		}
	    }
	}
	const excludeCount = excludeList.length;

	const primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991];

	// Set the hash table dimensions
	const hashSize = hashTableFillFactor*excludeCount;
	const hashSqrt = Math.ceil(Math.sqrt(hashSize));
	if (hashSqrt > primes[primes.length-2]) {
	    // We can't use values bigger than 16777216
	    // mod() begins to give problems for too big of integers
	    this.hashDim = {width: 983, height: 991};
	} else {
	    // Find the first two primes greater than or equal to excludeSize
	    for (let i = 0; i < primes.length-1; i++) {
		if (primes[i] >= hashSqrt) {
		    this.hashDim = {width: primes[i], height: primes[i+1]};
		    break;
		}
	    }
	}
	this.hashDim.size = this.hashDim.width * this.hashDim.height;
	console.log('exclusion hash table size:', this.hashDim.width, this.hashDim.height);
	console.log(`There are ${excludeCount} active exclusions`);
	const hashLoadPercent = 100.0*excludeCount/this.hashDim.size;
	console.log(`Hash table load: ${hashLoadPercent.toFixed(2)}%`);
    }

    // Create a new hash table with optimized parameters
    makeHashTable(width, height) {
	// Optimize the hash table parameters for the least collisions
	const hashTableAttempts = 15;
	this.stopCollisions = 20; // Don't allow more than some number of collisions

	// Create a new hash table size
	if (width == null || height == null) {
	    this.sizeHashTable();
	}
	
	// Note that "primes" (used for the hash table texture dimensions) and "bigPrimes" (used for the multipliers)
	// must be disjoint. (a == m is bad for the linear congruential generator).
	// Also, 991*16381 < 2^24 (needed for exact handling of integers with 32-bit floats)
	const bigPrimes = [15013, 15017, 15031, 15053, 15061, 15073, 15077, 15083, 15091, 15101, 15107, 15121, 15131, 15137, 15139, 15149, 15161, 15173, 15187, 15193, 15199, 15217, 15227, 15233, 15241, 15259, 15263, 15269, 15271, 15277, 15287, 15289, 15299, 15307, 15313, 15319, 15329, 15331, 15349, 15359, 15361, 15373, 15377, 15383, 15391, 15401, 15413, 15427, 15439, 15443, 15451, 15461, 15467, 15473, 15493, 15497, 15511, 15527, 15541, 15551, 15559, 15569, 15581, 15583, 15601, 15607, 15619, 15629, 15641, 15643, 15647, 15649, 15661, 15667, 15671, 15679, 15683, 15727, 15731, 15733, 15737, 15739, 15749, 15761, 15767, 15773, 15787, 15791, 15797, 15803, 15809, 15817, 15823, 15859, 15877, 15881, 15887, 15889, 15901, 15907, 15913, 15919, 15923, 15937, 15959, 15971, 15973, 15991, 16001, 16007, 16033, 16057, 16061, 16063, 16067, 16069, 16073, 16087, 16091, 16097, 16103, 16111, 16127, 16139, 16141, 16183, 16187, 16189, 16193, 16217, 16223, 16229, 16231, 16249, 16253, 16267, 16273, 16301, 16319, 16333, 16339, 16349, 16361, 16363, 16369, 16381, 16411, 16417, 16421, 16427, 16433, 16447, 16451, 16453, 16477, 16481, 16487, 16493, 16519, 16529, 16547, 16553, 16561, 16567, 16573, 16603, 16607, 16619, 16631, 16633, 16649, 16651, 16657, 16661, 16673, 16691, 16693, 16699, 16703, 16729, 16741, 16747, 16759, 16763, 16787, 16811, 16823, 16829, 16831];

	// The hash table will store the exclusions in a (atomIndex0 atomIndex1 epsilon Rmin)
	// With atomIndex0 < atomIndex1
	this.hashA = new Float32Array(2);
	let bestA = new Float32Array(2);
	let bestCollisions = this.stopCollisions;
	let worstCollisions = 0;
	let sumCollisions = 0;

	// Optimize the hash table parameters for the least collisionsp
	for (let i = 0; i < hashTableAttempts; i++) {
	    // Select random hash function factors
	    const pi0 = Math.floor(bigPrimes.length*Math.random());
	    const pi1 = Math.floor(bigPrimes.length*Math.random());
	    this.hashA[0] = bigPrimes[pi0];
	    this.hashA[1] = bigPrimes[pi1];
	    this.hashAtoms();

	    sumCollisions += this.maxCollisions;
	    if (this.maxCollisions < bestCollisions) {
		bestCollisions = this.maxCollisions;
		bestA[0] = this.hashA[0];
		bestA[1] = this.hashA[1];
	    }
	    if (this.maxCollisions > worstCollisions)
		worstCollisions = this.maxCollisions;
	    //console.log(`hashA ${this.hashA[0]} ${this.hashA[1]} col ${this.maxCollisions}`);
	}
	const meanMaxCollisions = sumCollisions/hashTableAttempts;

	// We've finished. Set to the best value
	this.hashA[0] = bestA[0];
	this.hashA[1] = bestA[1];
	// Regenerate the hash table using the best multipliers
	this.hashAtoms();
	// The shaders are SIMD, so we need to hardcode the maximum number
	// of collisions during any search
	
	if (this.maxCollisions >= this.stopCollisions) {
	    console.log(`ERROR! Reached maximum number of hash table collisions ${this.stopCollisions}. Hash table is incomplete. It is probably overloaded.`);	    
	}
	console.log('Mean maximum collisions in attempted exclusion hash tables:', meanMaxCollisions);
	console.log('Worst maximum collisions in attempted exclusion hash tables:', worstCollisions);
	console.log('Maximum collisions in optimized exclusion hash table:', this.maxCollisions);
	console.log('Multipliers for optimized exclusion hash table:', this.hashA[0], this.hashA[1]);
    }

    // This function creates the hash table from scratch
    hashAtoms() {
	// Create the hash data array
	this.hash = new Float32Array(4*this.hashDim.size);
	
	// Loop through the atoms
	const excludePerAtom = this.termsPerAtom[TermInfo.exclude];
	this.maxCollisions = 0;
	for (let ai = 0; ai < this.size; ai++) {
	    // Loop through the exclusion terms
	    const atomStart = 4*excludePerAtom*ai; // starting index in exclusion data
	    for (let ti = 0; ti < excludePerAtom; ti++) {
		const termStart = atomStart + 4*ti;
		// Each exclusion occupies one RGBA pixel
		const atomIndex1 = this.exclude[termStart];
		const epsilon = this.exclude[termStart+2];
		const RMin = this.exclude[termStart+3];

		// Hash only when atomIndex0 < atomIndex1 to not duplicate
		// Hash only active exclusions
		if (ai < atomIndex1 && RMin != 0) {
		    const hashPar = [ai, atomIndex1, epsilon, RMin];
		    const collisions = this.hashInsert(hashPar);
		    //console.log('exclusions', ai, atomIndex1, epsilon, RMin, 'col', collisions);
		}
	    } // End of exclusion term loop
	} // End of atom loop
	return this.maxCollisions;
    }

    // This must be done after the exclusion terms are inserted (this.exclude is valid)
    insertMoleculeIntoHashTable(map) {
	const excludePerAtom = this.termsPerAtom[TermInfo.exclude];
	for (let ai = 0; ai < map.length; ai++) {
	    const atomIndex0 = map[ai];
	    
	    // Read through the exclusions associated with this atom
	    const atomStart =  4*excludePerAtom*atomIndex0;
	    for (let ti = 0; ti < excludePerAtom; ti++) {
		const termStart = atomStart + 4*ti;
		// Each exclusion occupies one RGBA pixel
		const atomIndex1 = this.exclude[termStart];
		const epsilon = this.exclude[termStart+2];
		const RMin = this.exclude[termStart+3];

		// Hash only when atomIndex0 < atomIndex1 to not duplicate
		// Hash only active exclusions
		if (atomIndex0 < atomIndex1 && RMin != 0) {
		    const hashPar = [atomIndex0, atomIndex1, epsilon, RMin];
		    this.hashInsert(hashPar);
		}
	    } // End of exclusion term loop
	} // End of atom loop
	return this.maxCollisions;
    }

    // The atom indices must first be conditioned for the linear congruential generator
    startHash(atomIndex0, atomIndex1) {
	// Make sure we have good starting numbers (x % this.hashDim.width == 0 is bad)
	const i0 = (atomIndex0 % (this.hashDim.width-1)) + 1;
	const i1 = (atomIndex1 % (this.hashDim.height-1)) + 1;
	return this.nextHash(i0, i1);
    }

    // Hashing function based on a linear congruential generator
    nextHash(i0, i1) {
	// Use a linear congruential generator get new hash table coordinates
	const ret0 = (this.hashA[0] * i0) % this.hashDim.width;
	const ret1 = (this.hashA[1] * i1) % this.hashDim.height;
	return [ret0, ret1];
    }

    // Get exclusion or special LJ parameters from the hash table
    hashRead(atomIndex0, atomIndex1) {
	// Hash the atom indices
	let [i0, i1] = this.startHash(atomIndex0, atomIndex1);
	let hashStart = 4*(i0 + i1*this.hashDim.width);
	let hashPar = this.hash.slice(hashStart,hashStart+4);
	
	// Check the first position in the table
	console.log(`Looked for ${atomIndex0} ${atomIndex1} at ${i0} ${i1}`);
	if (hashPar[0] == atomIndex0 && hashPar[1] == atomIndex1) {
	    console.log(`Found ${atomIndex0} ${atomIndex1} at ${i0} ${i1}`);
	    // Found!
	    return hashPar;
	}

	// Hash again
	for (let collide = 1; collide <= this.maxCollisions; collide++) {
	    // Hash the last hash function result
	    [i0, i1] = this.nextHash(i0, i1);
	    console.log('i0 i1', i0, i1);
	    hashStart = 4*(i0 + i1*this.hashDim.width);
	    hashPar = this.hash.slice(hashStart,hashStart+4);
	    console.log(`Looked for ${atomIndex0} ${atomIndex1} at ${i0} ${i1}  after ${collide} collisions`);
	    if (hashPar[0] == atomIndex0 && hashPar[1] == atomIndex1) {
		console.log(`Found ${atomIndex0} ${atomIndex1} at ${i0} ${i1} after ${collide} collisions`);
		return hashPar;
	    }
	}

	// Not found
	// Return null
	return null;
    }
    
    // Insert exclusion or special L-J into the hash table
    // Returns the number of collisions
    // insertPar has the form [atomIndex0, atomIndex1 epsilon RMin]
    hashInsert(insertPar) {
	const atomIndex0 = insertPar[0];
	const atomIndex1 = insertPar[1];
	
	// Hash the atom indices
	let [i0, i1] = this.startHash(atomIndex0, atomIndex1);
	let hashStart = 4*(i0 + i1*this.hashDim.width);
	let hashPar = this.hash.slice(hashStart,hashStart+4);
	
	let collisions = 0;
	// Active exclusions have RMin != 0
	while (hashPar[3] != 0.0 && collisions < this.stopCollisions) {
	    // We've had a collision
	    collisions++; 
	    // Hash the last hash function result
	    [i0, i1] = this.nextHash(i0, i1);
	    hashStart = 4*(i0 + i1*this.hashDim.width);
	    hashPar = this.hash.slice(hashStart,hashStart+4);
	}

	// This is an empty space
	//console.log(`inserting ${atomIndex0} ${atomIndex1} at ${i0} ${i1} with ${collisions} collisions`);
	this.hash[hashStart] = insertPar[0];
	this.hash[hashStart+1] = insertPar[1];
	this.hash[hashStart+2] = insertPar[2];
	this.hash[hashStart+3] = insertPar[3];

	// Increase maxCollisions if necessary
	if (collisions > this.maxCollisions) this.maxCollisions = collisions;
	return collisions;
    }
    
    
    // Create a PDB string 
    writePDB() {
	const record = 'ATOM  ';
	const segPrefix = 'F';
	const alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
	const temp1 = '    ';
	const temp2 = '  1.00';
	const beta = '  0.00'
	const temp3 = '      ';
	const tempEnd = ' ';

	// Make a list of active atoms
	const atoms = [];
	for (let ai = 0; ai < this.size; ai++) {
	    if (this.pos[4*ai+3] != 0) {
		atoms.push(ai);
	    }
	}
	// Sort by monomer number (select[4*i+1])
	atoms.sort((ai0, ai1) => this.select[4*ai0+1] - this.select[4*ai1+1]);
	
	let PDBString = '';
	let serial = 0;
	let resId = 0;
	for (const i of atoms) {
	    // Ignore inactive atoms
	    if (this.pos[4*i+3] == 0) continue;
	    serial++;

	    // Get properties
	    const type = this.typeList[i];
	    const name = this.nameList[i];
	    const resName = this.resNameList[i];
	    const fragment = this.select[4*i];
	    const monomer = this.select[4*i+1];

	    // Position
	    const r = [];
	    for (let c = 0; c < 3; c++) {
		r[c] = this.pos[4*i+c].toFixed(3);
	    }

	    // Fields
	    const si = `     ${serial} `.slice(-6);
	    let name1 = null;
	    if (name.length < 4) {
		name1 = ` ${name}    `.slice(0,4);
	    } else {
		name1 = name.slice(0,4);
	    }
	    const chain = alphabet[fragment % alphabet.length];
	    const resId = `    ${monomer}`.slice(-4,);
	    const resName1 = ` ${resName}     `.slice(0,5);
	    const sx = `       ${r[0]}`.slice(-8);
	    const sy = `       ${r[1]}`.slice(-8);
	    const sz = `       ${r[2]}`.slice(-8);
	    const segName = `${segPrefix}${fragment}    `.slice(0,4);
	    const element = `${type[0]} `.slice(0,1);
	    PDBString += `${record}${si}${name1}${resName1}${chain}${resId}${temp1}${sx}${sy}${sz}${temp2}${beta}${temp3}${segName}${tempEnd}${element}\n`;
	}
	return PDBString;
    }

    static getColor(type) {
	const color = [];
	color.H =  [0.80, 0.80, 0.80];
	color.B =  [1.00, 0.60, 0.60];
	color.C =  [0.20, 0.62, 0.32];
	color.N =  [0.20, 0.20, 0.80];
	color.O =  [0.85, 0.15, 0.15];
	color.F =  [0.75, 1.00, 0.00];
	//color.Na = [0.84, 0.67, 0.59];
	color.Na = [0.95, 0.50, 0.09];
	color.Si = [0.20, 0.57, 0.60];
	color.P =  [0.78, 0.33, 0.09];
	color.S =  [0.77, 0.77, 0.24];
	color.Se = [0.39, 0.17, 0.05];
	color.Cl = [0.00, 1.00, 0.00];
	color.K  = [0.59, 0.67, 0.84];
	color.Br = [0.54, 0.00, 0.12];
	color.I =  [0.65, 0.00, 0.65];
	//color.M =  [0.55, 0.55, 0.55];
	color.M =  [0.88, 0.10, 0.88];
	color.Ne =  [0.8, 0.4, 0.8];
	color.unknown = [0.2, 0.2, 0.2];
	
	if (type.match(/^CLG/) || type == 'CLA') {
	    // CGenFF chlorine or CHARMM chloride ion
	    return color.Cl;
	} else if (type.match(/^SEG/) || type.match(/^Se/)) {
	    // CGenFF selenium
	    return color.Se;
	} else if (type.match(/^BR/) || type.match(/^Br/)) {
	    // Probably bromine and not boron
	    return color.Br;
	} else if (type.match(/^I/)) {
	    // Could catch indium
	    return color.I;
            // CHARMM ion type names
	} else if (type == 'SOD' || type.match(/^Na/)) {
	    return color.Na;
	} else if (type == 'POT' || type.match(/^K/)) {
	    return color.K;
	} else if (type == 'LIT' || type.match(/^Li/)) {
	    return color.M;
	} else if (type == 'MG' || type.match(/^Mg/)) {
	    return color.M;
	} else if (type == 'RUB' || type.match(/^Rb/)) {
	    return color.M;
	} else if (type == 'CES' || type.match(/^Cs/)) {
	    return color.M;
	} else if (type == 'BAR' || type.match(/^Ba/)) {
	    return color.M;
	} else if (type == 'ZN' || type.match(/^Zn/)) {
	    return color.M;
	} else if (type == 'CAD' || type.match(/^Cd/)) {
	    return color.M;
	} else if (type == 'FE' || type.match(/^Fe/)) {
	    return color.M;
	} else if (type == 'Ni') {
	    return color.M;
	} else if (type == 'CAL') {
	    return color.M;
	} else if (type == 'HE') {
	    // probably won't catch hydrogen
	    return color.Ne;
	} else if (type.match(/^SI/) || type.match(/^Si/)) {
	    // This could catch some sulfurs
	    return color.Si;
	} else if (type == 'NE' || type.match(/^Ne/)) {
	    // probably won't catch nitrogen
	    return color.Ne;
	} else if (type.match(/^A/)) {
	    return color.Ne;
	} else if (type.match(/^H/)) {
	    return color.H
	} else if (type.match(/^B/)) {
	    return color.B;
	} else if (type.match(/^C/)) {
	    return color.C;
	} else if (type.match(/^N/)) {
	    return color.N;
	} else if (type.match(/^O/)) {
	    return color.O;
	} else if (type.match(/^F/)) {
	    return color.F;
	} else if (type.match(/^P/)) {
	    return color.P;
	} else if (type.match(/^S/)) {
	    return color.S;
	} else {
	    return color.unknown;
	}
    }

    static guessMass(type) {
	if (type.match(/^CLG/) || type == 'CLA') {
	    // CGenFF chlorine or CHARMM chloride ion
	    return 35.453;
	} else if (type.match(/^BR/) || type.match(/^Br/)) {
	    // Probably bromine and not boron
	    return 79.904;
	} else if (type.match(/^SEG/) || type.match(/^Se/)) {
	    // CGenFF selenium
	    return 78.96;
	} else if (type.match(/^I/)) {
	    // Could catch indium
	    return 126.9045;
	    // CHARMM ion type names
	} else if (type == 'SOD' || type.match(/^Na/)) {
	    return 22.9897;
	} else if (type == 'POT' || type.match(/^K/)) {
	    return 39.0983;
	} else if (type == 'LIT' || type.match(/^Li/)) {
	    return 6.941;
	} else if (type == 'MG' || type.match(/^Mg/)) {
	    return 24.305;
	} else if (type == 'RUB' || type.match(/^Rb/)) {
	    return 85.4678;
	} else if (type == 'CES' || type.match(/^Cs/)) {
	    return 132.9055;
	} else if (type == 'BAR' || type.match(/^Ba/)) {
	    return 137.327;
	} else if (type == 'ZN' || type.match(/^Zn/)) {
	    return 65.39;
	} else if (type == 'CAD' || type.match(/^Cd/)) {
	    return 112.411;
	} else if (type == 'FE' || type.match(/^Fe/) || type == 'M') {
	    return 55.845;
	} else if (type == 'Ni') {
	    return 58.6934;
	} else if (type == 'CAL') {
	    return 40.078;
	} else if (type == 'HE' || type.match(/^He/)) {
	    return 1.0079;
	} else if (type == 'NE' || type.match(/^Ne/)) {
	    return 20.1797;
	} else if (type.match(/^A/)) {
	    return 39.948;
	} else if (type.match(/^Si/) || type.match(/^SI/)) {
	    return 28.0855;
	} else if (type.match(/^H/)) {
	    return 1.0079;
	} else if (type.match(/^B/)) {
	    return 10.811;
	} else if (type.match(/^C/)) {
	    return 12.0107;
	} else if (type.match(/^N/)) {
	    return 14.0067;
	} else if (type.match(/^O/)) {
	    return 15.9994;
	} else if (type.match(/^F/)) {
	    return 18.9984;
	} else if (type.match(/^P/)) {
	    return 30.9738;
	} else if (type.match(/^S/)) {
	    return 32.065;
	} else {
	    return 12.0;
	}
    }
    
}


///////////////////////////////////////////////////////////////////////////////
// A class for the molecular system
class SimulationSystem {
    constructor(coordTexData, nonbondTexData, box, timestep) {
	this.timestep = timestep; // in femtoseconds
	this.step = 0; // current time step
	this.temper = 310.0; // in kelvin
	this.box = box;
	this.boltz = 0.001987191; // k_Boltzmann in kcal/mol/K (NAMD value)
	this.coulomb = 332.0636; // Coulomb constant in kcal/mol*angstrom/e^2 (NAMD value)
	this.timeFactor = 48.88821; // units 'sqrt(dalton*angstrom^2/(kcal_mol))' 'fs' (NAMD value)
	this.columns = 4; // Textures are RGBA
	this.langevinDamping = 0.5; // in ps^-1
	this.radius = 4.0; // Roughly the size of an atom

	// Derived values
	this.dt = this.timestep/this.timeFactor; // timestep in sim. units
	this.dt2 = this.dt*this.dt;
	this.kT = this.boltz*this.temper; // in kcal/mol
	this.origin = [-0.5*this.box[0], -0.5*this.box[1], -0.5*this.box[2]];
	this.destination = [-this.origin[0], -this.origin[1], -this.origin[2]];
	this.wall = [this.destination[0]-this.radius, this.destination[1]-this.radius, this.destination[2]-this.radius];
	this.wallSpring = 2.0;
	this.restraintSpring = 10.0;
	// Correctly set all the values derived from langevinDamping
	this.setLangevin(this.temper, this.langevinDamping);
	
	// Have multiple frames of positions, velocities, and forces
	this.fields = coordTexData.length;
	this.capacity = Math.floor(coordTexData.length/this.columns);
	this.id = [...Array(this.capacity).keys()];
	this.pos = [];
	this.vel = [];
	for (let t = 0; t <= 1; t++) {
	    this.pos[t] = new Float32Array(this.fields);
	    this.vel[t] = new Float32Array(this.fields);
	}

	// Set the masses
	this.mass = new Float32Array(this.capacity);
	for (let i = 0; i < nonbondTexData.length; i+=this.columns) {
	    const atomIndex = Math.floor(i/this.columns);
	    this.mass[atomIndex] = nonbondTexData[i];
	}

	// Copy coordinates into pos[0]
	for (let i = 0; i < coordTexData.length; i++) {
	    this.pos[0][i] = coordTexData[i];
	}

	// Set the velocities and previous positions
	this.initVelocities(this.temper);
    }

    setLangevin(langevinTemper, langevinDamping) {
	this.temper = langevinTemper; // in kelvin
	this.kT = this.boltz*this.temper; // in kcal/mol
	this.langevinDamping = langevinDamping; // in ps^-1
	// langevinGamma is the langevinDamping in simulation units, sqrt(kcal_mol/(dalton*Å^2))
	this.langevinGamma = this.langevinDamping/1000.0*this.timeFactor;
	// Unitness time factor for the BBK integrator
	this.langevinRatio = 0.5*this.langevinGamma*this.dt; 
	// This is the multiplier for the random force
	// Units of kcal_mol/simTime^2 = kcal_mol^2/(dalton*Å^2) 
	// This ends up getting multiplied by dt or dt^2 in the integrator,
	// so the conditional avoids spurious infinities
	this.langevinRandom = (this.dt==0.0) ? 0.0 : (2.0*this.langevinGamma*this.kT/this.dt);

	// Set the velocities and previous positions
	this.initVelocities(this.temper);
    }
    
    getPos(atomIndex, t = 0) {
	const j = atomIndex*this.columns;
	return [this.pos[t][j], this.pos[t][j+1], this.pos[t][j+2]];
    }
    setPos(atomIndex, r, t = 0) {
	const j = atomIndex*this.columns;
	this.pos[t][j] = r[0];
	this.pos[t][j+1] = r[1];
	this.pos[t][j+2] = r[2];
	// Just copy the remaining columns
	// There is usually a fourth column that says whether an atom is active
	for (let c = 3; c < this.columns; c++)
	    this.pos[t][j+c] = this.pos[0][j+c];
    }

    getVel(atomIndex, t = 0) {
	const j = atomIndex*this.columns;
	return [this.vel[t][j], this.vel[t][j+1], this.vel[t][j+2]];
    }
    setVel(atomIndex, r, t = 0) {
	const j = atomIndex*this.columns;
	this.vel[t][j] = r[0];
	this.vel[t][j+1] = r[1];
	this.vel[t][j+2] = r[2];
    }

    initVelocities(temperature) {
	// Temperature is assumed to be in kelvin
	for (let i = 0; i < this.capacity; i++) {
	    const j = i*this.columns;
	    // Ignore inactive atoms
	    if (this.pos[0][j+3] == 0.0) continue;
	    
	    const vel0 = Math.sqrt(this.boltz*temperature/this.mass[i]);
	    const v = [vel0*this.randomNormal(), vel0*this.randomNormal(), vel0*this.randomNormal()];

	    const r0 = [this.pos[0][j], this.pos[0][j+1], this.pos[0][j+2]];
	    // Setting previous positions allows for either standard Verlet or velocity Verlet
	    let r1 = new Array(3);
	    for (let c = 0; c < 3; c++) {
		// This allows for starting standard Verlet, but the estimation of the
		// velocity is quite poor since the force isn't included
		r1[c] = r0[c] - v[c]*this.dt;
	    }
	    this.setVel(i, v, 0);
	    this.setVel(i, v, 1);
	    this.setPos(i, r1, 1);
	}
    }

    randomNormal() {
	return Math.sqrt(-2.0*Math.log(1.0 - Math.random())) * Math.cos(2.0*Math.PI*Math.random());
    }
    
    wrap(d) {
	return [this.wrapReal(d[0], this.box[0]), this.wrapReal(d[1], this.box[1]), this.wrapReal(d[2], this.box[2])];
    }

    // Wrap so the the result is -L/2 <= x < L/2
    wrapReal(x,len) {
	return x - len*Math.floor(x/len + 0.5);
    }
}



///////////////////////////////////////////////////////////////////////////////
// Handle rendering of the 2D thumbnails
class MoleculeThumbnail {
    constructor(mol, overrideName) {
	// Make an array with the positions, color, and radius of the atoms
	this.mol = mol;
	this.atomArray = [];
	this.screenScale = 0.95;
	this.name =  overrideName || mol.name;

	// Get the atom positions
	for (let ai = 0; ai < mol.index_type.length; ai++) {
	    const atom = {};
	    const type = mol.par_type[mol.index_type[ai]];
	    const ljRadiusIndex = TermInfo.parFields[TermInfo.lj]*mol.index_lj[ai] + 1;
	    const RMin = (mol.par_lj[ljRadiusIndex] < 2.4) ? 2.4 : mol.par_lj[ljRadiusIndex];
	    const color = AtomTextureData.getColor(type);

	    // Here we flip the y-axis because the L-amino acids looked like D
	    atom.pos = [mol.atom_coord[3*ai], -mol.atom_coord[3*ai+1], mol.atom_coord[3*ai+2]];
	    atom.color = color.map(v => Math.floor(255*v)); // [0, 1] -> [0, 255]
	    // LJ RMin is the equilibrium distance between two atoms. The radius is half that.
	    // But we also don't want the full radius
	    atom.radius = 0.5*RMin;
	    this.atomArray.push(atom);
	}

	// Get the bonds
	this.bondList = [];
	if (mol.term_bond != null) {
	    for (let i = 0; i < mol.term_bond.length; i+=2) {
		let ai0 = mol.term_bond[i];
		let ai1 = mol.term_bond[i+1];

		// Don't include bonds that are very long
		const maxBondLength = 6.0; // hopefully in Å
		// Calculate the distance
		let sumSq = 0.0;
		for (let c = 0; c < 3; c++) {
		    sumSq += (this.atomArray[ai0].pos[c] - this.atomArray[ai1].pos[c])**2;
		}
		const dist = Math.sqrt(sumSq);
		

		if (dist < maxBondLength) {
		    // Push the bond terms
		    // These my include Urey-Bradley terms, which we will remove.
		    const dx = this.atomArray[ai0].pos[2] + this.atomArray[ai1].pos[2];
		    const bondCenterZ = 0.5*(this.atomArray[ai0].pos[2] + this.atomArray[ai1].pos[2]);
		    this.bondList.push([ai0, ai1, bondCenterZ]);
		}
	    }
	}

	// Remove Urey-Bradley terms from the bond list
	if (mol.term_angle != null) {
	    const chemBondList = [];
	    for (const bond of this.bondList) {
		// Check if the atoms in this bond are 1,3 atoms of an angle
		let isUreyBradley = false;
		for (let i = 0; i < mol.term_angle.length; i+=3) {
		    let ai0 = mol.term_angle[i];
		    let ai2 = mol.term_angle[i+2];

		    if ((ai0 == bond[0] && ai2 == bond[1]) || (ai0 == bond[1] && ai2 == bond[0])) {
			// This is not a chemical bond, but instead a Urey-Bradley term
			isUreyBradley = true;
			break;
		    }
		}

		if (!isUreyBradley) {
		    chemBondList.push(bond);
		}
	    }
	    this.bondList = chemBondList;
	}
	    
	// Sort by the z value of the bond's center
	this.bondList.sort((a, b) => a[2] - b[2]);
	
	// Get the minimum and maximum coordinates in x and y
	let minX = this.atomArray[0].pos[0];
	let maxX = this.atomArray[0].pos[0];
	let minY = this.atomArray[0].pos[1];
	let maxY = this.atomArray[0].pos[1];
	for (const atom of this.atomArray) {
	    minX = Math.min(minX, atom.pos[0] - 0.5*atom.radius);
	    minY = Math.min(minY, atom.pos[1] - 0.5*atom.radius);
	    maxX = Math.max(maxX, atom.pos[0] + 0.5*atom.radius);
	    maxY = Math.max(maxY, atom.pos[1] + 0.5*atom.radius);
	}

	// Get the sizes in x and y and the maximum dimension
	const sizeX = maxX - minX;
	const sizeY = maxY - minY;
	let maxSize = Math.max(sizeX, sizeY);
	// This will happen for a single-atom molecule
	if (maxSize < 1e-6) {
	    maxSize = atom.radius;
	}
	
	// Scale and shift the molecules so the positions go from -0.5 to 0.5
	const cenX = 0.5*(minX + maxX);
	const cenY = 0.5*(minY + maxY);
	for (let i = 0; i < this.atomArray.length; i++) {
	    this.atomArray[i].pos[0] = (this.atomArray[i].pos[0] - cenX)/maxSize;
	    this.atomArray[i].pos[1] = (this.atomArray[i].pos[1] - cenY)/maxSize;
	    this.atomArray[i].radius /= maxSize;
	}
	//console.log('scale', minX, minY, maxX, maxY);
	//console.log(this.name, ':', this.atomArray);
    }

    draw(context2d) {
	//context2d.clearRect(0, 0, context2d.canvas.width, context2d.canvas.height);
	context2d.fillStyle = 'black';
	context2d.fillRect(0, 0, context2d.canvas.width, context2d.canvas.height);

	if (this.bondList.length >= 1) {
	    for (const bond of this.bondList) {
		this.drawBond(context2d, bond);
	    }
	} else {
	    for (const atom of this.atomArray) {
		this.drawSphere(context2d, atom, 0.4);
	    }
	}
    }

    drawBond(ctx, bond, bondScale = 0.1) {
	// Scale the positions
	const pix = this.screenScale*Math.min(ctx.canvas.width, ctx.canvas.height);
	const atom0 = this.atomArray[bond[0]];
	const atom1 = this.atomArray[bond[1]];
	const x0 = pix*atom0.pos[0] + 0.5*ctx.canvas.width;
	const y0 = pix*atom0.pos[1] + 0.5*ctx.canvas.height;
	const x1 = pix*atom1.pos[0] + 0.5*ctx.canvas.width;
	const y1 = pix*atom1.pos[1] + 0.5*ctx.canvas.height;
	const xm = 0.5*(x0 + x1);
	const ym = 0.5*(y0 + y1);

	// radii
	const rad0 = pix*atom0.radius*bondScale;
	const rad1 = pix*atom1.radius*bondScale;
	const radM = 0.5*(rad0 + rad1);
	
	// Get the atoms' colors
	const color0 = `rgb(${atom0.color[0]}, ${atom0.color[1]}, ${atom0.color[2]})`;
	const color1 = `rgb(${atom1.color[0]}, ${atom1.color[1]}, ${atom1.color[2]})`;

	// Determine the points
	const dx = x1 - x0;
	const dy = y1 - y0;
	const dist = Math.sqrt(dx*dx + dy*dy);
	// The vector along the bond
	const ex = dx/dist;
	const ey = dy/dist;
	// The orthogonal direction
	const nx = -ey;
	const ny = ex;

	// First segment
	{
	    const p0 = [x0 + rad0*nx, y0 + rad0*ny];
	    const p1 = [xm + radM*nx, ym + radM*ny];
	    const p2 = [xm - radM*nx, ym - radM*ny];
	    const p3 = [x0 - rad0*nx, y0 - rad0*ny];
	    ctx.fillStyle = color0;
	    ctx.strokeStyle = color0;
	    ctx.beginPath(); 
	    ctx.moveTo(...p0); 
	    ctx.lineTo(...p1);
	    ctx.lineTo(...p2);
	    ctx.lineTo(...p3);
	    ctx.lineTo(...p0);
	    ctx.fill();
	    ctx.stroke();
	}

	// Second segment
	{
	    const p0 = [x1 + rad1*nx, y1 + rad1*ny];
	    const p1 = [xm + radM*nx, ym + radM*ny];
	    const p2 = [xm - radM*nx, ym - radM*ny];
	    const p3 = [x1 - rad1*nx, y1 - rad1*ny];
	    ctx.fillStyle = color1;
	    ctx.strokeStyle = color1;
	    ctx.beginPath(); 
	    ctx.moveTo(...p0); 
	    ctx.lineTo(...p1);
	    ctx.lineTo(...p2);
	    ctx.lineTo(...p3);
	    ctx.lineTo(...p0);
	    ctx.fill();
	    ctx.stroke();
	}

	this.drawSphere(ctx, atom0, 0.25);
	this.drawSphere(ctx, atom1, 0.25);
    }
    
    drawSphere(ctx, atom, radiusScale = 0.6) {
	// Positions and color
	const pix = this.screenScale*Math.min(ctx.canvas.width, ctx.canvas.height);
	const x = pix*atom.pos[0] + 0.5*ctx.canvas.width;
	const y = pix*atom.pos[1] + 0.5*ctx.canvas.height;
	const rad = pix*atom.radius*radiusScale;
	const highlight = atom.color.map( v => (4*(v+50) > 255) ? 255 : 4*(v+50) );
	const shadow = atom.color.map( v => Math.floor(0.65*v) );
	const stop0 = `rgb(${highlight[0]}, ${highlight[1]}, ${highlight[2]})`;
	//const stop1 = `rgb(${atom.color[0]}, ${atom.color[1]}, ${atom.color[2]})`;
	const stop1 = `rgb(${shadow[0]}, ${shadow[1]}, ${shadow[2]})`;
	//console.log('draw', x, y, rad, stop0, stop1);
	
	// Positions for color gradient
	const grX0 = x - 0.3*rad;
	const grY0 = y - 0.3*rad;
	const grX1 = x - 0.3*rad;
	const grY1 = y - 0.3*rad;

	// Create the gradient
	//ctx.fillStyle = '';
	const gr = ctx.createRadialGradient(grX0,grY0,0.1*rad,grX1,grY1,0.95*rad);
	gr.addColorStop(0,stop0);
	gr.addColorStop(1,stop1);
	ctx.fillStyle = gr;

	// Draw the 
	ctx.beginPath();
	ctx.arc(x,y,rad, 0, Math.PI*2, true);
	ctx.closePath();
	ctx.fill();
	// ctx.fillStyle = "#000000";
        // ctx.font = "12px Arial";
	// ctx.fillText(body.id.toString(),x,y);
    };
    
}
