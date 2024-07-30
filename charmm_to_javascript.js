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

// charmm_to_javascript.js
// Classes for reading CHARMM psf, pdb, parameter, and stream files
'use strict';

let fs = require('fs');
let path = require('path');
let process = require('process'); // We can get process.argv

class Atom {
    constructor(id, type, mass, charge, seg, res, resName, name) {
	this.id = id;
	this.type = type;
	this.charge = charge;
	this.mass = mass;
	this.seg = seg;
	this.res = res;
	this.resName = resName;
	this.name = name;
	this.x = 0.0;
	this.y = 0.0;
	this.z = 0.0;
    }
}

const KIND = {
    lj: 0,
    bond: 1,
    angle: 2,
    dihedral: 3,
    improper: 4,
    crossterm: 5,
    bespoke: 6,
    exclusion: 7,
}

const KIND_NUM0 = KIND.bespoke;
const KIND_NUM = Object.keys(KIND).length;
const KIND_NAME = new Array(KIND_NUM);
Object.keys(KIND).forEach(key => KIND_NAME[KIND[key]] = key);


class MolecularStructure {
    constructor(psf, pdb, xsc, parameterFileList) {
	let psfData = fs.readFileSync(psf, 'utf8');
	let pdbData = fs.readFileSync(pdb, 'utf8');
	let xscData = fs.readFileSync(xsc, 'utf8');

	// Read the dimensions
	console.log('\nReading the dimensions from the XSC file:');
	this.boxDim = this.readXscOrtho(xscData);
	console.log('Dimensions:', this.boxDim.join(' '));

	// PSF data
	this.atoms = [];
	this.bonds = [];
	this.angles = [];
	this.dihedrals = [];
	this.impropers = [];
	this.crossterms = [];
	
	// Reading the PSF data
	console.log(`\nReading the PSF: ${psf}`);
	this.readPsfAtoms(psfData);
	console.log('atoms ' + this.atoms.length);
	const sectionList = ['bond', 'angle', 'dihedral', 'improper', 'crossterm'];
	for (const section of sectionList) {
	    this.readPsfSection(psfData, section);
	}
	console.log('bonds ' + this.bonds.length);
	console.log('angles ' + this.angles.length);
	console.log('dihedrals ' + this.dihedrals.length);
	console.log('impropers ' + this.impropers.length);
	console.log('crossterms ' + this.crossterms.length);

	// Read coordinates from the PDB
	console.log(`\nReading the PDB: ${pdb}`);
	const pdbNum = this.readPdbCoordinates(pdbData);
	console.log(`Read ${pdbNum} atom positions`);
	if (pdbNum != this.atoms.length) {
	    console.log(`ERROR! PSF ${psf} and PDB ${pdb} have a different number of atoms (${this.atoms.length} != ${pdbNum})`);
	    process.exit(1);
	}

	// Get the list of types
	const typeList = [];
	for (const atom of this.atoms) typeList.push(atom.type);
	
	// Get the unique types
	this.uniqueTypeList = [...new Set(typeList)].sort();
	const uniqueTypeList1 = [];
	// Put it in a list to be consistent with the other types of interaction
	for (const uniqueType of this.uniqueTypeList) uniqueTypeList1.push([uniqueType]);
	
	// Put these into indexable array
	this.terms = [];
	this.terms[KIND.lj] = [...this.atoms.keys()].map(i => [i]);
	this.terms[KIND.bond] = this.bonds;
	this.terms[KIND.angle] = this.angles;
	this.terms[KIND.dihedral] = this.dihedrals;
	this.terms[KIND.improper] = this.impropers;
	this.terms[KIND.crossterm] = this.crossterms;

	// This will give index of the unique parameter corresponding to each bonded term
	this.termParIndices = new Array(KIND_NUM0);
	for (let k = 0; k < KIND_NUM0; k++) {
	    this.termParIndices[k] = new Array(this.terms[k].length);
	}

	// What parameters will we need?
	console.log('\nOrganizing needed parameters:');
	// Unique interaction types
	this.unique = new Array(KIND_NUM);
	for (let k = 0; k < KIND_NUM0; k++) {
	    this.unique[k] = this.neededParameters(this.terms[k], this.termParIndices[k], k);
	    console.log(`${KIND_NAME[k]} types ${this.unique[k].length}`);
	}
	for (let k = KIND_NUM0; k < KIND_NUM; k++) this.unique[k] = [];

	// Initialize lists for available parameters from the parameter file
	this.availTypes = new Array(KIND_NUM);
	this.availParameters = new Array(KIND_NUM);
	for (let k = 0; k < KIND_NUM; k++) {
	    this.availTypes[k] = [];
	    this.availParameters[k] = [];
	}
	
	// Read the parameters from the parameter files
	// The results are added to this.availTypes and this.availParameters
	console.log('\nReading the CHARMM-style parameter files');
	for (const parameterFile of parameterFileList) {
	    let prmData = fs.readFileSync(parameterFile, 'utf8');
	    // The parameter file name is just given for diagnostics
	    this.readParameters(prmData, this.uniqueTypeList, parameterFile);
	}

	// How many parameters were read?
	for (let k = 0; k < KIND_NUM; k++) {
	    //console.log(`Found ${this.availTypes[k].length} potentially useful ${KIND_NAME[k]} parameters (first length: ${this.availParameters[k][0].length}) in the parameter files`);
	    if (this.availTypes[k].length > 0) 
		console.log(`Found ${this.availTypes[k].length} potentially useful ${KIND_NAME[k]} parameters (first length: ${this.availParameters[k][0].length})`);
	}

	// const k = KIND.dihedral;
	// for (let i = 0; i < this.availParameters[k].length; i++) {
	//     console.log(this.availParameters[k][i], this.availTypes[k][i]);
	// }

	// Match parameters to terms
	// The result is that this.parameters should be fully populated
	// this.parameters are lists of parameter values that are linked to by this.termParIndices
	this.linkParameters();

	// Dihedrals can have multiple terms with different multiplicities
	// Expand these into multiple terms in the terms array
	this.expandDihedrals();
	
	// Organize bespoke Lennard-Jones parameters
	this.parameters[KIND.bespoke] = [];
	this.terms[KIND.bespoke] = [];
	this.termParIndices[KIND.bespoke] = [];
	this.findBespokePairs(); // This is a slow step, O(N^2)

	// List of 1-2 and 1-3 bonded atoms
	console.log('\nGenerating exclusions');
	this.terms[KIND.exclusion] = this.generateExclusions(this.bonds);
	console.log(`Generated ${this.terms[KIND.exclusion].length} 1-3 exclusions`);
	// Exclusions allow for special nonbonded interactions (dummy, epsilon, Rmin/2)
	// No non-bonded interactions for 1-2 and 1-3
	// The signature for an exclusion is eps=0.0, RMin=-1.0
	this.parameters[KIND.exclusion] = [[0.0, -1.0]];
	this.unique[KIND.exclusion] = [['X']];
	this.termParIndices[KIND.exclusion] = new Array(this.terms[KIND.exclusion].length);
	this.termParIndices[KIND.exclusion].fill(0);

	// Special 1-4 LJ terms
	this.generateSpecialExclusions14(this.bonds);
	    
	// Organize the masses
	console.log('\nAtom properties');
	const massList = this.atoms.map(item => parseFloat(item.mass));
	[this.uniqueMassList, this.atomMassIndices] = this.compressProperty(massList);
	console.log('Unique masses: ' + this.uniqueMassList.length);
	console.log('Most common mass: ' + this.uniqueMassList[0])

	// Organize the charges
	const chargeList = this.atoms.map(item => parseFloat(item.charge));
	[this.uniqueChargeList, this.atomChargeIndices] = this.compressProperty(chargeList);
	console.log('Unique charges: ' + this.uniqueChargeList.length);
	console.log('Most common charge: ' + this.uniqueChargeList[0])	

	// Organize the types
	const typeList2 = this.atoms.map(item => item.type);
	[this.uniqueTypeList2, this.atomTypeIndices] = this.compressProperty(typeList);
	console.log('Unique types: ' + this.uniqueTypeList2.length);
	console.log('Most common type: ' + this.uniqueTypeList2[0])

	// Organize the resNames
	const resNameList = this.atoms.map(item => item.resName);
	[this.uniqueResNameList, this.atomResNameIndices] = this.compressProperty(resNameList);
	console.log('Unique residue names: ' + this.uniqueResNameList.length);
	console.log('Most common residue name: ' + this.uniqueResNameList[0])

	// Organize the atom names
	const nameList = this.atoms.map(item => item.name);
	[this.uniqueNameList, this.atomNameIndices] = this.compressProperty(nameList);
	console.log('Unique atom names: ' + this.uniqueNameList.length);
	console.log('Most common atom name: ' + this.uniqueNameList[0])
	
	// Convert Urey-Bradley angles to bonds
	// We should be careful that these UB bonds are not considered as bonds for exclusions
	// Doing this after generating the exclusions avoids this
	this.convertUreyBradleyAnglesToBonds(this.bonds);

	// It has not yet been compressed, so it's still possible to map the types in this.unique to this.parameters
	this.compressed = false;
    }

    // Given that the we need to subtract the standard Lennard-Jones for exclusions,
    // it is convenient to merge them with special 1-4 and bespoke Lennard-Jones
    mergeBespokeToExclusion() {
	console.log('\nMerging bespoke terms with exclusion terms');
	console.log(`Initially ${this.parameters[KIND.exclusion].length} exclusion parameters`);
	const offset = this.parameters[KIND.exclusion].length;
	for (let p = 0; p < this.parameters[KIND.bespoke].length; p++) {
	    let par = this.parameters[KIND.bespoke][p];
	    let types = this.unique[KIND.bespoke][p];
	    this.parameters[KIND.exclusion].push(par);
	    this.unique[KIND.exclusion].push(types);
	}

	// Put the bespoke parameters in the exclusion list
	let bespokeCount = 0;
	console.log(`Initially ${this.terms[KIND.exclusion].length} exclusion terms`);
	for (let t = 0; t < this.terms[KIND.bespoke].length; t++) {
	    // Check the corner case of both bespoke and exclusion or special 1-4
	    // The exclusion or special 1-4 should override the bespoke interaction
	    let term = this.terms[KIND.bespoke][t];
	    let j0 = this.search(this.parameters[KIND.exclusion], term.reverse());
	    // reverse mutates
	    let j1 = this.search(this.parameters[KIND.exclusion], term.reverse());

	    if (j0 < 0 && j1 < 0) {
		this.terms[KIND.exclusion].push(term);
		this.termParIndices[KIND.exclusion].push(this.termParIndices[KIND.bespoke][t] + offset);
		bespokeCount++;
	    }
	}
	console.log(`Merged ${bespokeCount} of ${this.terms[KIND.bespoke].length} bespoke Lennard-Jones terms into the exclusion term array`);
	console.log(`Any remaining bespoke terms were overridden by exclusions or special 1-4 Lennard-Jones terms`);
	console.log(`Finally ${this.parameters[KIND.exclusion].length} exclusion parameters`);
	console.log(`Finally ${this.terms[KIND.exclusion].length} exclusion terms`);

	// Clear the bespoke lists
	this.unique[KIND.bespoke] = [];
	this.parameters[KIND.bespoke] = [];
	this.terms[KIND.bespoke] = [];
	this.termParIndices[KIND.bespoke] = [];
	console.log('Cleared bespoke arrays');
    }

    // Merge impropers into the dihedrals
    mergeImproperToDihedral() {
	console.log('\nMerging impropers with regular dihedrals');
	console.log(`Initially, there are ${this.parameters[KIND.improper].length} impropers`);
	console.log(`Initially, there are ${this.parameters[KIND.dihedral].length} dihedrals`);

	// Create the new dihedral types
	const oldToNew = [];
	const startIndex = this.parameters[KIND.dihedral].length;
	for (let i = 0; i < this.parameters[KIND.improper].length; i++) {
	    let par = this.parameters[KIND.improper][i];
	    // Give impropers an multiplicity parameter of 0
	    this.parameters[KIND.dihedral].push([par[0], 0, par[1]]);
	    this.unique[KIND.dihedral].push(this.unique[KIND.improper][i]);
	    oldToNew[i] = startIndex + i;
	}

	// Add the improper terms to the dihedral term list
	for (let ti = 0; ti < this.terms[KIND.improper].length; ti++) {
	    const term = this.terms[KIND.improper][ti];
	    const parIndex = this.termParIndices[KIND.improper][ti];
	    this.terms[KIND.dihedral].push(term);
	    this.termParIndices[KIND.dihedral].push(oldToNew[parIndex]);
	}
	console.log(`After merging, there are ${this.parameters[KIND.dihedral].length} dihedrals`);
	
	// Delete the impropers
	this.unique[KIND.improper] = [];
	this.parameters[KIND.improper] = [];
	this.terms[KIND.improper] = [];
	this.termParIndices[KIND.improper] = [];
	console.log('Cleared improper arrays');
    }

    deleteZeroParameters() {
	console.log('\nDeleting parameters with zero force constants');
	let deletedSomething = false;

	for (let k = 0; k < KIND_NUM; k++) {
	    // Don't delete exclusions or Lennard-Jones interactions
	    if (k == KIND.exclusion || k == KIND.lj) continue;

	    // Find any zero force constants
	    const newParameters = [];
	    const deleteParIndices = [];
	    const oldToNew = new Array(this.parameters[k].length);
	    for (let i = 0; i < this.parameters[k].length; i++) {
		if (this.parameters[k][i][0] == 0.0) {
		    // Zero parameter
		    //console.log(`Zero ${KIND_NAME[k]} parameter ${i}: ${this.parameters[k][i].join(' ')}`);
		    oldToNew[i] = -1; // These are deleted
		    deleteParIndices.push(i);
		} else {
		    oldToNew[i] = newParameters.length;
		    newParameters.push(this.parameters[k][i]);
		}
	    }
	    //console.log(`Found ${deleteParIndices.length} ${KIND_NAME[k]} parameters with zero force constants`);

	    // Delete terms with zero force constants
	    const newTerms = [];
	    const newTermParIndices = [];
	    let deleteTermNum = 0;
	    for (let ti = 0; ti < this.terms[k].length; ti++) {
		let parIndex = this.termParIndices[k][ti];
		if (oldToNew[parIndex] >= 0) {
		    newTerms.push(this.terms[k][ti]);
		    newTermParIndices.push(oldToNew[parIndex]);
		} else {
		    deleteTermNum++;
		}
	    }

	    // Set the new terms and parameters
	    this.parameters[k] = newParameters;
	    this.terms[k] = newTerms;
	    this.termParIndices[k] = newTermParIndices;
	    if (deleteTermNum > 0) {
		console.log(`Deleted ${deleteParIndices.length} ${KIND_NAME[k]} parameters and ${deleteTermNum} ${KIND_NAME[k]} terms with zero force constants`);
		deletedSomething = true;
	    }
	    //console.log(`Deleted ${deleteTermNum.length} ${KIND_NAME[k]} terms linked to these parameters`);
	}

	if (!deletedSomething) {
	    console.log(`No terms with zero force constants were found. No terms were deleted.`);
	}
    }
    
    compressProperty(propertyList) {
	let uniqueList = [...new Set(propertyList)];

	// Link indices to values
	let indexList = new Array(propertyList.length);
	for (let i = 0; i < indexList.length; i++) {
	    let j = uniqueList.indexOf(propertyList[i]);
	    indexList[i] = j;
	}

	// Sort by number of usages
	let counts = new Array(uniqueList.length);
	counts.fill(0);
	for (const index of indexList) counts[index]++;

	// Make an array with count and index
	let indicesCounts = new Array(uniqueList.length); 
	for (let i = 0; i < indicesCounts.length; i++) indicesCounts[i] = [i, counts[i]];
	// Sort it by the count in descending order
	indicesCounts.sort( function(a,b) {return b[1] - a[1];} );

	// Reorder the parameter lists
	const sortUniqueList = new Array(uniqueList.length);
	const oldToNew = new Array(uniqueList.length);
	for (let i = 0; i < indicesCounts.length; i++) {
	    sortUniqueList[i] = uniqueList[indicesCounts[i][0]];
	    oldToNew[indicesCounts[i][0]] = i;
	}

	// Reassign the indices
	const sortIndexList = indexList.map( index => oldToNew[index] );
	return [sortUniqueList, sortIndexList];
    }
    
    readPsfAtoms(psfData) {
	const psf_lines = psfData.split(/\r?\n/);
	if (psf_lines[0].slice(0,3) != 'PSF') {
	    console.log('ERROR! Tried to read CHARMM PSF data with missing header (may not be a PSF)');
	    process.exit(1);
	}

	let reading = false;
	for (const line of psf_lines) {
	    if (reading) {		    
		let tokens = line.trim().split(/\s+/);
		// Quit if we left the atoms section.
		if (tokens.length >= 2 && tokens[1].match(/!N/)) { break }
		// Check that this is a valid line.
		if (tokens.length < 8) { continue }

		// Get the atom data.
		// Our ids are zero-based, rather than 1 based as in the psf.
		let id = parseInt(tokens[0], 10) - 1; 
		let seg = tokens[1];
		let res = tokens[2];
		let resName = tokens[3];
		let name = tokens[4];
		let type = tokens[5];
		let charge = tokens[6];
		let mass = tokens[7];

		// Add the new atom to the list.
		let atom = new Atom(id, type, mass, charge, seg, res, resName, name);
		this.atoms.push(atom);
	    }

	    // Find the atom records.
	    if (line.match(/!NATOM/)) {
		let tokens = line.trim().split(/\s+/);
		if (tokens[1].match(/^!NATOM/))  reading = true;
	    }
	}
    }

    readPsfSection(psfData, section) {
	// What kind of bonded terms are we looking for?
	let sectionName = 'NONE';
	let bodies = 1;
	let itemsPerLine = 1;
	let array = [];
	if (section.match(/^bond/)) {
	    sectionName = '!NBOND';
	    bodies = 2;
	    itemsPerLine = 4;
	    array = this.bonds;
	} else if (section.match(/^angle/)) {
	    sectionName = '!NTHETA';
	    bodies = 3;
	    itemsPerLine = 3;
	    array = this.angles;
	} else if (section.match(/^dihedral/)) {
	    sectionName = '!NPHI';
	    bodies = 4;
	    itemsPerLine = 2;
	    array = this.dihedrals;
	} else if (section.match(/^improper/)) {
	    sectionName = '!NIMPHI';
	    bodies = 4;
	    itemsPerLine = 2;
	    array = this.impropers;
	} else if (section.match(/^crossterm/)) {
	    sectionName = '!NCRTERM';
	    bodies = 8;
	    itemsPerLine = 1;
	    array = this.crossterms;
	} else {
	    console.log('ERROR! Unknown psf section ' + section)
	    return;
	}

	let reading = false;
	for (const line of psfData.split(/\r?\n/)) {
	    if (reading) {
		let tokens = line.trim().split(/\s+/);
		// Quit if we left this section.
		if (tokens.length >= 2 && tokens[1].match(/!N/)) { break }
		// Ignore blank lines
		if (line.length == 0 || tokens.length == 0) { continue }
		// The lines need not be complete line (itemsPerLine items)
		// But it must have at a least one complete term
		// and a number of items that is a multiple of bodies
		if (tokens.length < bodies || tokens.length % bodies != 0) {
		    console.log(`Warning! Invalid line in ${sectionName} section of PSF file: '${line}'`);
		}

		// Read the items.
		let n = 0;
		let l = [];
		for (const tok of tokens) {
		    // Get the serial number.
		    // PSF files are one-based, but
		    // we will keep the arrays zero-based;
		    l.push(parseInt(tok)-1);

		    if (l.length == bodies) {
			// This is a complete item.
			// Add it to the array.
			array.push(l);
			l = [];
		    }
		}
	    }

	    // Find the atom records.
	    if (line.match(/!N/)) {
		let tokens = line.trim().split(/\s+/);
		if (tokens[1].match('^'+sectionName)) reading = true;
	    }
	}
    }
    
    readPdbCoordinates(pdbData) {
	let pdbLines = pdbData.split(/\r?\n/);
	let id = 0;
	for (const line of pdbLines) {
	    if ( line.slice(0,4) == 'ATOM' || line.slice(0,6) == 'HETATM') {
		const x = parseFloat(line.slice(30,38));
		const y = parseFloat(line.slice(38,46));
		const z = parseFloat(line.slice(46,54));

		if (id < this.atoms.length) {
		    this.atoms[id].x = x;
		    this.atoms[id].y = y;
		    this.atoms[id].z = z;
		}
		id++;
	    }
	}
	return id;
    }

    // Assumes orthogonal (a, 0, 0), (0, b, 0), (0, 0, c) structure
    readXscOrtho(xscData) {
	let xscLines = xscData.split(/\r?\n/);
	let lx = 0;
	let ly = 0;
	let lz = 0;
	for (const line of xscLines) {
	    let l = line.trim();
	    if (l.length == 0) continue;
	    if (l.match(/^#/)) continue;
	    let tokens = l.split(/\s+/);
	    lx = parseFloat(tokens[1]);
            ly = parseFloat(tokens[5]);
            lz = parseFloat(tokens[9]);
	    //console.log(l);
	}
	return [lx, ly, lz];
    }

    
    // Determine what parameters are needed from the unique sets 
    neededParameters(termArray, termParIndex, kind) {
	const uniqueArray = [];

	for (let i = 0; i < termArray.length; i++) { 
	    const types = [];
	    const term = termArray[i];
	    for (const atomIndex of term) {
		types.push(this.atoms[atomIndex].type);
	    }
	    
	    // Impose a canonical order.
	    const first = this.atoms[term[0]].type;
	    const last = this.atoms[term[term.length-1]].type;
	    // All of the terms, except possibly impropers, can be reversed without change
	    // We sort them alphabetically
	    if (kind != KIND.improper && last < first) types.reverse();

	    // Have we already found this bond?
	    const j = this.search(uniqueArray, types);
	    if (j<0) {
		// This is a new bond type.
		uniqueArray.push(types);
		termParIndex[i] = uniqueArray.length-1;
	    } else {
		termParIndex[i] = j;
	    }
	}
	return uniqueArray;
    }

    // Convert a number to fixed digit representation, while deleting trailing zeros
    fixedDigits(n, digits) {
	const str = n.toFixed(digits);
	const decimalPos = str.indexOf('.');
	let cut = str.length;
	for (let i = str.length-1; i >= decimalPos+2; i--) {
	    if (str[i] != '0') {
		// We stop at the first nonzero character
		break;
	    }
	    cut = i;
	}
    
	return str.slice(0,cut);
    }
    
    // Search for a matching list in a list of lists
    // An index of -1 indicates that a match was not found
    search(haystackList, needleList) {
	for (let j = 0; j < haystackList.length; j++) {
	    let haystackItem = haystackList[j];

	    // The two lists must have the same lengths to be equal.
	    if (haystackItem.length == needleList.length) {
		let equal = 1;
		// All items of the two lists must be equal.
		for (let i = 0; i < haystackItem.length; i++) {
		    if (haystackItem[i] != needleList[i]) {
			equal = 0;
			break;
		    }
		}
		if ( equal ) return j;
	    }
	}
	return -1;
    }

    // Search for a matching list in a list of lists
    // An index of -1 indicates that a match was not found
    searchAll(haystackList, needleList) {
	const ret = [];
	for (let j = 0; j < haystackList.length; j++) {
	    let haystackItem = haystackList[j];

	    // The two lists must have the same lengths to be equal.
	    if (haystackItem.length == needleList.length) {
		let equal = 1;
		// All items of the two lists must be equal.
		for (let i = 0; i < haystackItem.length; i++) {
		    if (haystackItem[i] != needleList[i]) {
			equal = 0;
		    }
		}
		if ( equal ) ret.push(j);
	    }
	}
	return ret;
    }

    
    // Search for a matching list in a list of lists
    // An index of -1 indicates that a match was not found
    searchWild(haystackList, needleList, wildValue) {
	for (let j = 0; j < haystackList.length; j++) {
	    let haystackItem = haystackList[j];

	    // The two lists must have the same lengths to be equal.
	    if (haystackItem.length == needleList.length) {
		let equal = 1;
		// All items of the two lists must be equal.
		for (let i = 0; i < haystackItem.length; i++) {
		    if (haystackItem[i] != needleList[i] && haystackItem[i] != wildValue) {
			equal = 0;
			break;
		    }
		}
		if ( equal ) return j;
	    }
	}
	return -1;
    }

    // Read parameters from CHARMM-format prm or str data
    // fileName is included for diagnostics only
    // results are added to this.availTypes and this.availParameters
    // uniqueTypeList allows us to exclude definitely useless parameters
    readParameters(prmData, uniqueTypeList, fileName) {
	let section = ''; // What section are we in?
	let bodies = 0; // How many bodies is the term (number of types on the line)
	let kind = -1; // What "KIND" enumeration?
	let pars = []; // How many parameters (it's a list since sometimes there are variations)
	let reading = false;
	const prmLines = prmData.split(/\r?\n/);
	for (let line of prmLines) {
	    // Ignore comments.
	    line = line.trim();
	    if (line.length == 0) continue;
	    if (line.match(/^!/) || line.match(/^\*/) || line.match(/^HBOND/)) continue;

	    // Throw out any part of the line after the comment character
	    const cut = line.indexOf('!');
	    if (cut >= 0) line = line.slice(0, cut).trim();

	    // Split into tokens
	    let tokens = line.split(/\s+/);
	    
	    // First we check if we've hit a section header
	    // If not, we check if we need this set of parameter values
	    // If so, we add it to the appropriate list
	    if (line.match(/^ATOMS/) || line.match(/^MASS/) || line.match(/^RESI/) || line.match(/^read/) || line.match(/^END/)) {
		section = '';
	    } else if (line.match('^BONDS')) {
		section = 'BONDS';
		kind = KIND.bond;
		bodies = 2;
		pars = [2];
	    } else if (line.match('^ANGLES')) {
		section = 'ANGLES';
		kind = KIND.angle;
		bodies = 3;
		pars = [2,4]; // can include Urey-Bradley
	    } else if (line.match('^DIHEDRALS')) {
		section = 'DIHEDRALS';
		kind = KIND.dihedral;
		bodies = 4;
		pars = [3];
	    } else if (line.match('^IMPROPER')) {
		section = 'IMPROPERS';
		kind = KIND.improper;
		bodies = 4;
		pars = [3];
	    } else if (line.match('^NONBONDED')) {
		section = 'NONBONDED';
		kind = KIND.lj;
		bodies = 1;
		pars = [3,6]; // Can include special 1-4 parameters
	    } else if (line.match('^CMAP')) {
		section = 'CMAP';
		kind = KIND.crossterm;
		bodies = 8;
		pars = []; // The number varies based on the n
	    } else if (line.match('^NBFIX')) {
		section = 'NBFIX';
		kind = KIND.bespoke;
		bodies = 2;
		pars = [2];
	    } else if (section.length > 0 && kind == KIND.crossterm && !isNaN(tokens[0]) && reading) {
		// Keep collecting the CMAP map values associated with the last parameter
		let last = this.availParameters[kind].length - 1;
		let par = tokens.map( item => parseFloat(item) );
		this.availParameters[kind][last].push(...par); }
	    else if (section.length > 0 && kind == KIND.crossterm && !isNaN(tokens[0])) {
		// This is a crossterm that we definitely aren't using
		// Do nothing
	    } else if (section.length > 0) {
		// This is a regular line of a section we are reading.
		// Find if the types in this line match the needed parameters.
		if (tokens.length == 1) continue;
		// Make sure this line has enough tokens.
		if (tokens.length < bodies + 1) {
		    console.log(`Warning! Invalid CHARMM ${section} parameter line in file ${fileName}: ${tokens}`);
		    continue;
		}

		// Does this line have any matching types?
		let types = tokens.slice(0,bodies)
		let match = false;
		for (const type of types) {
		    if (uniqueTypeList.includes(type)) {
			match = true;
			break;
		    }
		}

		// Push the types and the parameters
		// Ignore lines that don't include any of the types we are interested in
		if (match) {
		    this.availTypes[kind].push(types);
		    let par = tokens.slice(bodies).map( item => parseFloat(item) );

		    // Check that we have the right number of parameters
		    if (kind != KIND.crossterm && !pars.includes(par.length)) {
			console.log(`ERROR! Invalid CHARMM ${section} parameter line in file ${fileName}: ${tokens}`);
			process.exit(1);
		    }
		    
		    // Modify the presentation of the parameters
		    if (kind == KIND.lj) {
			// CHARMM has an unused first field and Rmin/2 instead of Rmin
			// Also epsilon is negative for some reason
			// Remove the unused parameter 
			if (par.length == 6) {
			    // special 1-4 LJ terms
			    par = [Math.abs(par[1]), 2.0*par[2], Math.abs(par[4]), 2.0*par[5]];
			} else {
			    par = [Math.abs(par[1]), 2.0*par[2]];
			}			    
		    } else if (kind == KIND.bespoke) {
			// CHARMM epsilon is recorded as a negative value
			par = [Math.abs(par[0]), par[1]];
		    } else if (kind == KIND.improper) {
			// There's an unused field
			par = [par[0], par[2]];
		    }
		    
		    this.availParameters[kind].push(par);
		    reading = true;
		} else {
		    reading = false;
		}
	    } // End conditional for line kind
	} // End line loop
    } // End readParameters function
    
    linkParameters() {
	// Try to find the parameters corresponding to all terms
	console.log('\nSearching for parameters');
	this.parameters = [];
	for (let k = 0; k < KIND_NUM0; k++) {
	    this.parameters[k] = [];
	    for (const uniqueTypes of this.unique[k]) {
		const types = uniqueTypes.slice();
		let index_list = this.searchAll(this.availTypes[k], types);
		// It wasn't found. Try to reverse the types
		if (index_list.length == 0)
		    index_list = this.searchAll(this.availTypes[k], types.reverse());
		// It wasn't found. Try with wildcards
		if (index_list.length == 0) {
		    let index = this.searchWild(this.availTypes[k], types, 'X');
		    if (index >= 0) {
			console.log(`Matched ${KIND_NAME[k]} ${types} to wildcard ${this.availTypes[k][index]}`);
			index_list = [index];
		    }
		}
		// It wasn't found. Try with wildcards reversed.
		if (index_list.length == 0) {
		    let index = this.searchWild(this.availTypes[k], types.reverse(), 'X');
		    if (index >= 0) {
			console.log(`Matched ${KIND_NAME[k]} ${types} to wildcard ${this.availTypes[k][index]}`);
			index_list = [index];
		    }
		}

		// We can't find this parameter
		if (index_list.length == 0) {
		    console.log(`ERROR! No parameters for ${KIND_NAME[k]} ${types}`);
		    process.exit(1);
		} else {
		    // Dihedrals can have multiple sets of parameters
		    if (k == KIND.dihedral) {
			let par_set = [];
			for (const index of index_list)
			    par_set = par_set.concat(this.availParameters[k][index]);

			if (index_list.length > 1) 
			    console.log(`Found multiple multiplicities for dihedral ${types}: ${par_set}`);
			this.parameters[k].push(par_set);
		    } else {
			const last = index_list[index_list.length-1];
			const par = this.availParameters[k][last];
			this.parameters[k].push(par);
			if (index_list.length > 1) 
			    console.log(`Warning! Multiple sets of parameters for ${KIND_NAME[k]} ${types}. Using the newest set: ${par}`);
		    }
		}
	    }
	    console.log(`Found parameters for all ${this.unique[k].length} unique ${KIND_NAME[k]} terms`);
	}
    }

    expandDihedrals() {
	console.log('\nExpanding dihedrals with multiple terms');
	const k = KIND.dihedral;

	// Clone the parameters and unique types
	const newParameters = [];
	const newUnique = [];
	for (let p = 0; p < this.parameters[k].length; p++) {
	    newParameters[p] = this.parameters[k][p];
	    newUnique[p] = this.unique[k][p];
	}
	
	// Expand the parameter and unique arrays
	const expandIndexList = []; // Original parameter index 
	const expandNewListList = []; // List of lists of new indices
	for (let p = 0; p < this.parameters[k].length; p++) {
	    const par = this.parameters[k][p];
	    const unique = this.unique[k][p];
	    const parCount = par.length/3;

	    const newIndices = [];
	    for (let i = 1; i < parCount; i++) {
		const parK = par[3*i];
		const parMulti = par[3*i+1];
		const parDelta = par[3*i+2];
		
		newIndices.push(newParameters.length);
		newParameters.push([parK, parMulti, parDelta]);
		newUnique.push(unique);
	    }
	    if (parCount > 1) {
		// Cut the original parameter list down to the first three
		newParameters[p] = par.slice(0,3);

		// Keep track of which parIndices need to be expanded
		expandIndexList.push(p);
		// The new parameters are added at the end of the arrays
		expandNewListList.push(newIndices);
		console.log(`Expanding dihedral ${unique} to ${parCount} terms`);
	    }
	}

	// Set the new parameter and unique lists
	console.log(`Initially ${this.parameters[k].length} unique dihedral types`);
	this.parameters[k] = newParameters;
	this.unique[k] = newUnique;
	console.log(`Now ${this.parameters[k].length} unique dihedral types`);

	// Now we need to add new terms for the expanded dihedrals
	console.log(`Initially ${this.terms[k].length} dihedral terms`);

	// Clone the terms and termParIndices
	const newTerms = [];
	const newTermParIndices = [];
	for (let ti = 0; ti < this.terms[k].length; ti++) {
	    newTerms[ti] = this.terms[k][ti];
	    newTermParIndices[ti] = this.termParIndices[k][ti];
	}
	
	for (let ti = 0; ti < this.terms[k].length; ti++) {
	    const parIndex = this.termParIndices[k][ti];
	    const j = expandIndexList.indexOf(parIndex);

	    if (j >= 0) {
		// We need to expand this one
		for (const expandIndex of expandNewListList[j]) {
		    // Just copy this term
		    newTerms.push(this.terms[k][ti]);
		    // Add the new parameter indices
		    newTermParIndices.push(expandIndex);
		}
	    }
	}

	// Set the new terms
	this.terms[k] = newTerms;
	this.termParIndices[k] = newTermParIndices;
	console.log(`Now ${this.terms[k].length} dihedral terms`);
    }
    
    generateExclusions(bonds) {
	// Generate a bond list for each atom
	const bondList = new Array(this.atoms.length);
	for (let a = 0; a < this.atoms.length; a++) bondList[a] = [];
	for (let b = 0; b < bonds.length; b++) {
	    const a0 = bonds[b][0];
	    const a1 = bonds[b][1];

	    bondList[a0].push(a1);
	    bondList[a1].push(a0);
	}

	// Generate lists containing nearest and next nearest neighbors
	const neighList = new Array(this.atoms.length);
	const excludeList = new Array(this.atoms.length);
	for (let a = 0; a < this.atoms.length; a++) {
	    // Concatenate the bondList of the neighbors
	    // Notice that Array() would be bad if there is only one bond
	    neighList[a] = [...bondList[a]];
	    for (const neighIndex of bondList[a]) {
		neighList[a].push(...bondList[neighIndex]);
	    }
	    excludeList[a] = [...new Set(neighList[a])].sort((a, b) => (a - b));
	}

	// Generate pairs
	const exclusions13 = [];
	for (let a = 0; a < this.atoms.length; a++) {
	    for (let n = 0; n < excludeList[a].length; n++) {
		let a1 = excludeList[a][n];
		// We don't want to double count or self-interaction
		if (a < a1) {
		    exclusions13.push([a, a1]);
		}
	    }
	}
	//console.log(exclusions13);
	return exclusions13;
    }

    generateSpecialExclusions14(bonds) {
	console.log('\nGenerating special 1-4 LJ interactions');
	console.log(`Initially ${this.terms[KIND.exclusion].length} exclusion terms with ${this.parameters[KIND.exclusion].length} unique parameters`);

	// Find which atoms have extra Lennard-Jones parameters
	const typeList14 = [];
	const parList14 = [];
	for (let i = 0; i < this.parameters[KIND.lj].length; i++) {
	    if (this.parameters[KIND.lj][i].length == 4) {
		// Splice removes these parameters from the lj array 
		const special = this.parameters[KIND.lj][i].splice(2);
		const type = this.unique[KIND.lj][i];
		typeList14.push(type[0]);
		parList14.push(special);
	    } else if(this.parameters[KIND.lj][i].length != 2) {
		console.log(`Warning! The length of the LJ parameter set for ${this.unique[KIND.lj][i]} is not 2 or 4 (${this.parameters[KIND.lj][i].length})`);
	    }
	}
	console.log(`Found ${typeList14.length} special 1-4 LJ interaction types`);

	// Find atoms with this type
	let specialCount = 0;
	const unique14 = [];
	const parameters14 = [];
	const terms14 = [];
	const termParIndices14 = [];
	console.log('Special 1-4 types: ' + typeList14.join(' '));
	for (let a = 0; a < this.atoms.length; a++) {
	    const type = this.atoms[a].type;
	    const index14 = typeList14.indexOf(type);

	    // This atom doesn't have 1-4 special LJ terms
	    if (index14 < 0) continue;
	    
	    // Get the special 1-4 parameters for this atom
	    const ljPar = parList14[index14];
	    
	    // Find their 1-4 interactions
	    let neighList14 = this.get14(a, bonds);
	    for (const neigh of neighList14) {
		// Check if we've already seen these two atoms
		const term = [a,neigh].sort();
		const termIndex = this.search(terms14, term);
		// Make sure we don't add the term more than once
		if (termIndex >= 0) continue;
		
		// Add the exclusion term
		terms14.push(term);
		
		// Get the parameters for the neighbor
		const neighType = this.atoms[neigh].type;

		// Have we already seen this pair of types?
		const types = [type, neighType].sort();
		let uniqueIndex = this.search(unique14, types);		    
		if (uniqueIndex >= 0) {
		    // We've already seen this pair of types
		    termParIndices14.push(uniqueIndex);
		} else {
		    // We haven't seen this pair of types yet.
		    // EDGE CASE: First we look to see if the neighbor type is also in the 1-4 list
		    // Actually this edge case is common in CHARMM proteins
		    const neighIndex14 = typeList14.indexOf(neighType);
		    let neighPar = null;
		    if (neighIndex14 >= 0) {
			// This is an edge case (neighbor also has special 1-4)
			neighPar = parList14[neighIndex14];
			console.log('Combining two special 1-4 LJ parameters using Lorentz-Bethelot rules:', type, neighType);
		    } else {
			// This is the usual case (neighbor also has standard LJ)
			// Find the neighbor type in the LJ parameters
			const neighParIndex = this.search(this.unique[KIND.lj], [neighType]);
			neighPar = this.parameters[KIND.lj][neighParIndex];
			console.log(`Combining special 1-4 LJ parameters for ${type} with using Lorentz-Bethelot rules:`, type, neighType);
		    }

		    // Calculate the mixed parameters using Lorentz-Berthelot rules
		    const mixEps = this.fixedDigits(Math.sqrt(neighPar[0] * ljPar[0]), 8);
		    const mixRMin = this.fixedDigits(0.5*(neighPar[1] + ljPar[1]), 8);
		    const mixPar = [mixEps, mixRMin];
		    // Add these parameters to the 1-4 unique types list
		    termParIndices14.push(parameters14.length);
		    unique14.push(types);
		    parameters14.push(mixPar);
		}
		specialCount++;
	    } // End loop over 1-4 neighbors
	} // End loop over atoms


	// We need to merge our 1-4 data into unique and parameters arrays for exclusions
	const oldToNew = new Array(parameters14.length);
	for (let i = 0; i < parameters14.length; i++) {
	    oldToNew[i] = this.parameters[KIND.exclusion].length;
	    this.unique[KIND.exclusion].push(unique14[i]);
	    this.parameters[KIND.exclusion].push(parameters14[i]);
	}

	// We need to merge our 1-4 data into term and termParIndices arrays for exclusions
	for (let i = 0; i < terms14.length; i++) {
	    this.terms[KIND.exclusion].push(terms14[i]);
	    const parIndex = oldToNew[termParIndices14[i]];
	    this.termParIndices[KIND.exclusion].push(parIndex);
	}

	console.log(`Generated ${specialCount} special exclusions for ${unique14.length} pairs of types`);
	console.log(`Now ${this.terms[KIND.exclusion].length} exclusion terms with ${this.parameters[KIND.exclusion].length} unique parameters`);
	
	return specialCount;
    }

    getNeighList(atomIndex, bonds) {
	// Get the 1-2 atom indices
	let neighList = [];
	for (let b = 0; b < bonds.length; b++) {
	    if (bonds[b][0] == atomIndex) {
		neighList.push(bonds[b][1]); 
	    } else if (bonds[b][1] == atomIndex) {
		neighList.push(bonds[b][0]); 
	    }
	}
	return neighList;
    }
    
    get14(atomIndex, bonds) {
	let neighList12 = this.getNeighList(atomIndex, bonds);
	let seenList = neighList12.concat(atomIndex);

	// Generate the raw 1-3 neighbor list
	let rawList13 = [];
	for (const neigh of neighList12)
	    rawList13.push(...this.getNeighList(neigh, bonds));

	// Remove repeated neighbors
	let neighList13 = [];	    
	for (const neigh of rawList13) {
	    if (!seenList.includes(neigh)) neighList13.push(neigh);
	}
	// Add the 1-3 neighbors to the seen list
	seenList.push(...neighList13);
	
	// Generate the raw 1-4 neighbor list
	let rawList14 = [];
	for (const neigh of neighList13)
	    rawList14.push(...this.getNeighList(neigh, bonds));

	// Remove repeated neighbors
	let neighList14 = [];	    
	for (const neigh of rawList14) {
	    if (!seenList.includes(neigh)) neighList14.push(neigh);
	}

	return neighList14;
    }

    convertUreyBradleyAnglesToBonds() {
	console.log('\nConverting Urey-Bradley angles to bonds');
	console.log('Initial unique bond types: ' + this.parameters[KIND.bond].length);
	console.log('Initial bond terms: ' + this.terms[KIND.bond].length);
	
	let ureyBradleyCount = 0;
	let newBondTerms = 0;
	// Loop over the unique parameter index
	for (let ai = 0; ai < this.parameters[KIND.angle].length; ai++) {
	    if (this.parameters[KIND.angle][ai].length == 4) {
		// This is a Urey-Bradley angle
		// Get the UB parameters and remove them from the angle
		ureyBradleyCount++;
		const ub = this.parameters[KIND.angle][ai].splice(2);
		// Create a new set of bond parameters for this
		this.parameters[KIND.bond].push(ub);
		let typeList = this.unique[KIND.angle][ai];
		this.unique[KIND.bond].push([typeList[0], typeList[2]]);
		
		// Find all angle terms with this angle type
		const angleTerms = [];
		for (let ti = 0; ti < this.termParIndices[KIND.angle].length; ti++) {
		    if (this.termParIndices[KIND.angle][ti] == ai)
			angleTerms.push(this.terms[KIND.angle][ti]);
		}
		for (const angle of angleTerms) {
		    // Add the corresponding 1-3 bonded terms between atom 0 and 2
		    const bond = [angle[0], angle[2]];
		    this.terms[KIND.bond].push(bond);
		    // We just pushed this parameter, so it is the last in the unique parameter list
		    this.termParIndices[KIND.bond].push(this.parameters[KIND.bond].length - 1);
		    newBondTerms++;
		}
	    }
	}
	console.log(`Found ${ureyBradleyCount} Urey-Bradley unique angle types`);
	console.log(`Generated ${newBondTerms} new bond terms to represent them`);
	console.log('Final unique bond types: ' + this.parameters[KIND.bond].length);
	console.log('Final bond terms: ' + this.terms[KIND.bond].length);
    }


    findBespokePairs() {
	// Add the bespoke parameters
	console.log('\nAssigning bespoke (NBFIX) parameters');
	const bespokeTypeList = [];
	for (let i = 0; i < this.availTypes[KIND.bespoke].length; i++) {
	    const type0 = this.availTypes[KIND.bespoke][i][0];
	    const type1 = this.availTypes[KIND.bespoke][i][1];
	    if (this.uniqueTypeList.includes(type0) && this.uniqueTypeList.includes(type1)) {
		const typeList = [type0, type1].sort();
		const ti = this.search(bespokeTypeList, typeList)
		if (ti < 0) {
		    bespokeTypeList.push(typeList);
		    this.parameters[KIND.bespoke].push(this.availParameters[KIND.bespoke][i]);
		    this.unique[KIND.bespoke].push(typeList);
		}
	    }
	}
	console.log(`Found ${this.parameters[KIND.bespoke].length} active set of bespoke Lennard-Jones parameters`);

	// Find pairs that match the bespoke parameters
	for (let i = 0; i < this.atoms.length; i++) {
	    for (let j = i; j < this.atoms.length; j++) {
		let types = [];
		let atomIndices = [];
		// Sort the types
		if (this.atoms[i].type <= this.atoms[j].type) {
		    atomIndices = [i, j];
		    types = [this.atoms[i].type, this.atoms[j].type];
		} else {
		    atomIndices = [j, i];
		    types = [this.atoms[j].type, this.atoms[i].type];
		}
		
		const ti = this.search(bespokeTypeList, types);
		if (ti >= 0) {
		    this.terms[KIND.bespoke].push(atomIndices);
		    this.termParIndices[KIND.bespoke].push(ti);
		}
	    }
	}
	console.log(`Found ${this.terms[KIND.bespoke].length} bespoke Lennard-Jones pairs`);	
    }

    compress() {
    	this.compressed = true;

	console.log('\nCompressing identical sets of parameter values');
	for (let k = 0; k < KIND_NUM; k++) {
	    // Find identical parameters
	    let uniqueParameters = [];
	    let uniqueMap = [];
	    for (let pi = 0; pi < this.parameters[k].length; pi++) {
		let j = this.search(uniqueParameters, this.parameters[k][pi]);

		if (j < 0) {
		    uniqueParameters.push(this.parameters[k][pi]);
		    uniqueMap.push([pi]);
		} else {
		    uniqueMap[j].push(pi);
		}
	    }
	    
	    if (this.parameters[k].length != uniqueMap.length) {
		console.log(`Compressing ${this.parameters[k].length} sets of ${KIND_NAME[k]} parameter values into ${uniqueMap.length} unique sets`);
	    }

	    // Get identical types.
	    for (const indices of uniqueMap) {
		const typeIdentical = indices.map( item => this.unique[k][item] );
		//console.log(typeIdentical);
	    }

	    // Set the new parameters
	    this.parameters[k] = uniqueParameters;
	    
	    // Set the parameter indices to use the new unique values
	    for (let ti = 0; ti < this.termParIndices[k].length; ti++) {
		let oldParIndex = this.termParIndices[k][ti];
		
		// Find the unique index for this parameter index
		let newParIndex = -1;
		for (let ui = 0; ui < uniqueMap.length; ui++) {
		    if (uniqueMap[ui].includes(oldParIndex)) {
			newParIndex = ui;
			break;
		    }
		}
		
		if (newParIndex < 0) {
		    console.log(`ERROR! In compress(): parameter index ${oldParIndex} not found`);
		    process.exit(1);
		} else {
		    this.termParIndices[k][ti] = newParIndex;
		}
	    }
	} // End kind loop for compression

	// Sort parameters by usage so that more popular parameters use fewer characters
	this.sortParameters();

	// How many parameters
	console.log('\nNumber of parameters of each type after compression');
	for (let k = 0; k < KIND_NUM; k++) {
	    console.log(`${KIND_NAME[k]} parameters: ${this.parameters[k].length}`);
	}
    }

    // Sort parameters by usage so that more popular parameters use fewer characters
    sortParameters() {
	this.compressed = true;
	console.log('\nSorting parameters by usage');
	
	for (let k = 0; k < KIND_NUM; k++) {
	    if (this.parameters[k].length == 0) continue;
	    
	    let counts = new Array(this.parameters[k].length);
	    counts.fill(0);

	    // Count number of usages
	    for (const parIndex of this.termParIndices[k]) counts[parIndex]++;

	    // Make an array with count and index
	    let indicesCounts = new Array(this.parameters[k].length); 
	    for (let i = 0; i < indicesCounts.length; i++) indicesCounts[i] = [i, counts[i]];
	    // Sort it by the count in descending order
	    indicesCounts.sort( function(a,b) {return b[1] - a[1];} );
	    console.log('Most common ' + KIND_NAME[k] + ' parameter ' + indicesCounts[0][0] + ' (' +  indicesCounts[0][1] + ' times)');

	    // Reorder the parameter lists
	    const newParameters = new Array(this.parameters[k].length);
	    const oldToNew = new Array(this.parameters[k].length);
	    for (let i = 0; i < indicesCounts.length; i++) {
		newParameters[i] = this.parameters[k][indicesCounts[i][0]];
		oldToNew[indicesCounts[i][0]] = i;
	    }
	    this.parameters[k] = newParameters;

	    // Reassign the term indices
	    this.termParIndices[k] = this.termParIndices[k].map( parIndex => oldToNew[parIndex] );
	}
	console.log('Parameter indices have been renumbered to give lower numbers to those used the most');
    }

    // The most common fragments should have the smallest ids
    compressFragments() {
	// console.log('length monomer', this.atomMonomerList.length);
	// console.log('length fragment', this.atomFragmentList.length);
	if (this.atomMonomerList != null) {
	    let uniqueMonomers = null;	
	    [uniqueMonomers, this.atomMonomerList] = this.compressProperty(this.atomMonomerList);
	    console.log('Unique monomers: ' + uniqueMonomers.length);
	}

	if (this.atomFragmentList != null) {
	    let uniqueFragments = null;
	    [uniqueFragments, this.atomFragmentList] = this.compressProperty(this.atomFragmentList);
	    console.log('Unique fragments: ' + uniqueFragments.length);
	}
    }
    
    // Create fragment list by residue
    generateMonomersByResidue() {
	console.log('\nGenerating monomer list by residue');
	//Make monomer list
	let monomer = 0;
	let lastSeg = this.atoms[0].seg;
	let lastRes = this.atoms[0].res;
	const monomerList = [];
	for (const atom of this.atoms) {
	    if (atom.seg != lastSeg || atom.res != lastRes) {
		monomer++;
		lastRes = atom.res;
		lastSeg = atom.seg;
	    }
	    monomerList.push(monomer);
	}
	console.log(`Found ${monomer+1} monomers (residues)`);

	// Don't store a monomer list if there is only one monomer
	if (monomer != 0) {
	    this.atomMonomerList = monomerList;
	} else {
	    console.log(`The monomer list will not be written since there is only one monomer.`);
	}
    }

    // Create fragment list by making a bond graph
    // Fragments are zero-based like everything else
    generateFragmentsByBondGraph() {
	console.log('\nGenerating fragment list using the bond graph');
	
	// We need to clone the bond array since getGraphComponents mutates the list
	const bondList = [];
	for (const bond of this.bonds) {
	    const b = [...bond];
	    bondList.push(b);
	}
	const subgraphList = this.getGraphComponents(bondList);
	let nextFragment = subgraphList.length;
	const molecularFragments = subgraphList.length;
	console.log(`Found ${molecularFragments} molecular fragments (at least two atoms)`);

	// Find the largest cluster
	// Find the mean cluster size
	let sizeSum = 0.0;
	let maxComponentSize = 0;
	for (let gi = 0; gi < subgraphList.length; gi++) {
	    const size = subgraphList[gi].length;
	    sizeSum += size;
	    if (size > maxComponentSize) {
		maxComponentSize = size;
	    }
	}
	if (subgraphList.length > 0) {
	    console.log(`Mean multi-atom fragment size: ${sizeSum/subgraphList.length}`);
	    console.log(`Largest multi-atom fragment: ${maxComponentSize}`);
	}

	// Make a list of fragment indices for each atom
	// Solo atoms will not be found and need their own fragment
	const fragList = [];
	for (let i = 0; i < this.atoms.length; i++) {
	    // Which, subgraph contains this atom
	    let fragmentId = -1;
	    for (let gi = 0; gi < subgraphList.length; gi++) {
		const j = subgraphList[gi].indexOf(i);
		// Was this atom found?
		if (j >= 0) {
		    fragmentId = gi;
		    break;
		}
	    }

	    // If the atom was not found, it is a solo atom (not in the bond list)
	    // Give it a unique fragment ID
	    if (fragmentId < 0) {
		fragmentId = nextFragment;
		nextFragment++;
	    } 

	    // A fragment number is pushed for all atoms
	    fragList.push(fragmentId);
	}
	console.log(`Found ${nextFragment - molecularFragments} free atoms or ions`);
	console.log(`Found a total of ${nextFragment} fragments`);
	//console.log('Length of atomFragmentList', this.atomFragmentList.length);

	if (nextFragment > 1) {
	    this.atomFragmentList = fragList;
	} else {
	    console.log(`The fragment list will not be written since there is only one fragment.`);
	}
    }

    // Merge neighbor lists until only we are left with graph components
    // A component is a connected subgraph that is not part of any larger connected subgraph
    // Mutates the input list
    getGraphComponents(neighListList) {
	// Initial size
	let n = neighListList.length;
	
	// Continue to merge the lists until they stop changing
	let n0 = -1;
	while (n != n0) {
	    n0 = n;
	    neighListList = this.mergeFirstOverlap(neighListList);
	    n = neighListList.length;
	}
	return neighListList;
    }

    
    // Take a list of lists of atoms that neighbor each other
    // Merge the first list found that overlaps
    mergeFirstOverlap(neighListList) {
	// Compare all pairs of neighbor lists
	// We don't have to worry about the changing size
	// of the lists because we return after merging.
	for (let i = 0; i < neighListList.length-1; i++) {
	    for (let j = i + 1; j < neighListList.length; j++) {
		if (this.checkOverlap(neighListList[i], neighListList[j])) {
		    // They overlap. Merge the lists.
		    const lij = [...new Set([...neighListList[i] ,...neighListList[j]])];
		    // Remove the list higher in the array (j)
		    neighListList.splice(j, 1);
		    // Insert merged list at list position lower in the array (i)
		    neighListList[i] = lij;
		    return neighListList;
		}
	    }
	}
	return neighListList;
    }
    
    // Do two lists have any of the same elements
    checkOverlap(aList, bList) {
	for (const a of aList) {
	    if (bList.includes(a)) {
		return true;
	    }
	}
	return false;
    }
    
    
    writeWithTypes(outFile) {
	if (this.compressed) {
	    console.log('Warning! Cannot write file with types for compressed data.');
	    return;
	}
	    
	let data = '';

	for (let k = 0; k < KIND_NUM; k++) {
	    data += 'BEGIN ' + KIND_NAME[k] + ' types\n'
	    // Loop over the unique parameters
	    for (let pi = 0; pi < this.parameters[k].length; pi++) {
		let typeList = [];
		for (const type of this.unique[k][pi]) typeList.push(type);
		data += pi + ' ' + this.parameters[k][pi].join(' ') + ' ' + typeList + '\n';
	    }
	    data += '\n';
	}

	for (let k = 0; k < KIND_NUM; k++) {
	    data += 'BEGIN ' + KIND_NAME[k] + ' terms\n'
	    // Loop over the interaction terms 
	    for (let ti = 0; ti < this.terms[k].length; ti++) {
		let typeList = [];
		for (const atomIndex of this.terms[k][ti]) typeList.push(this.atoms[atomIndex].type);
		
		// Put the atom indices first and then the  parameter index
		data += this.terms[k][ti].join(' ') + ' '  + this.termParIndices[k][ti] + ' ' + typeList + '\n';
	    }
	    data += '\n';
	}

	fs.writeFileSync(outFile, data);
    }

     writeSimple(outFile) {
	let data = '';

	for (let k = 0; k < KIND_NUM; k++) {
	    data += 'BEGIN ' + KIND_NAME[k] + ' types\n'
	    // Loop over the unique parameters
	    for (let pi = 0; pi < this.parameters[k].length; pi++) {
		data += pi + ' ' + this.parameters[k][pi].join(' ') + '\n';
	    }
	    data += '\n';
	}

	for (let k = 0; k < KIND_NUM; k++) {
	    data += 'BEGIN ' + KIND_NAME[k] + ' terms\n'
	    // Loop over the interaction terms 
	    for (let ti = 0; ti < this.terms[k].length; ti++) {
		let typeList = [];
		for (const atomIndex of this.terms[k][ti]) typeList.push(this.atoms[atomIndex].type);
		
		// Put the atom indices first and then the  parameter index
		data += this.terms[k][ti].join(' ') + ' '  + this.termParIndices[k][ti] + ' ' + typeList + '\n';
	    }
	    data += '\n';
	}

	fs.writeFileSync(outFile, data);
     }

    printSummary() {
	console.log('\nSUMMARY');
	console.log('=======');
	console.log(`Number of atoms: ${this.atoms.length}`);
	console.log(`Number of unique types: ${this.uniqueTypeList.length}`);
	console.log(`Number of unique masses: ${this.uniqueMassList.length}`);
	console.log(`Number of unique charges: ${this.uniqueChargeList.length}`);

	console.log('\nNumber of terms');
	for (let k = 0; k < KIND_NUM; k++) {
	    console.log(`${KIND_NAME[k]} ${this.terms[k].length}`);
	}

	console.log('\nNumber of parameter types');
	for (let k = 0; k < KIND_NUM; k++) {
	    console.log(`${KIND_NAME[k]} types ${this.parameters[k].length}`);
	}
	console.log('=======\n');
    }
    
    radix64(integer) {
	const digits ='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/';
	if (integer == 0) return 0;
	let ret = '';

	// Allow for negative numbers
	if (integer < 0) ret += '-';
	
	
	let value = Math.floor(Math.abs(integer));
	let d = value%64;
	
	while(value > d) {
	    ret += digits[d];
	    value = Math.floor(value/64);
	}

	return ret;
    }

    writeJSFunction(functionName, outFile) {
	console.log('\nWriting structure data');
	const sp = '  ';
	let data = '';

	// Write atom data
	console.log(`Number of atoms: ${this.atoms.length}`);
	data += `function ${functionName}() {\n`
 	data += `${sp}const mol = [];\n\n`;

	data += sp + `mol.name = '${outFile}';\n\n`;

	data += sp + `mol.num = ${this.atoms.length};\n\n`;
	
	data += sp + 'mol.box = ['
	data += this.boxDim.join(',');
	data += '];\n\n';

	// Write the monomer (residue) id numbers, or let the other code figure it out?
	if (this.atomMonomerList != null) {
	    data += sp + 'mol.atom_monomer = [';
	    data += this.atomMonomerList.join(',');
	    data += '];\n\n'
	}

	// Write the fragment id numbers, or let the other code figure it out?
	if (this.atomFragmentList != null) {
	    data += sp + 'mol.atom_fragment = [';
	    data += this.atomFragmentList.join(',');
	    data += '];\n\n'
	}

	let coords = [];
	for (const a of this.atoms) {
	    coords.push(this.fixedDigits(a.x, 5));
	    coords.push(this.fixedDigits(a.y, 5));
	    coords.push(this.fixedDigits(a.z, 5));
	}
	console.log(`Number of coordinates: ${coords.length}`);
	data += sp + 'mol.atom_coord = ['
	data += coords.join(',');
	data += '];\n\n'

	data += `${sp}// TYPE string\n`;
	data += sp + 'mol.par_type = ['
	data += this.uniqueTypeList2.map( item => `'${item}'` ).join(',');
	data += '];\n\n'
	
	data += `${sp}// MASS m(Da)\n`;
	data += sp + 'mol.par_mass = ['
	data += this.uniqueMassList.join(',');
	data += '];\n\n'
	
	data += `${sp}// CHARGE q(e)\n`;
	data += sp + 'mol.par_charge = ['
	data += this.uniqueChargeList.join(',');
	data += '];\n\n'

	data += `${sp}// RESNAME string\n`;
	data += sp + 'mol.par_resname = ['
	data += this.uniqueResNameList.map( item => `'${item}'` ).join(',');
	data += '];\n\n'

	data += `${sp}// NAME string\n`;
	data += sp + 'mol.par_name = ['
	data += this.uniqueNameList.map( item => `'${item}'` ).join(',');
	data += '];\n\n'
	
	// Write the parameter values
	for (let k = 0; k < KIND_NUM; k++) {
	    // Don't write anything for empty parameters
	    if (this.parameters[k] == 0) continue;
	    // Don't bonther with crossterms, which are unsupported at this time
	    if (k == KIND.crossterm) continue;
	    
	    if (k == KIND.lj) {
		data += sp + `// par_${KIND_NAME[k]} fields 2 epsilon(kcal/mol) Rmin(Å) V_LennardJones = epsilon*((Rmin/r)^12 - 2*(Rmin/r)^6 \n`;
		data += sp + '// Note that we use Rmin and not Rmin/2 like CHARMM\n';
	    } else if (k == KIND.bond) {
		data += sp + `// par_${KIND_NAME[k]} fields 2 Kb(kcal/mol) b0(Å) V_bond = Kb*(|ri-rj| - b0)^2 \n`;
	    }  else if (k == KIND.angle) {
		data += sp + `// par_${KIND_NAME[k]} fields 2 Ka(kcal/mol/rad^2) theta0(degrees) V_angle = Ka*(theta - theta0)^2 \n`;
	    }  else if (k == KIND.dihedral) {
		data += sp + `// par_${KIND_NAME[k]} fields 3 Kd(kcal/mol) multiplicity delta(degrees) V_dihedral = Kd*(1 + cos(multiplicity*chi - delta)) \n`;
	    } else if (k == KIND.improper) {
		data += sp + `// par_${KIND_NAME[k]} fields 2 Ki(kcal/mol) psi0(degrees) V_improoper = Ki*(psi - psi0)^2 \n`;
	    } else if (k == KIND.crossterm) {
		data += sp + `// par_${KIND_NAME[k]} fields 1+grid_dimension^2 grid_dimension values(kcal/mol)... \n`;
	    }  else if (k == KIND.bespoke) {
		data += sp + `// par_${KIND_NAME[k]} fields 2 epsilon(kcal/mol) Rmin(Å) V_LennardJones = epsilon*((Rmin/r)^12 - 2*(Rmin/r)^6 \n`;
	    }  else if (k == KIND.exclusion) {
		data += sp + `// par_${KIND_NAME[k]} fields 2 epsilon(kcal/mol) Rmin(Å) V_LennardJones = epsilon*((Rmin/r)^12 - 2*(Rmin/r)^6 \n`;
		data += sp + '// Regular exclusions have epsilon = 0.0 and RMin = -1.0, while bespoke and special 1-4 parameters have nonzero epsilon and positive RMin values\n';
		data += sp + '// Note that we use Rmin and not Rmin/2 like CHARMM\n';
		data += sp + '// Note that this is combination of exclusions and special 1-4 Lennard-Jones parameters\n';
	    }
	    
	    data += sp + `mol.par_${KIND_NAME[k]} = [`;
	    data += this.parameters[k].join(',');
	    data += '];\n\n'
	}

	data += sp + 'mol.index_type = ['
	data += this.atomTypeIndices.join(',');
	data += '];\n\n'

	data += sp + 'mol.index_mass = ['
	data += this.atomMassIndices.join(',');
	data += '];\n\n'
	
	data += sp + 'mol.index_charge = ['
	data += this.atomChargeIndices.join(',');
	data += '];\n\n'

	data += sp + 'mol.index_resname = ['
	data += this.atomResNameIndices.join(',');
	data += '];\n\n'

	data += sp + 'mol.index_name = ['
	data += this.atomNameIndices.join(',');
	data += '];\n\n'
	
	// Write the parameter index for each term
	for (let k = 0; k < KIND_NUM; k++) {
	    // Don't write anything for empty terms
	    if (this.termParIndices[k].length == 0) continue;
	    // Don't bother with crossterms, which are unsupported at this time
	    if (k == KIND.crossterm) continue;
	    
	    data += `${sp}mol.index_${KIND_NAME[k]} = [`;
	    let indices = [];
	    data += this.termParIndices[k].join(',');
	    data += '];\n\n'
	    console.log(`Number of parameter indices for ${KIND_NAME[k]} terms: ${this.termParIndices[k].length}`);
	}

	// Write the atom indices for each term
	for (let k = 0; k < KIND_NUM; k++) {
	    // Don't write anything for empty terms
	    if (this.terms[k] == 0) continue;
	    // Don't bother with crossterms, which are unsupported at this time
	    if (k == KIND.crossterm) continue;
	    
	    if ( k != KIND.lj ) {
		data += `${sp}mol.term_${KIND_NAME[k]} = [`;
		let indices = [];
		for (const t of this.terms[k]) {
		    for (const v of t) indices.push(v);
		}
		data += indices.join(',');
		data += '];\n\n'
		console.log(`Number of atom indices for ${KIND_NAME[k]} terms: ${indices.length}`);
	    }
	}
	data += sp + 'return mol;\n'
	data += '}\n\n'
	
	fs.writeFileSync(outFile, data);
    }

    // Write as a JSON-like file
    write(moleculeName, outFile) {
	console.log('\nWriting structure data');
	const obj = {};

	obj.name = moleculeName;
	obj.file = outFile;
        obj.num = this.atoms.length;
	obj.box = this.boxDim.map(dim => Number(dim.toFixed(6)));

	// Write the monomer (residue) id numbers
	if (this.atomMonomerList != null) {
	    obj.atom_monomer = this.atomMonomerList;
	}

	// Write the fragment id numbers, or let the other code figure it out?
	if (this.atomFragmentList != null) {
	    obj.atom_fragment = this.atomFragmentList;
	}

	let coords = [];
	for (const a of this.atoms) {
	    coords.push(Number(this.fixedDigits(a.x, 5)));
	    coords.push(Number(this.fixedDigits(a.y, 5)));
	    coords.push(Number(this.fixedDigits(a.z, 5)));
	}
	console.log(`Number of coordinates: ${coords.length}`);
	obj.atom_coord = coords;

	// Write the parameter values
	obj.par_type = this.uniqueTypeList2;
	obj.par_mass = this.uniqueMassList;
	obj.par_charge = this.uniqueChargeList;
	obj.par_resname = this.uniqueResNameList;
	obj.par_name = this.uniqueNameList;
	for (let k = 0; k < KIND_NUM; k++) {
	    // Don't write anything for empty parameters
	    if (this.parameters[k] == null || this.parameters[k].length == 0) continue;
	    // Don't bonther with crossterms, which are unsupported at this time
	    if (k == KIND.crossterm) continue;

	    const key = `par_${KIND_NAME[k]}`;
	    //obj[key] = this.parameters[k].flat().map(p => parseFloat(p));
	    // Unfortunately flat() isn't available on nodejs v10 and below
	    const flatParams = [];
	    for (const paramSet of this.parameters[k]) {
		for (const paramValue of paramSet) {
		    flatParams.push(parseFloat(paramValue));
		}
	    }
	    obj[key] = flatParams;
	    console.log(`Number of parameter sets for ${KIND_NAME[k]} terms: ${this.parameters[k].length}`);
	}

	// For atom properties, for each atom
	// we have an index specifying which parameter value to use
	obj.index_type = this.atomTypeIndices;
	obj.index_mass = this.atomMassIndices;
	obj.index_charge = this.atomChargeIndices;
	obj.index_resname = this.atomResNameIndices;
	obj.index_name = this.atomNameIndices;

	// For terms, we write the parameter index for each term
	for (let k = 0; k < KIND_NUM; k++) {
	    // Don't write anything for empty terms
	    if (this.termParIndices[k] == null || this.termParIndices[k].length == 0) continue;
	    // Don't bother with crossterms, which are unsupported at this time
	    if (k == KIND.crossterm) continue;

	    const key = `index_${KIND_NAME[k]}`;
	    obj[key] = this.termParIndices[k];
	    console.log(`Number of parameter indices for ${KIND_NAME[k]} terms: ${this.termParIndices[k].length}`);
	}

	// Now write the terms
	// This consists of atom indices for each term
	for (let k = 0; k < KIND_NUM; k++) {
	    // Don't write anything for empty terms
	    if (this.terms[k] == null || this.terms[k].length == 0) continue;
	    // Don't bother with crossterms, which are unsupported at this time
	    if (k == KIND.crossterm) continue;
	    // Mixing LJ terms are 1-body terms, so no list of terms is needed
	    if (k == KIND.lj) continue;

	    const key = `term_${KIND_NAME[k]}`;
	    let indices = [];
	    for (const t of this.terms[k]) {
		for (const v of t) indices.push(v);
	    }
	    obj[key] = indices;
	    console.log(`Number of atom indices for ${KIND_NAME[k]} terms: ${indices.length}`);
	}
	
	const data = `MOLECULE['${moleculeName}']=${JSON.stringify(obj)};`;
	fs.writeFileSync(outFile, data);
    }
}


if (process.argv.length < 8) {
    console.log(`Usage: ${process.argv[0]} ${process.argv[1]} psfFile pdbFile xscFile prmFile0 [prmFile1]... moleculeName [fragments|monomers] outFile`);
    console.log(`    - The presence of "fragments" means that fragment numbers will be assigned based on the bond network. Note that searching the bond network can be slow.`);
    console.log(`    - The presence of "monomers" means that fragment numbers will be assigned based the segment and resid in the psf.`);
    console.log(`    - With neither fragments nor monomers, no fragment data will be written`);
    process.exit(1);
}

const psf = process.argv[2];
const pdb = process.argv[3];
const xsc = process.argv[4];
let prmList = [];
let moleculeName = '';
let mapFragments = false;
let fragmentType = 'none';
if (process.argv[process.argv.length-2] == 'fragments' || process.argv[process.argv.length-2] == 'monomers') {
    mapFragments = true;
    fragmentType = process.argv[process.argv.length-2];
    prmList = process.argv.slice(5,-3);
    moleculeName = process.argv[process.argv.length-3];
} else {
    prmList = process.argv.slice(5,-2);
    moleculeName = process.argv[process.argv.length-2];
}
const outFile = process.argv[process.argv.length-1];

console.log('');
console.log('psf:', psf);
console.log('pdb:', pdb);
console.log('xsc:', xsc);
console.log('prmList:', prmList);
console.log('moleculeName:', moleculeName);
console.log('mapFragments:', mapFragments);
console.log('fragmentType:', fragmentType);
console.log('outFile:', outFile);

const mol = new MolecularStructure(psf, pdb, xsc, prmList);
mol.printSummary();

// Generate residue monomer ids
mol.generateMonomersByResidue();
// Generate the fragment ids
if (mapFragments) {
    if (fragmentType == 'monomers') {
	// Make the fragments and the monomers the same
	mol.atomFragmentList = [...mol.atomMonomerList];
    } else {
	// Generate the bond graph, which can be slow
	mol.generateFragmentsByBondGraph();
    }
}
mol.compressFragments();

// Write the intermediate structure before compression (so types are still available)
let j = outFile.lastIndexOf('.');
const outFile1 = outFile.slice(0,j) + '.dat';
mol.writeWithTypes(outFile1);

// Merge impropers into dihedrals
mol.mergeImproperToDihedral();

// Merge the bespoke terms
mol.mergeBespokeToExclusion();

// Delete zero parameters
mol.deleteZeroParameters();

// Compress the terms
mol.compress();
mol.printSummary();

// Write the JavaScript function form
// j = outFile.lastIndexOf('.');
// const outFile2 = outFile.slice(0,j) + '.0.js';
// mol.writeJSFunction(`create_molecule_${moleculeName}`, outFile2);

// Write the JSON-like file
mol.write(moleculeName, outFile);
