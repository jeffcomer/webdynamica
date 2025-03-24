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


;(function() {
    // Note that you must set SIM_PARAMETERS and MOLECULE
    'use strict';

    // Enumeration of different simulation states
    const SimState = {none: 0, running: 1, minimizing: 2};
    const SimStateName = ['paused', 'running', 'minimizing'];

    // Default simulation parameters
    const defaultParameters = {
	// Model and physical parameters
	baseMoleculeName: 'default',
	timestep: 1.0, // in femtoseconds
	temperature: 310.0, // in kelvin
        thermostatOn: true, // boolean
	langevinDamping: 0.5, // in ps^-1
	fieldOfView: 25.0, // in degrees
	wallSpring: 2.0, // in kcal mol^-1 Å^-2
	restraintSpring: 10.0, // in kcal mol^-1 Å^-2

	// Interface parameters
	additionalCapacity: 300, // extra atoms in texture (allows adding more molecules)
	initialState: SimState.running, // Running
	
	// Display parameters
	cameraPosFactor: 0.4, // how far is camera in Lz units
	lightPosFactor: 0.7, // how light is camera in Lz units
	fadeFactor: 0.7, // how light is camera in max dimension (max(L)) units
    };

    // Set the default parameters of nothing has been set
    for (const parKey of Object.keys(defaultParameters)) {
	if (SIM_PARAMETERS[parKey] == null) {
	    SIM_PARAMETERS[parKey] = defaultParameters[parKey];
	}
    }
    // Fallback molecule if no molecule is defined
    if (!Object.hasOwn(MOLECULE, 'default')) {
	MOLECULE['default']={"name":"benzene","file":"mol_benzene.webdyn.js","num":12,"box":[30,30,30],"atom_coord":[-1.375,0,-0.002,-2.455,0,-0.001,-0.688,1.191,0.002,-1.227,2.126,-0.001,-0.688,-1.191,-0.001,-1.227,-2.126,0.004,0.688,1.191,-0.004,1.227,2.126,0.005,0.688,-1.191,-0.002,1.227,-2.126,-0.003,1.375,0,-0.001,2.455,0,-0.004],"par_type":["CG2R61","HGR61"],"par_mass":[12.011,1.008],"par_charge":[-0.115,0.115],"par_resname":["BENZ"],"par_name":["CE1","CG","CD1","HD1","CD2","HD2","HG","HE1","CE2","HE2","CZ","HZ"],"par_lj":[0.07,3.9848,0.03,2.7164],"par_bond":[22,2.1525,340,1.08,305,1.375,35,2.4162],"par_angle":[30,120,40,120],"par_dihedral":[4.2,2,180,3.1,2,180,2.4,2,180],"par_exclusion":[0,-1],"index_type":[0,1,0,1,0,1,0,1,0,1,0,1],"index_mass":[0,1,0,1,0,1,0,1,0,1,0,1],"index_charge":[0,1,0,1,0,1,0,1,0,1,0,1],"index_resname":[0,0,0,0,0,0,0,0,0,0,0,0],"index_name":[1,6,2,3,4,5,0,7,8,9,10,11],"index_lj":[0,1,0,1,0,1,0,1,0,1,0,1],"index_bond":[1,2,1,2,1,2,1,2,1,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,3,3],"index_angle":[0,1,0,1,0,1,0,1,0,1,0,0,1,0,0,0,0,0],"index_dihedral":[1,0,1,0,0,2,0,2,1,0,1,0,0,0,2,1,0,1,0,2,0,0,2,2],"index_exclusion":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"term_bond":[0,1,2,0,2,3,4,0,4,5,6,2,6,7,8,4,8,9,10,6,10,8,10,11,0,5,0,3,2,7,2,1,4,9,4,1,6,11,6,3,8,11,8,5,10,9,10,7,0,8,0,6,2,10,2,4,4,10,6,8],"term_angle":[0,4,5,0,4,8,0,2,3,0,2,6,2,6,7,2,6,10,2,0,1,2,0,4,4,8,9,4,8,10,4,0,1,6,10,11,6,10,8,6,2,3,8,10,11,8,4,5,10,8,9,10,6,7],"term_dihedral":[0,4,8,10,0,4,8,9,0,2,6,10,0,2,6,7,1,0,2,6,1,0,2,3,1,0,4,8,1,0,4,5,2,6,10,8,2,6,10,11,2,0,4,8,2,0,4,5,3,2,0,4,3,2,6,10,3,2,6,7,4,8,10,6,4,8,10,11,4,0,2,6,5,4,8,10,5,4,8,9,6,10,8,9,7,6,10,8,7,6,10,11,9,8,10,11],"term_exclusion":[0,1,0,2,0,3,0,4,0,5,0,6,0,8,1,2,1,4,2,3,2,4,2,6,2,7,2,10,3,6,4,5,4,8,4,9,4,10,5,8,6,7,6,8,6,10,6,11,7,10,8,9,8,10,8,11,9,10,10,11]};
    }

    // Make available all molecules loaded into MOLECULE
    // This is used for making molecule thumbnails (selection box) and
    // determining how big the bond, angle, dihedral, exclusion textures need to be.
    SIM_PARAMETERS.moleculeNames = Object.keys(MOLECULE);
    // Converts from resname to a displayed name
    SIM_PARAMETERS.resNameMap = {'AR': 'argon', 'BAMI': 'benzamidinium(+)', 'BENZ': 'benzene', '3CB': 'benzoate(–)', 'C12': 'dodecane', 'BGLC': 'glucose', 'C': 'graphene', 'methane': 'METH', 'N2': 'dinitrogen', 'O2': 'dioxygen', 'GLY': 'glycine', 'PRO': 'proline', 'THR': 'threonine', 'ALA': 'alanine', 'GLU': 'glutamate(–)', 'PHE': 'phenylalanine', 'LYS': 'lysine(+)', 'ARG': 'arginine(+)', 'SER': 'serine', 'VAL': 'valine', 'TRP': 'tryptophan', 'TYR': 'tyrosine', 'TIP3': 'water', 'SPCF': 'water', 'C8': 'octane', 'NAFT': 'naphthalene', 'CLA': 'chloride ion(–)', 'SOD': 'sodium ion(+)', 'POT': 'potassium ion(+)', 'MG': 'magnesium ion(2+)', 'CAL': 'calcium ion(2+)', 'HSO4': 'bisulfate(–)', 'LIG2': 'lignin fragment', 'CELU': 'cellulose', 'QTZ': 'quartz', 'SIO2': 'silicon dioxide nanoparticle', 'ACET': 'acetate(–)', 'ACEH': 'acetic acid', 'IBPD': 'ibuprofen', 'METH': 'methane', 'TCLE': '1,1,1-trichloroethane', 'TBRE': '1,1,1-tribromoethane', '18C6': '18-crown-6'};
    // Converts from resname to an alternate displayed name
    SIM_PARAMETERS.altNameMap = {'C12': 'C12 normal alkane', 'methane': 'CH<sub>4</sub>', 'N2': 'N<sub>2</sub>', 'O2': 'O<sub>2</sub>', 'GLY': 'Gly', 'PRO': 'Pro', 'THR': 'Thr', 'ALA': 'Ala', 'GLU': 'Glu', 'PHE': 'Phe', 'LYS': 'Lys', 'ARG': 'Arg', 'SER': 'Ser', 'VAL': 'Val', 'TRP': 'Trp', 'TYR': 'Tyr', 'TIP3': 'H<sub>2</sub>O', 'BGLC': 'grape sugar, dextrose', 'BENZ': 'C<sub>6</sub>H<sub>6</sub>', 'BAMI': 'benzamidine (protonated)', '3CB': 'benzoic acid (deprotonated)', 'C8': 'C8 normal alkane', 'CLA': 'Cl<sup>–</sup>', 'SOD': 'Na<sup>+</sup>', 'POT': 'K<sup>+</sup>', 'MG': 'Mg<sup>2+</sup>', 'CAL': 'Ca<sup>2+</sup>', 'HSO4': 'HSO<sub>4</sub><sup>–</sup>, hydrogen sulfate', 'CELU': '&beta;(1→4) polymer of glucose', 'QTZ': 'SiO<sub>2</sub>', 'SIO2': 'SiO<sub>2</sub>', '18C6': '18-6 crown ether'};
    // These are the molecule names of residues that can be inserted (mol.name)
    SIM_PARAMETERS.insertNames = ['SPCF','zwit_aa_G','zwit_aa_A','zwit_PRO','zwit_aa_F','zwit_aa_Y','zwit_aa_W','zwit_aa_E','zwit_aa_K','zwit_aa_R','zwit_aa_S','zwit_aa_T','zwit_aa_V','glucose','methane','octane','dodecane','benzene','benzoate','benzamidine','NAFT','N2','O2','argon','sodium_ion','chloride_ion','potassium','magnesium','bisulfate','acetate','acetic_acid','ibuprofen','trichloroethane','tribromoethane','crown'];
    // How to map SVG group ids to molecules 
    // (Assuming the scene svg exists)
    SIM_PARAMETERS.sceneMoleculeMap = {'Air': 'air', 'Seawater': 'seawater', 'Wood': 'wood', 'Sand': 'sand'};


    // Set the default projection and camera
    class RenderInfo {
	constructor(gl, simSys, fieldOfViewDegrees, cameraPosFactor = 0.4, lightPosFactor = 0.7, fadeFactor = 0.8) {
	    // Rendering variables
	    this.shininess = 50.0;
	    this.ambient = 0.20;
	    this.fieldOfView = fieldOfViewDegrees*Math.PI/180.0; // in radians

	    // Find the minimum dimension 
	    const minDim = (simSys.box[0]<simSys.box[1]) ? simSys.box[0] : simSys.box[1];
	    const maxDim = (simSys.box[0]>simSys.box[1]) ? simSys.box[0] : simSys.box[1];
	    this.fadeDistance = fadeFactor*maxDim;
	    
	    this.cameraPos = [0.0, 0.0, cameraPosFactor*maxDim/Math.tan(0.5*this.fieldOfView)];
	    this.lightPos = [0.4*simSys.box[0], 0.45*simSys.box[1], lightPosFactor*simSys.box[2]];
	    this.targetPos = [0.0, 0.0, 0.0];
	    this.up = [0, 1, 0];
	    this.aspect = gl.canvas.clientWidth / gl.canvas.clientHeight;

	    // View
	    this.near = simSys.radius;
	    const maxDimZ = (simSys.box[2]>maxDim) ? simSys.box[2] : maxDim;
	    this.far = 5.0*maxDimZ;
	    this.cameraMatrix = m4.lookAt(this.cameraPos, this.targetPos, this.up);
	    this.viewMatrix = m4.inverse(this.cameraMatrix);
	    this.projectionMatrix = m4.perspective(this.fieldOfView, this.aspect, this.near, this.far);
	    this.viewProjectionMatrix = m4.multiply(this.projectionMatrix, this.viewMatrix);
	}
    }

    class MDTextures {
	constructor(gl, atomTexData, shaders, simSys) {
	    // Store the texture dimensions here also
	    this.width = atomTexData.width;
	    this.height = atomTexData.height;
	    this.size = this.width*this.height;
	    this.atomDim = {width: this.width, height: this.height, size: this.size};
	    this.bondDim = atomTexData.bondDim;
	    this.angleDim = atomTexData.angleDim;
	    this.dihedralDim = atomTexData.dihedralDim;
	    this.excludeDim = atomTexData.excludeDim;
	    this.hashDim = atomTexData.hashDim;

	    // Hash function (linear congruential generator) multipliers
	    // Note that this is stored as a reference, so it
	    // changes when the hash table is regenerated
	    this.hashA = atomTexData.hashA;
	    
	    // Force field textures
	    this.bond = createTexture(gl, atomTexData.bond, atomTexData.bondDim);
	    this.angle = createTexture(gl, atomTexData.angle, atomTexData.angleDim);
	    this.dihedral = createTexture(gl, atomTexData.dihedral, atomTexData.dihedralDim);
	    this.exclude = createTexture(gl, atomTexData.exclude, atomTexData.excludeDim);
	    this.nonbond = createTexture(gl, atomTexData.nonbond, atomTexData);
	    // Hash table for exclusions and special L-J interactions
	    this.hash = createTexture(gl, atomTexData.hash, atomTexData.hashDim);

	    // Restraint texture
	    this.restrain = createTexture(gl, atomTexData.restrain, atomTexData);

	    // Color and radius texture for rendering
	    this.color = createTexture(gl, atomTexData.color, atomTexData);

	    // Selection texture [fragment, monomer, material, radius]
	    this.select = createTexture(gl, atomTexData.select, atomTexData);

	    // Initial random numbers for random number textures
	    const randomTexData0 = new Float32Array(atomTexData.nonbond.length);
	    const randomTexData1 = new Float32Array(atomTexData.nonbond.length);
	    // Must be between 1 and RAND_M-1 (zero is bad)
	    // Choosing the m values correctly ensures a good distribution for the first step
	    for (let i = 0; i < randomTexData0.length; i++) {
		let m = shaders.RAND_M3;
		switch(i%4) {
		case 0:
		    m = shaders.RAND_M0;
		    break;
		case 1:
		    m = shaders.RAND_M1;
		    break;
		case 2:
		    m = shaders.RAND_M2;
		    break;
		case 3:
		    m = shaders.RAND_M3;
		    break;
		}

		// For each generator, the possible values go from 1 to m-1
		// (x mod m) == 0 will only produce 0 and must be avoided
		randomTexData0[i] = Math.floor(Math.random()*(m-2.0)) + 1;
		randomTexData1[i] = Math.floor(Math.random()*(m-2.0)) + 1;
	    }

	    // Position textures
	    const zeroData = new Float32Array(atomTexData.pos.length);
	    //for (let i = 0; i < zeroData.length; i++) zeroData[i] = (i%4==0)?1.0:0.0; 
	    this.pos0 = createTexture(gl, atomTexData.pos, atomTexData);
	    this.pos1 = createTexture(gl, simSys.pos[1], atomTexData);
	    this.extraPos = createTexture(gl, simSys.pos[0], atomTexData);
	    this.vel = createTexture(gl, simSys.vel[0], atomTexData);
	    this.halfVel = createTexture(gl, zeroData, atomTexData); // v(t + dt/2)
	    this.force = createTexture(gl, zeroData, atomTexData);
	    this.extraForce = createTexture(gl, zeroData, atomTexData);
	    this.random0 = createTexture(gl, randomTexData0, atomTexData);
	    this.random1 = createTexture(gl, randomTexData1, atomTexData);

	    // Create the frame buffers
	    this.posFB0 = createFramebuffer(gl, this.pos0);
	    this.posFB1 = createFramebuffer(gl, this.pos1);
	    this.extraPosFB = createFramebuffer(gl, this.extraPos);
	    this.velFB = createFramebuffer(gl, this.vel);
	    this.halfVelFB = createFramebuffer(gl, this.halfVel);
	    this.forceFB = createFramebuffer(gl, this.force);
	    this.extraForceFB = createFramebuffer(gl, this.force);
	    this.randomFB0 = createFramebuffer(gl, this.random0);
	    this.randomFB1 = createFramebuffer(gl, this.random1);

	    // Picking molecule texture
	    this.pickDim = {width: 1, height: 1};
	    this.pick = createTexture(gl, null, this.pickDim);
	    const depthBuffer = gl.createRenderbuffer();
	    gl.bindRenderbuffer(gl.RENDERBUFFER, depthBuffer);
	    gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, 1, 1);

	    // Create and bind the framebuffer for the picking texture
	    this.pickFB = gl.createFramebuffer();
	    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickFB);
	    // Make a 1 by 1 depth buffer
	    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.pick, 0);
	    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, depthBuffer);

	    // Render to this unit square quad
	    this.unitSquareBuffer = gl.createBuffer();
	    gl.bindBuffer(gl.ARRAY_BUFFER, this.unitSquareBuffer);
	    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1]), gl.STATIC_DRAW);

	    // Vertex buffer for the half-sphere model for rendering atoms
	    this.sphereIdBuffer = gl.createBuffer();
	    gl.bindBuffer(gl.ARRAY_BUFFER, this.sphereIdBuffer);
	    const sphereId = new Float32Array([...Array(shaders.SPHERE_VERTICES*atomTexData.size).keys()]);
	    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(sphereId), gl.STATIC_DRAW);
	}
    }
    
    function createTexture(gl, data, dim) {
	const tex = gl.createTexture();
	gl.bindTexture(gl.TEXTURE_2D, tex);
	gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, dim.width, dim.height, 0, gl.RGBA, gl.FLOAT, data);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
	return tex;
    }
    
    function createFramebuffer(gl, tex) {
	const fb = gl.createFramebuffer();
	gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
	gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
	return fb;
    }

    function prepareWebGL() {
	// Get a WebGL context
	const canvas = document.getElementById("canvas3d");
	const gl = canvas.getContext("webgl");
	if (!gl) {
	    alert('WebGL unavailable');
	}
	// Floating point textures?
	const ext1 = gl.getExtension('OES_texture_float');
	if (!ext1) {
	    alert('Need OES_texture_float');
	}
	// Render to floating point textures?
	const ext2 = gl.getExtension('WEBGL_color_buffer_float');
	if (!ext2) {
	    alert('Need WEBGL_color_buffer_float');
	}
	// Use textures in a vertex shader?
	if (gl.getParameter(gl.MAX_VERTEX_TEXTURE_IMAGE_UNITS) < 1) {
	    alert('Can not use textures in vertex shaders');
	}
	return gl;
    }

    // Find the maximum number of each term per atom
    // Generate the atom texture data structure
    // Returns the [atomTexData, simSys]
    // atomTexData is the AtomTextureData object and simSys is the SimulationSystem object
    function prepareAtoms(MOLECULE, simParams, baseMoleculeName) {
	// Check that the base molecule has been loaded
	if (!Object.hasOwn(MOLECULE, baseMoleculeName)) {
	    console.log(`ERROR! baseMolecule ${baseMoleculeName} has not been defined. Make sure the name is correct and that the script containing this molecule has been loaded in the HTML file.`)
	    return [];
	}
	const baseMolecule = MOLECULE[baseMoleculeName];
	
	// Make the molecule thumbnails for the selection panel
	console.log(`Base molecule "${baseMolecule.name}" has ${baseMolecule.num} atoms`);
	const thumbnails = {};
	for (const name of simParams.moleculeNames) {
	    if (!Object.hasOwn(MOLECULE, name)) {
		console.log(`Warning! Molecule "${name}" from moleculeNames has not been defined. Make sure the name is correct and that the script containing this molecule has been loaded in the HTML file.`)
	    } else {
		const mol = MOLECULE[name];
		// If the molecule has only one residue name, we will make a thumbnail associated with it
		if (mol.par_resname != null && mol.par_resname.length == 1) {
		    const resName = mol.par_resname[0];
		    thumbnails[resName] = new MoleculeThumbnail(mol, name);
		    //console.log('thumbnail', resName);
		}
		//console.log(`Molecule "${mol.name}" has ${mol.num} atoms`);
	    }
	}

	// Make the molecule thumbnails for the insertion panel
	const insertThumbnails = [];
	for (const name of simParams.insertNames) {
	    //if (!Object.hasOwn(MOLECULE, name)) {
		//console.log(`Warning! Molecule "${name}" from insertNames has not been defined. Make sure the name is correct and that the script containing this molecule has been loaded in the HTML file.`)
	    //}

	    // Silently ignore molecules that haven't been loaded.
	    if (Object.hasOwn(MOLECULE, name)) {
		const mol = MOLECULE[name];
		let resName = 'UNK';
		let moleculeName = name;
		if (mol.par_resname != null) resName = mol.par_resname[0];
		if (Object.hasOwn(simParams.resNameMap, resName)) {
		    // Choose a better displayed name
		    moleculeName = simParams.resNameMap[resName];
		}
		insertThumbnails.push(new MoleculeThumbnail(mol, moleculeName));
	    }
	}
	
	// Get the maximum number of terms per atom for each molecule
	const termsPerAtom = [];
	for (let k = 0; k <= TermInfo.exclude; k++) {
	    termsPerAtom[k] = AtomTextureData.maxTermsPerAtom(baseMolecule, k);
	}
	for (let k = 0; k <= TermInfo.exclude; k++) {
	    for (const name of simParams.moleculeNames) {
		// Check that this molecule exists
		if (!Object.hasOwn(MOLECULE, name)) {
		    console.log(`Warning! The molecule ${name} has not been defined and will not be loaded.`)
		} else {
		    const maxTerms = AtomTextureData.maxTermsPerAtom(MOLECULE[name], k);
		    if (maxTerms > termsPerAtom[k]) termsPerAtom[k] = maxTerms;
		}
	    }
	    console.log(`${TermInfo.name[k]} terms per atom: ${termsPerAtom[k]}`);
	}
	
	// Create an object with the atom texture data
	const atomTexData = new AtomTextureData(baseMolecule, termsPerAtom, simParams.additionalCapacity);

	// Construct the simulation system
	const simSys = new SimulationSystem(atomTexData.pos, atomTexData.nonbond, baseMolecule.box, simParams.timestep);
	if (simParams.wallSpring != null) simSys.wallSpring = simParams.wallSpring;
	if (simParams.restraintSpring != null) simSys.restraintSpring = simParams.restraintSpring;
	if (simParams.wall != null) simSys.wall = simParams.wall;
	simSys.setLangevin(simParams.temperature, simParams.langevinDamping);
	console.log('box size:', simSys.box[0], simSys.box[1], simSys.box[2]);
	console.log('wall position:', simSys.wall[0], simSys.wall[1], simSys.wall[2]);
	console.log('temperature:', simSys.temper);
	console.log('timestep in simulation units:', simSys.dt);
	console.log('langevinDamping in simulation units:', simSys.langevinGamma);
	console.log('langevinRatio (0.5*langevinGamma*dt):', simSys.langevinRatio);
	console.log('langevinRandom (2.0*langevinGamma*kT/dt):', simSys.langevinRandom);
	
	// Apply restraints to graphene
	const restrainNum = atomTexData.restrainBackground(simSys.restraintSpring);
	console.log(`Restraining ${restrainNum} atoms`);
	
	return [atomTexData, simSys, thumbnails, insertThumbnails];
    }

    
    function calcKineticEnergy(currPos, currVel, currNonbond, simSys) {
	const n = Math.floor(currPos.length/4);
	let activeNum = 0;
	let enerKin = 0.0;
	
	for (let ai = 0; ai < n; ai++) {
	    // Ignore inactive atoms
	    if (currPos[4*ai + 3] == 0.0) continue;
	    activeNum++;
	    
	    let sumSq = 0.0;
	    for (let c = 0; c < 3; c++) {
		// Simulation vel. units: Å/sqrt(Da Å^2/(kcal_mol)) = sqrt(kcal_mol/(Da Å))
		sumSq += currVel[4*ai + c]**2;
	    }
	    enerKin += 0.5*currNonbond[4*ai]*sumSq;
	}
	const temperKin = 2.0*enerKin/(3.0*simSys.boltz*activeNum);
	return [enerKin, temperKin];
    }

    function calcPotentialEnergy(currPos, currForce) {
	const n = Math.floor(currPos.length/4);
	let enerPot = 0.0;
	for (let i = 3; i < currForce.length; i+=4) {
	    if (currPos[i] == 0.0) continue;
	    enerPot += currForce[i];
	}
	return enerPot;
    }

    //////////////////////////////////////////////////////////////////////
    // User interface logic
    //////////////////////////////////////////////////////////////////////
    function getMouseMatrix(gl, mouseX, mouseY, renderInfo) {
	// Render to the canvas
	// compute the rectangle the near plane of our frustum covers
	const aspect = gl.canvas.clientWidth / gl.canvas.clientHeight;
	const top = Math.tan(renderInfo.fieldOfView * 0.5) * renderInfo.near;
	const bottom = -top;
	const left = aspect * bottom;
	const right = aspect * top;
	const width = Math.abs(right - left);
	const height = Math.abs(top - bottom);

	// Compute the portion of the near plane covers the 1 pixel under the mouse
	const pixelX = mouseX * gl.canvas.width / gl.canvas.clientWidth;
	const pixelY = gl.canvas.height - mouseY * gl.canvas.height / gl.canvas.clientHeight - 1;
	
	const subLeft = left + pixelX * width / gl.canvas.width;
	const subBottom = bottom + pixelY * height / gl.canvas.height;
	const subWidth = width / gl.canvas.width;
	const subHeight = height / gl.canvas.height;
	
	// Make a frustum for that 1 pixel
	const projectionMatrix = m4.frustum(subLeft, subLeft + subWidth, subBottom, subBottom + subHeight, renderInfo.near, renderInfo.far);
	const viewProjectionMatrix = m4.multiply(projectionMatrix, renderInfo.viewMatrix);    
	
	return viewProjectionMatrix; 
    }

    // Structure for following the mouse
    class MouseInfo {
	constructor(mouseScale = 50.0) {
	    this.x = -1; // mouse x coordinate in the canvas in pixels
	    this.y = -1; // mouse y coordinate in the canvas in pixels
	    this.atom = -1; // atom index
	    this.monomer = -1; // monomer (residue) index
	    this.material = -1; // type of material (0: normal, 1: water, 2: background)
	    this.fragment = -1; // fragment index
	    this.type = 'UNK'; // atom type
	    this.resName = 'none'; // residue name
	    this.moleculeName = 'none'; // molecule name
	    this.worldX = 0;
	    this.worldY = 0;
	    this.worldZ = 0;
	    this.scale = mouseScale;
            this.atomX = 0;
	    this.atomY = 0;
	    this.atomZ = 0;
	    this.down = false;
	    this.dragging = false;
	    this.dragOriginX = -1;
	    this.dragOriginY = -1;
	    this.minMove = 4; 
	}

	// Attempt to calculate the xy-position of the mouse pointer in the simulation world 
	// There is something wrong with the logic, since we don't get the correct positions
	// and have to use a fudge factor (this.scale)
	toWorld(gl, viewProjectionMatrix, overrideMouseX, overrideMouseY) {
	    const mouseX = overrideMouseX ?? this.x;
	    const mouseY = overrideMouseY ?? this.y;
	    
	    const M = viewProjectionMatrix;
	    // Convert to unit square coordinates
	    const mx = (2.0*mouseX - gl.canvas.clientWidth)/gl.canvas.clientWidth;
	    const my = -(2.0*mouseY - gl.canvas.clientHeight)/gl.canvas.clientHeight;

	    const worldX = (M[13]*M[4] - M[12]*M[5] + (M[4]*M[9] - M[5]*M[8])*this.atomZ + M[5]*mx - M[4]*my)/(M[0]*M[5] - M[1]*M[4]);
	    const worldY = (M[1]*M[12] - M[0]*M[13] + (M[1]*M[8] - M[0]*M[9])*this.atomZ + M[0]*my - M[1]*mx)/(M[0]*M[5] - M[1]*M[4]);

	    this.worldX = worldX;
	    this.worldY = worldY;
	    this.worldZ = this.atomZ;
	}
    }


    
    //////////////////////////////////////////////////////////////////////
    // Functions that run shaders
    //////////////////////////////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////////////////
    // Draw scene shader
    function drawAtoms(env, program, progLocs, posTex, selectInfo, renderInfo, frameBuffer, overrideWidth, overrideHeight) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object

	const canvasWidth = overrideWidth || gl.canvas.width;
	const canvasHeight = overrideHeight || gl.canvas.height;
	
	// Check for resizing of the canvas
	webglUtils.resizeCanvasToDisplaySize(gl.canvas);
	gl.viewport(0, 0, canvasWidth, canvasHeight);

	// Write to the canvas or a framebuffer (for framebuffer = null)
	gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuffer);
	
	// The sphere vertex id buffer 
	gl.bindBuffer(gl.ARRAY_BUFFER, env.textures.sphereIdBuffer);
	gl.enableVertexAttribArray(progLocs.id);
	gl.vertexAttribPointer(progLocs.id,1,gl.FLOAT,false,0,0);
	
	// Clear the canvas and depth buffer
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	// Turn on backface culling
	gl.enable(gl.CULL_FACE);
	// Enable the depth buffer
	gl.enable(gl.DEPTH_TEST);

	// Select the textures
	const colorTexIndex = 0;
	const posTexIndex = 1;
	const selectTexIndex = 2;
	gl.activeTexture(gl.TEXTURE0 + colorTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.color);
	gl.activeTexture(gl.TEXTURE0 + posTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, posTex);
	gl.activeTexture(gl.TEXTURE0 + selectTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.select);
	
	/////////////////////////
	// Draw the scene
	gl.useProgram(program);
	gl.uniform1i(progLocs.posTex, posTexIndex);
	gl.uniform1i(progLocs.colorTex, colorTexIndex);
	gl.uniform1i(progLocs.selectTex, selectTexIndex);
	gl.uniform1f(progLocs.selectId, selectInfo.id);
	gl.uniform4f(progLocs.selectMask, ...selectInfo.mask);
	gl.uniform1f(progLocs.selectScale, selectInfo.scale);
	gl.uniform1f(progLocs.hideMaterial, selectInfo.hideMaterial);
	gl.uniformMatrix4fv(progLocs.matrix, false, renderInfo.viewProjectionMatrix);
	gl.uniform1f(progLocs.shininess, renderInfo.shininess);
	gl.uniform1f(progLocs.ambient, renderInfo.ambient);
	gl.uniform1f(progLocs.fadeDistance, renderInfo.fadeDistance);
	gl.uniform3fv(progLocs.cameraPos, renderInfo.cameraPos);
	gl.uniform3fv(progLocs.lightPos, renderInfo.lightPos);
	gl.drawArrays(gl.TRIANGLES, 0, env.shaders.SPHERE_VERTICES*textures.size);

	// These settings can interfere with data extraction
	gl.disable(gl.CULL_FACE);
	gl.disable(gl.DEPTH_TEST);
    }

    // Calculate the atomic forces using a shader
    function runForceShader(env, posTex, destForceFrameBuffer, fast = false) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object

	// Render to the next force
	gl.bindFramebuffer(gl.FRAMEBUFFER, destForceFrameBuffer);
	// The size of our destination texture
	gl.viewport(0, 0,  textures.width, textures.height);

	const bondIndex = 0;
	const angleIndex = 1;
	const dihedralIndex = 2;
	const excludeIndex = 3;
	const nonbondIndex = 4;
	const restrainIndex = 5;
	const posIndex = 6;
	const hashIndex = 7;
	
	// Bind the textures
	gl.activeTexture(gl.TEXTURE0 + bondIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.bond);
	gl.activeTexture(gl.TEXTURE0 + angleIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.angle);
	gl.activeTexture(gl.TEXTURE0 + dihedralIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.dihedral);
	gl.activeTexture(gl.TEXTURE0 + excludeIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.exclude);
	gl.activeTexture(gl.TEXTURE0 + nonbondIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.nonbond);
	gl.activeTexture(gl.TEXTURE0 + restrainIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.restrain);
	gl.activeTexture(gl.TEXTURE0 + posIndex);
	gl.bindTexture(gl.TEXTURE_2D, posTex);
	gl.activeTexture(gl.TEXTURE0 + hashIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.hash);

	// Do we use the fast shader with near cutoff?
	let locs = null;
	let prog = null;
	if (fast) {
	    prog = shaders.forceFast;
	    locs = shaders.forceFastLocs;
	} else {
	    prog = shaders.force;
	    locs = shaders.forceLocs;
	}
	
	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(locs.dummy);
	gl.vertexAttribPointer(locs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the force calculation
	gl.useProgram(prog);
	gl.uniform1i(locs.bondTex, bondIndex);
	gl.uniform1i(locs.angleTex, angleIndex);
	gl.uniform1i(locs.dihedralTex, dihedralIndex);
	gl.uniform1i(locs.excludeTex, excludeIndex);
	gl.uniform1i(locs.nonbondTex, nonbondIndex);
	gl.uniform1i(locs.posTex, posIndex);
	gl.uniform1i(locs.hashTex, hashIndex);
	gl.uniform1f(locs.hashA0, textures.hashA[0]);
	gl.uniform1f(locs.hashA1, textures.hashA[1]);
	gl.uniform1i(locs.restrainTex, restrainIndex);
	gl.uniform3f(locs.box, ...env.simSys.box);
	gl.uniform3f(locs.wall, ...env.simSys.wall);
	gl.uniform1f(locs.wallSpring, env.simSys.wallSpring);
	gl.uniform1f(locs.coulomb, env.simSys.coulomb);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    // Update velocities using one of the "kick" shaders
    // This can be the first or second kick of the NVE velocity-Verlet
    // or Langevin velocity-Verlet shaders
    function runUpdateVelocityShader(env, program, progLocs, srcVelTex, forceTex, randomTex, destVelFB) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object

	// Render to the destination
	gl.bindFramebuffer(gl.FRAMEBUFFER, destVelFB);
	// The size of our destination texture
	gl.viewport(0, 0,  textures.width,  textures.height);
	
	// Bind the textures
	const srcVelIndex = 0;
	const forceIndex = 1;
	const randomIndex = 2;
	const nonbondIndex = 4;
	gl.activeTexture(gl.TEXTURE0 + srcVelIndex);
	gl.bindTexture(gl.TEXTURE_2D, srcVelTex);
	gl.activeTexture(gl.TEXTURE0 + forceIndex);
	gl.bindTexture(gl.TEXTURE_2D, forceTex);
	gl.activeTexture(gl.TEXTURE0 + randomIndex);
	gl.bindTexture(gl.TEXTURE_2D, randomTex);
	gl.activeTexture(gl.TEXTURE0 + nonbondIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.nonbond);

	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(progLocs.dummy);
	gl.vertexAttribPointer(progLocs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the Verlet "kick" shader
	gl.useProgram(program);
	gl.uniform1i(progLocs.velTex, srcVelIndex);
	gl.uniform1i(progLocs.forceTex, forceIndex);
	gl.uniform1i(progLocs.randomTex, randomIndex);
	gl.uniform1i(progLocs.nonbondTex, nonbondIndex);
	gl.uniform1f(progLocs.dt, env.simSys.dt);
	if (progLocs.langevinRatio != null) {
	    // The NVE shader doesn't have these
	    gl.uniform1f(progLocs.langevinRatio, env.simSys.langevinRatio);
	    gl.uniform1f(progLocs.langevinRandom, env.simSys.langevinRandom);
	}
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    // Update positions using on of the Velocity Verlet "drift" shader
    function runUpdatePositionShader(env, posTex, velTex, destPosFB) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object

	// Render to the next positions: r(t+dt)
	gl.bindFramebuffer(gl.FRAMEBUFFER, destPosFB);
	// The size of our destination texture
	gl.viewport(0, 0, textures.width, textures.height);
	
	// Bind the textures
	const nonbondIndex = 4;
	const velIndex = 5;
	const posIndex = 6;
	gl.activeTexture(gl.TEXTURE0 + posIndex);
	gl.bindTexture(gl.TEXTURE_2D, posTex);
	gl.activeTexture(gl.TEXTURE0 + velIndex);
	gl.bindTexture(gl.TEXTURE_2D, velTex);
	gl.activeTexture(gl.TEXTURE0 + nonbondIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.nonbond);

	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.driftLocs.dummy);
	gl.vertexAttribPointer(shaders.driftLocs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the Verlet "drift" shader
	gl.useProgram(shaders.drift);
	gl.uniform1i(shaders.driftLocs.posTex, posIndex);
	gl.uniform1i(shaders.driftLocs.velTex, velIndex);
	gl.uniform1i(shaders.driftLocs.nonbondTex, nonbondIndex);
	gl.uniform3f(shaders.driftLocs.box, ...env.simSys.box);
	gl.uniform1f(shaders.driftLocs.dt, env.simSys.dt);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    
    // Update the random numbers
    // Implicit parameters: atomTexInfo, unitSquareBuffer
    function runUpdateRandomShader(env, srcTex, destFrameBuffer) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object
	
	// Render to the next random number texture
	gl.bindFramebuffer(gl.FRAMEBUFFER, destFrameBuffer);
	// The size of our destination texture
	gl.viewport(0, 0, textures.width, textures.height);
	
	// Bind the textures
	const randomIndex = 7
	gl.activeTexture(gl.TEXTURE0 + randomIndex);
	gl.bindTexture(gl.TEXTURE_2D, srcTex);

	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.randomLocs.dummy);
	gl.vertexAttribPointer(shaders.randomLocs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the Verlet step shader
	gl.useProgram(shaders.random);
	gl.uniform1i(shaders.randomLocs.randomTex, randomIndex);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    // Displace selected atoms so that the picked atom is at the world position of the mouse
    // The world position of the mouse is the mouse position on the 2D rendered projection
    // transformed back into the three-dimensional simulation space.
    // Since the mouse position only gives x and y, the z position is left free using axisMask
    // Displace atoms using the x channel from the selectTex texture to mark atoms
    // Typically we'll select using the fragment texture
    // Implicit parameters: atomTexInfo, unitSquareBuffer 
    function runDisplaceRelativeShader(env, posTex, selectInfo, desiredPos, destPosFrameBuffer) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object
	
	// Displace only in x and y, leaving z unchanged
	const axisMask = [1, 1, 0];
	
	// Render to the next random number texture
	gl.bindFramebuffer(gl.FRAMEBUFFER, destPosFrameBuffer);
	// The size of our destination texture
	gl.viewport(0, 0, textures.width, textures.height);
		
	// Bind the textures
	const posTexIndex = 1;
	const selectTexIndex = 2;
	gl.activeTexture(gl.TEXTURE0 + posTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, posTex);
	gl.activeTexture(gl.TEXTURE0 + selectTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.select);

	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.displaceRelativeLocs.dummy);
	gl.vertexAttribPointer(shaders.displaceRelativeLocs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the velocity reset shader
	gl.useProgram(shaders.displaceRelative);
	gl.uniform1i(shaders.displaceRelativeLocs.posTex, posTexIndex);
	gl.uniform1i(shaders.displaceRelativeLocs.selectTex, selectTexIndex);
	gl.uniform1f(shaders.displaceRelativeLocs.selectId, selectInfo.id);
	gl.uniform4f(shaders.displaceRelativeLocs.selectMask, ...selectInfo.mask);
	gl.uniform1f(shaders.displaceRelativeLocs.pickedAtomIndex, selectInfo.atom);
	gl.uniform3f(shaders.displaceRelativeLocs.desiredPos, ...desiredPos);
	gl.uniform3f(shaders.displaceRelativeLocs.axisMask, ...axisMask);
	gl.uniform3f(shaders.displaceRelativeLocs.box, ...env.simSys.box);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }
    
    // Use a shader to reset the velocities according to a Maxwell-Boltzmann distribution
    // at the given temperature
    // Implicit parameters: atomTexInfo, unitSquareBuffer
    function runResetVelocityShader(env, randomTex, kT, destVelFrameBuffer) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object
	
	// Render to the next random number texture
	gl.bindFramebuffer(gl.FRAMEBUFFER, destVelFrameBuffer);
	// The size of our destination texture
	gl.viewport(0, 0,  textures.width, textures.height);
	
	// Bind the textures
	const nonbondTexIndex = 4
	const randomTexIndex = 7
	gl.activeTexture(gl.TEXTURE0 + nonbondTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.nonbond);
	gl.activeTexture(gl.TEXTURE0 + randomTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, randomTex);

	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.initVelLocs.dummy);
	gl.vertexAttribPointer(shaders.initVelLocs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the velocity reset shader
	gl.useProgram(shaders.initVel);
	gl.uniform1i(shaders.initVelLocs.nonbondTex, nonbondTexIndex);
	gl.uniform1i(shaders.initVelLocs.randomTex, randomTexIndex);
	gl.uniform1f(shaders.initVelLocs.kT, kT);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    // Move the atoms downhill
    function runDownhillShader(env, posTex, forceTex, stepSize, destPosFrameBuffer) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object
	
	// Render to the next random number texture
	gl.bindFramebuffer(gl.FRAMEBUFFER, destPosFrameBuffer);
	// The size of our destination texture
	gl.viewport(0, 0, textures.width, textures.height);

	// Bind the textures
	const posTexIndex = 6;
	const forceTexIndex = 0;
	gl.activeTexture(gl.TEXTURE0 + posTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, posTex);
	gl.activeTexture(gl.TEXTURE0 + forceTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, forceTex);

	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.downhillLocs.dummy);
	gl.vertexAttribPointer(shaders.downhillLocs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the shader program
	gl.useProgram(shaders.downhill);
	gl.uniform1i(shaders.downhillLocs.posTex, posTexIndex);
	gl.uniform1i(shaders.downhillLocs.forceTex, forceTexIndex);
	gl.uniform1f(shaders.downhillLocs.stepSize, stepSize);
	gl.uniform3f(shaders.downhillLocs.box, ...env.simSys.box);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    // Displace the selection
    function runDisplaceShader(env, srcPosTex, selectInfo, displacement, destPosFrameBuffer) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object
	
	// Render to the next random number texture
	gl.bindFramebuffer(gl.FRAMEBUFFER, destPosFrameBuffer);
	// The size of our destination texture
	gl.viewport(0, 0, textures.width, textures.height);

	// Bind the textures
	const posTexIndex = 6;
	const selectTexIndex = 2;
	gl.activeTexture(gl.TEXTURE0 + posTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, srcPosTex);
	gl.activeTexture(gl.TEXTURE0 + selectTexIndex);
	gl.bindTexture(gl.TEXTURE_2D, textures.select);
	
	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.displaceLocs.dummy);
	gl.vertexAttribPointer(shaders.displaceLocs.dummy, 2, gl.FLOAT, false, 0, 0);

	// Run the shader program
	gl.useProgram(shaders.displace);
	gl.uniform1i(shaders.displaceLocs.posTex, posTexIndex);
	gl.uniform1i(shaders.displaceLocs.selectTex, selectTexIndex);
	gl.uniform1f(shaders.displaceLocs.selectId, selectInfo.id);
	gl.uniform4f(shaders.displaceLocs.selectMask, ...selectInfo.mask);
	gl.uniform3f(shaders.displaceLocs.displace, ...displacement);
	gl.uniform3f(shaders.displaceLocs.box, ...env.simSys.box);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    
    // Just copy a texture
    function runCopyShader(env, srcTex, destFrameBuffer) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object
	
	// Render to the next random number texture
	gl.bindFramebuffer(gl.FRAMEBUFFER, destFrameBuffer);
	// The size of our destination texture
	gl.viewport(0, 0, textures.width, textures.height);

	// Bind the textures
	const texIndex = 6;
	gl.activeTexture(gl.TEXTURE0 + texIndex);
	gl.bindTexture(gl.TEXTURE_2D, srcTex);

	// Draw to the unit square 
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.copyLocs.dummy);
	gl.vertexAttribPointer(shaders.copyLocs.dummy, 2, gl.FLOAT, false, 0, 0);
	
	// Run the shader program
	gl.useProgram(shaders.copy);
	gl.uniform1i(shaders.copyLocs.posTex, texIndex);
	gl.drawArrays(gl.TRIANGLES, 0, 6);  // Draw 2 triangles (6 vertices)
    }

    
    /////////////////////////////////////////
    // Extract data from the texture
    function extractData(env, texture, width, height) {
	const gl = env.gl; // WebGL context
	const textures = env.textures; // MDTextures object
	const shaders = env.shaders; // MDShaders object
	
	// WebGL1 only allows to render to a 4-byte (RGBA) frame buffer
	// Destination texture has to be 4 times bigger to hold x, y, z, w components of four (4-byte) floats
	const destWidth = 4*width;
	const destHeight = height;
	
	// Render to this destination canvas
	gl.canvas.width = destWidth;
	gl.canvas.height = destHeight;
	gl.viewport(0, 0, destWidth, destHeight); // This is necessary!
	// Render to the canvas
	gl.bindFramebuffer(gl.FRAMEBUFFER, null);

	// Activate the texture
	const texIndex = 0;
	gl.activeTexture(gl.TEXTURE0 + texIndex);
	gl.bindTexture(gl.TEXTURE_2D, texture);

	// Render from the second texture
	// Run the second shader
	gl.useProgram(shaders.extract);
	gl.bindBuffer(gl.ARRAY_BUFFER, textures.unitSquareBuffer);
	gl.enableVertexAttribArray(shaders.extractLocs.dummy);
	gl.vertexAttribPointer(shaders.extractLocs.dummy, 2, gl.FLOAT, false, 0, 0);
	gl.uniform1i(shaders.extractLocs.dataTex, texIndex); 
	gl.uniform2f(shaders.extractLocs.destTexDim, destWidth, destHeight);
	gl.drawArrays(gl.TRIANGLES, 0, 6);

	// Get the raw values from the texture
	// each pixel (RGBA) becomes a single 32-bit float (thus, multiply by 4)
	const raw = new Uint8Array(4 * destWidth * destHeight);
	gl.readPixels(0, 0, destWidth, destHeight, gl.RGBA, gl.UNSIGNED_BYTE, raw);
	const dataArray = new Float32Array(raw.buffer);

	// Revert to the correct size
	webglUtils.resizeCanvasToDisplaySize(gl.canvas);
	gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
	
	return dataArray;
    }


    function setSceneColors(sceneElementColors, highlightGroupId) {
	const sceneMap = document.getElementById('scenesvg');
	if (sceneMap == null) return false;
	
	const highRgb = [100, 200, 240];
	
	for (const id of Object.keys(sceneElementColors)) {
	    const element = document.getElementById(id);
	    const parentId = element.parentElement.id;

	    if (parentId == highlightGroupId) {
		// set the color to a highlight color
		const fill = sceneElementColors[id];
		const rgb = fill.substring(4, fill.length-1).split(',');

		const newRgb = [];
		for (let i = 0; i < highRgb.length; i++) {
		    newRgb[i] = Math.floor(0.4*parseInt(rgb[i]) + 0.6*highRgb[i]);
		}
		element.style.fill = `rgb(${newRgb[0]},${newRgb[1]},${newRgb[2]})`
	    } else {
		// Set the color to its original value
		element.style.fill = sceneElementColors[id];
	    }
	}
	return true;
    }
    
    
    //////////////////////////////////////////////////////////////////////
    // The main part of the program
    //////////////////////////////////////////////////////////////////////
    function main() {
	// These are used in the shaders to avoid numerical problems
	const atomMinDist = 1e-4; // Assuming angstroms and high precision
	// units 'sqrt(k*370*K/(hydrogen*u))' 'Å/fs' = 0.01747031
	const atomMaxDisplace = 2.0; // Assuming angstroms and high precision

	// Simulation state
	let currState = SIM_PARAMETERS.initialState;
	let lastState = SimState.none;
	let thermostatOn = true;
	let hideMaterial = -1; // Don't hide anything
	
	/////////////////////////////////////////////////////////////
	// Prepare the WebGL shaders and textures
        /////////////////////////////////////////////////////////////
	const gl = prepareWebGL();

	// External variables set by initSystem()
	let atomTexData = null;
	let simSys = null;
	let thumbnails = null;
	let insertThumbnails = null;
	let shaders = null;
	let textures = null;
	let currTexInfo = null;
	let nextTexInfo = null;
	let renderInfo = null;
	let env = null;
	let step = 0;
	let dynamicsStep = 0;
	let frame = 0;

	// Initialize the mouse/pointer/touch data
	let mouseInfo = new MouseInfo();
	
	function initSystem(baseMoleculeName) {
	    // Prepare the texture data
	    [atomTexData, simSys, thumbnails, insertThumbnails] = prepareAtoms(MOLECULE, SIM_PARAMETERS, baseMoleculeName);

	    // Initialize the shader code
	    shaders = new MDShaders(gl, atomTexData, atomMinDist, atomMaxDisplace);

	    // Initialize the textures and frame buffers
	    textures = new MDTextures(gl, atomTexData, shaders, simSys);

	    // Structure for cycling textures
	    currTexInfo = {posFB: textures.posFB0, pos: textures.pos0, randomFB: textures.randomFB0, random: textures.random0};
	    nextTexInfo = {posFB: textures.posFB1, pos: textures.pos1, randomFB: textures.randomFB1, random: textures.random1};

	    // Initialize the view
	    renderInfo = new RenderInfo(gl, simSys, SIM_PARAMETERS.fieldOfView, SIM_PARAMETERS.cameraPosFactor, SIM_PARAMETERS.lightPosFactor, SIM_PARAMETERS.fadeFactor);

	    // Pack these objects together to avoid function calls with too many parameters.
	    // - Maybe the more idiomatic way in JavaScript is to define the functions here
	    // so these variables would be in scope.
	    // - But this makes me uncomfortable because it isn't clear
	    // what the inputs to the function are.
	    // - However, we have to update this when we alter the textures
	    env = {gl: gl, shaders: shaders, textures: textures, simSys: simSys};

	    // Simulation step counters
	    step = 0;
	    dynamicsStep = 0;
	    frame = 0;
	}
	// Create the initial system based on the baseMolecule
	initSystem(SIM_PARAMETERS.baseMoleculeName);

	
	/////////////////////////////////////////////////////////////
	// User interface functions
        /////////////////////////////////////////////////////////////

	// Prepare the svg scene for interaction, IF IT EXISTS
	// For example, the beach scene allows you to click the water
	// and see a molecular representation of seawater, or click
	// the palm tree and see a molecular representation of wood
	const sceneMap = document.getElementById('scenesvg');
	const sceneElementColors = {};
	// The scene map may or may not exist
	if (sceneMap != null) {
	    const sceneRegions = sceneMap.childNodes;
	    // Check for clicks on the scene
	    for (const region of sceneRegions) {
		if (region.nodeName == 'g') {
		    region.addEventListener('mouseover', function() {enterSceneRegion(region)});
		    region.addEventListener('mouseleave', function() {exitSceneRegion(region)});
		    region.addEventListener('mousedown', function() {clickSceneRegion(region)});

		    // Get the original colors of all child objects
		    // We'll reset them after highlighting
		    for (const child of region.childNodes) {
			if (child.nodeName == 'path' || child.nodeName == 'rect') {
			    //console.log(child.id, child.nodeName, child.style.fill);
			    sceneElementColors[child.id] = child.style.fill;
			}
		    }
		}
	    }
	    setSceneColors(sceneElementColors);
	}

	// Handle movements of the mouse on the WebGL canvas or touchmove events
	function mouseMove(event) {
	    // Stop window from scrolling
	    event.preventDefault();

	    // Is this a touch or a mouse move?
	    // Get the mouse position
	    let moveX = null;
	    let moveY = null;
	    let touchEvent = false;
	    const rect = gl.canvas.getBoundingClientRect();
	    const bodyRect = document.body.getBoundingClientRect()
	    if(event.type.includes('touch')) {
		const { touches, changedTouches } = event.originalEvent ?? event;
		const touch = touches[0] ?? changedTouches[0];
		touchEvent = true;
		mouseInfo.x = touch.pageX - rect.left - window.pageXOffset;
		mouseInfo.y = touch.pageY - rect.top - window.pageYOffset;
		//mouseStatus.innerHTML = 'touchmove';
	    } else if (event.type.includes('mouse')) {
		mouseInfo.x = event.clientX - rect.left;
		mouseInfo.y = event.clientY - rect.top;
	    }
	    mouseInfo.toWorld(gl, renderInfo.viewProjectionMatrix);


	    // Dragging displacement
	    const dx = mouseInfo.x - mouseInfo.dragOriginX;
	    const dy = mouseInfo.y - mouseInfo.dragOriginY;

	    // Require a minimum mouse displacement to start dragging
	    if (mouseInfo.down && (Math.abs(dx) > mouseInfo.minMove || Math.abs(dy) > mouseInfo.minMove)) {
		// Don't drag if no atom is selected or it is the background material
		if (validSelection()) {
		    // Set the state to dragging
		    mouseInfo.dragging = true;
		    if (currState == SimState.running ) {
			// Minimization is a bit slower
			runSteps = Math.ceil(runSteps/2);
		    }
     		    currState = SimState.minimizing;
     		    updatePlayButton();

		    // Move the molecular fragment
		    // This is where we move the atoms to match the mouse
		    const desiredPos = [mouseInfo.scale*mouseInfo.worldX, mouseInfo.scale*mouseInfo.worldY,
					mouseInfo.worldZ];
		    const selectInfo = {id: mouseInfo.fragment, mask: [1, 0, 0, 0], atom: mouseInfo.atom};
		    runDisplaceRelativeShader(env, currTexInfo.pos, selectInfo, desiredPos, nextTexInfo.posFB);
		    //runDisplaceShader(env, currTexInfo.pos, selectInfo, [0.1, 0, 0], nextTexInfo.posFB); // test
		    swapPositionBuffers();
		}
	    }

	    // If this is a new touch event, pause and set the drag origin
	    // It is important that this block is before the touch drag block below
	    if (touchEvent && !mouseInfo.down && mouseInfo.x >= 0 && mouseInfo.y >= 0 && mouseInfo.x < rect.width && mouseInfo.y < rect.height) {
		// Get the atom under the mouse
		mousePoint();

		// Pause the simulation
		lastState = currState;
		currState = SimState.none;
		updatePlayButton();

		mouseInfo.down = true;
		mouseInfo.dragOriginX = mouseInfo.x;
		mouseInfo.dragOriginY = mouseInfo.y;
		mouseInfo.dragWorldX = mouseInfo.worldX;
		mouseInfo.dragWorldY = mouseInfo.worldY;
		mouseInfo.dragWorldZ = mouseInfo.atomZ;
	    }

	    // Is the mouse outside the element? This can happen with touch screens
	    if (mouseInfo.x < 0 || mouseInfo.y < 0 || mouseInfo.x >= rect.width || mouseInfo.y >= rect.height) {
		if (mouseInfo.dragging) {
		    // Delete the molecule that was dragged out of the canvas
		    deleteSelection();
		}
		mouseEnd();
	    }
	}

	// The mouse or touch pointer has left the WebGL canvas or the mouse button has been released
	function mouseEnd() {
	    // If we were dragging, run some minimization
	    if (mouseInfo.dragging) {
		minimize(200);
		mouseInfo.dragging = false;
	    }
	    
	    // Continue whatever we were doing
	    if (mouseInfo.down) {
		mouseInfo.down = false;
		currState = lastState;
		updatePlayButton();
	    }
	    mouseInfo = new MouseInfo(mouseInfo.scale);
	    updateSelectionPanel();
	}

	// Here we determine what the mouse is pointing at in the WebGL canvas
	// This requires running a WebGL shader to render a single pixel under the mouse
	function mousePoint() {
	    // Get the mouse position
	    const mouseMatrix = getMouseMatrix(gl, mouseInfo.x, mouseInfo.y, renderInfo);
	    const renderInfoPick = new RenderInfo(gl, simSys, SIM_PARAMETERS.fieldOfView, SIM_PARAMETERS.cameraPosFactor);
	    renderInfoPick.viewProjectionMatrix = mouseMatrix;

	    // What atom is the mouse on?
	    // Render a single pixel to the pick buffer, coloring by (x, y, z, atomIndex+1)
	    // shaders.drawPositionsFS has flat coloring with float colors of (x, y, z, atomIndex+1)
	    const selectInfo = {id: mouseInfo.fragment, mask: [1, 0, 0, 0], atom: mouseInfo.atom, scale: 1.0, hideMaterial: hideMaterial};
	    drawAtoms(env, shaders.drawPositions, shaders.drawPositionsLocs, currTexInfo.pos, selectInfo, renderInfoPick, textures.pickFB, 1, 1);
	    const pickData = extractData(env, textures.pick, 1, 1);
	    // drawPositions program puts the atom index in w component
	    mouseInfo.atom = Math.floor(pickData[3] + 0.5) - 1; // Convert back to 0-based index
	    mouseInfo.atomX = pickData[0];
	    mouseInfo.atomY = pickData[1];
	    mouseInfo.atomZ = pickData[2];
	    // Calculate a fudge factor to convert the mouse position (in canvas pixels)
	    // to the simulation world coordinates.
	    // This is an empirical kludge since I haven't figured out how to revert the projection correctly
	    if (Math.abs(mouseInfo.atomX) > 5.0) {
		mouseInfo.scale = mouseInfo.atomX/mouseInfo.worldX;
	    } else if (Math.abs(mouseInfo.atomY) > 5.0) {
		mouseInfo.scale = mouseInfo.atomY/mouseInfo.worldY;
	    }
	    
	    // Set the fragment, monomer, type, resName of the selected atom
	    if (mouseInfo.atom >= 0 && mouseInfo.atom < atomTexData.typeList.length) {
		// The fragment, monomer, type, and resName are taken from JavaScript arrays
		// using the atom index
		mouseInfo.fragment = atomTexData.select[4*mouseInfo.atom];
		mouseInfo.monomer = atomTexData.select[4*mouseInfo.atom+1];
		mouseInfo.material = atomTexData.select[4*mouseInfo.atom+2];
		mouseInfo.type = atomTexData.typeList[mouseInfo.atom];
		mouseInfo.resName = atomTexData.resNameList[mouseInfo.atom];
	    } else {
		// Set values that represent that nothing is selected
		mouseInfo.fragment = -1;
		mouseInfo.monomer = -1;
		mouseInfo.material = -1;
		mouseInfo.type = 'UNK';
		mouseInfo.resName = 'none';
	    }
	    
	    updateSelectionPanel();
	}

	function validSelection() {
	    // Don't select if there is no atom
	    if (mouseInfo.atom < 0 || mouseInfo.fragment < 0 || mouseInfo.material < 0 || mouseInfo.material == MaterialEnum.background) {
		return false;
	    } else {
		return true;
	    }
	}

	/////////////////////////////////////////////////////////////
	// Handle interaction with a scene
	// For example, the beach scene allows you to click the water
	// and see a molecular representation of seawater, or click
	// the palm tree and see a molecular representation of wood
	function enterSceneRegion(sceneGroup) {
	    setSceneColors(sceneElementColors, sceneGroup.id);
	}
	function exitSceneRegion(sceneGroup) {
	    setSceneColors(sceneElementColors, null);
	}
	function clickSceneRegion(sceneGroup) {
	    const map = SIM_PARAMETERS.sceneMoleculeMap;
	    if (Object.hasOwn(map, sceneGroup.id)) {
		const baseMoleculeName = map[sceneGroup.id];
		console.log(`Clicked ${sceneGroup.id}, setting system to ${baseMoleculeName}`);
		if (Object.hasOwn(MOLECULE, baseMoleculeName)) {
		    // Reinitialize everything!
		    initSystem(baseMoleculeName);
		    // This is an inelegant hardcoded thermostat setting for air
		    thermostatOn = (baseMoleculeName != 'air');
		    updateThermostatCheck();
		} else {
		    console.log(`ERROR! baseMolecule ${baseMoleculeName} has not been defined. Make sure the name is correct and that the script containing this molecule has been loaded in the HTML file.`)
		}
	    } else {
		console.log(`ERROR! No molecule in sceneMoleculeMap for scene id ${sceneGroup.id}`);
	    }

	    // Set the system name
	    document.getElementById('systemName').textContent = sceneGroup.id;
	}

	/////////////////////////////////////////////////////////////
	// Mouse events
	gl.canvas.addEventListener('mousemove', mouseMove);
	// On click, pause the simulation and store the clicked mouse position
	// Don't start moving atoms until the mouse moves while down
	gl.canvas.addEventListener('mousedown', (e) => {
	    const rect =gl.canvas.getBoundingClientRect();
	    mouseInfo.x = e.clientX - rect.left;
	    mouseInfo.y = e.clientY - rect.top;

	    // Pause the simulation
	    lastState = currState;
	    currState = SimState.none;
	    updatePlayButton();
	    
	    mouseInfo.down = true;
	    mouseInfo.dragOriginX = mouseInfo.x;
	    mouseInfo.dragOriginY = mouseInfo.y;
	    mouseInfo.dragWorldX = mouseInfo.worldX;
	    mouseInfo.dragWorldY = mouseInfo.worldY;
	    mouseInfo.dragWorldZ = mouseInfo.atomZ;
	});
	// We've released the mouse, so continue what we were doing
	gl.canvas.addEventListener('mouseup', mouseEnd);
	
	// Mouse has left the viewport
	gl.canvas.addEventListener('mouseout', (e) => {
	    if (mouseInfo.dragging) {
		// Delete the molecule that was dragged out of the canvas
		deleteSelection();
	    }
	    mouseEnd();
	});

	
	/////////////////////////////////////////////////////////////
	// Touch events
	gl.canvas.addEventListener('touchmove', mouseMove, {passive: false});
	gl.canvas.addEventListener('touchstart', mouseMove);
	gl.canvas.addEventListener('touchcancel', mouseEnd);
	gl.canvas.addEventListener('touchend', mouseEnd);

	/////////////////////////////////////////////////////////////
	// Handle keyboard input
	window.addEventListener('keydown', function(event) {
	    const dispMag = 2.0; // in Å
	     switch (event.key) {
             case 'w':
		 // up
		 displaceSel([0, dispMag, 0]);
                 break;

             case 'a':
		 // left
		 displaceSel([-dispMag, 0, 0]);
                 break;
		 
             case 's':
		 // down
		 displaceSel([0, -dispMag, 0]);
                 break;
		 
             case 'd':
		 // right
		 displaceSel([dispMag, 0, 0]);
                 break;
            
             case 'f':
		 // toward
		 displaceSel([0, 0, dispMag]);
		 break;

	     case 'b':
		 // away
		 displaceSel([0, 0, -dispMag]);
		 break;

	     case 'Delete':
	     case 'Backspace':
		 deleteSelection();
		 break;
             }
	});
	
	// Keyboard displacements
	function displaceSel(displace) {
	    // Don't allow moving the background
	    if (mouseInfo.material == MaterialEnum.background) return;
	    
	    // Do the displacement in multiple stages
	    const dispMag = Math.sqrt(displace[0]**2 + displace[1]**2 + displace[2]**2);
	    const stages = Math.ceil(dispMag/0.5);
	    const stageDisplace = [displace[0]/stages, displace[1]/stages, displace[2]/stages];

	    // Displace and minimize the energy
	    const selectInfo = {id: mouseInfo.fragment, mask: [1, 0, 0, 0], atom: mouseInfo.atom};
	    for (let s = 0; s < stages; s++) {
		runDisplaceShader(env, currTexInfo.pos, selectInfo, stageDisplace, nextTexInfo.posFB);
		swapPositionBuffers();
		minimize(100); // a little quicker
	    }
	}
	

	/////////////////////////////////////////////////////////////
	// Buttons, checkboxes, and input fields
	
	// Handle play button
	window.togglePlay = function() {
	    if (currState == SimState.running) {
		currState = SimState.none;
	    } else {
		currState = SimState.running;
	    }
	    updatePlayButton();
	}
	function updatePlayButton() {
	    if (currState == SimState.running) {
		// Pause symbol
		document.getElementById("buttonPlay").innerHTML = "||";
		// Looks bad on mobile
		// document.getElementById("buttonPlay").innerHTML = "&#9612;&nbsp;&#9612;";
	    } else {
		// Play symbol
		document.getElementById("buttonPlay").innerHTML = "▶";
	    }
	    
	    if (currState == SimState.minimizing) {
		document.getElementById("buttonMinimize").innerHTML = "Minimizing...";
	    } else {
		document.getElementById("buttonMinimize").innerHTML = "Minimize energy";
	    }
	}
	updatePlayButton();

	// Handle thermostat
	window.toggleThermostat = function() {
	    thermostatOn = !thermostatOn;
	    updateThermostatCheck();
	}
	function updateThermostatCheck() {
	    const checkbox =  document.getElementById('checkThermostat');
	    checkbox.checked = thermostatOn;
	}
	updateThermostatCheck();

	function isReal(str) {
	    if (typeof str != 'string') return false;
	    return !isNaN(str) && !isNaN(parseFloat(str));
	}
	function updateThermostatTemperature() {
	    const inputTemperStr = document.getElementById('inputThermostat').value;
	    // Check that it's a number
	    if (isReal(inputTemperStr)) {
		const temper = parseFloat(inputTemperStr);
		simSys.setLangevin(temper, SIM_PARAMETERS.langevinDamping);
	    }
	}
	// Set the initial temperature in the input form
	document.getElementById('inputThermostat').value = simSys.temper;

	// Run a shader to reset the velocities
	window.pressResetVelocities = function() {
	    advanceRandomNumbers();
	    runResetVelocityShader(env, currTexInfo.random, simSys.kT, textures.velFB);
	}

	// Handle the insertion panel
	let insertIndex = 0;
	function updateInsertPanel() {
	    if (insertIndex < insertThumbnails.length) {
		document.getElementById('insertName').textContent = insertThumbnails[insertIndex].name;
		insertThumbnails[insertIndex].draw(insertContext2d);
	    }
	}
	window.pressPreviousMolecule = function() {
	    insertIndex--;
	    if (insertIndex < 0) insertIndex = insertThumbnails.length-1;
	    updateInsertPanel();
	}
	window.pressNextMolecule = function() {
	    insertIndex++;
	    if (insertIndex >= insertThumbnails.length) insertIndex = 0;
	    updateInsertPanel();
	}
	window.pressInsertMolecule = function() {
	    if (insertIndex < insertThumbnails.length) {
		const thumb = insertThumbnails[insertIndex];
		insertMolecule(thumb.mol, thumb.name);
	    }
	}

	// Write the PDB text
	window.pressPDB = function() {
	    // Get the newest positions
	    const currPos = extractData(env, currTexInfo.pos, textures.width, textures.height);
	    atomTexData.setPositions(currPos);
	    const pdbData = document.getElementById('pdbData');
	    pdbData.innerHTML = "PDB_DATA:\n" + atomTexData.writePDB();
	}

	// Display water or not
	window.toggleShowWater = function() {
	    hideMaterial = (hideMaterial == MaterialEnum.water) ? -1 : MaterialEnum.water;
	    updateShowWaterCheck();
	}
	function updateShowWaterCheck() {
	    const checkbox =  document.getElementById('checkShowWater');
	    checkbox.checked = (hideMaterial != MaterialEnum.water);
	}
	updateShowWaterCheck();

	// Minimization button
	window.pressMinimize = function() {
	    if (currState == SimState.running ) {
		// Minimization is a bit slower
		runSteps = Math.ceil(runSteps/2);
	    }
	    
	    if ( currState != SimState.minimizing ) {
		currState = SimState.minimizing;
	    } else {
		currState = SimState.none;
	    }
	    updatePlayButton()
	}

	// Show a 2D diagram of the selected molecule
	function updateSelectionPanel() {
	    // Set the name of the selected residue
	    if (Object.hasOwn(SIM_PARAMETERS.resNameMap, mouseInfo.resName)) {
		mouseInfo.moleculeName = SIM_PARAMETERS.resNameMap[mouseInfo.resName];
		document.getElementById('moleculeName').innerHTML = mouseInfo.moleculeName;
	    } else {
		mouseInfo.moleculeName = mouseInfo.resName;
		document.getElementById('moleculeName').innerHTML = `<i>${mouseInfo.resName}</i>`;
	    }
	    if (Object.hasOwn(SIM_PARAMETERS.altNameMap, mouseInfo.resName)) {
		document.getElementById('alternateName').innerHTML = `other names: ${SIM_PARAMETERS.altNameMap[mouseInfo.resName]}`;
	    } else {
		document.getElementById('alternateName').innerHTML = '&nbsp;';
	    }
	    
	    // Draw the selected residue in the 2d canvas
	    if (Object.hasOwn(thumbnails, mouseInfo.resName)) {
		thumbnails[mouseInfo.resName].draw(selectContext2d);
	    } else {
		selectContext2d.clearRect(0, 0, selectContext2d.canvas.width, selectContext2d.canvas.height);
	    }
	}
	
	/////////////////////////////////////////////////////////////
	// Functions that modify the system used by the user interface
        /////////////////////////////////////////////////////////////

	// Run a shader to minimize energy
	// This gets called in many different places
	function minimize(steps = 200) {
	    let stepSize = 0.05;
	    const minStepSize = 1e-4;
	    const stepScale = Math.pow(minStepSize, 1.0/steps);
	    for (let i = 0; i < steps; i++) {
		// Calculate the forces
		// false: We use the slower version for minimization since there is no near cutoff
		runForceShader(env, currTexInfo.pos, textures.forceFB, false);
		// Take a step downhill
		runDownhillShader(env, currTexInfo.pos, textures.force, stepSize, nextTexInfo.posFB);
		// Swap to the next positions
		swapPositionBuffers();
		// Take smaller steps each time
		stepSize *= stepScale;
	    }
	    advanceRandomNumbers();
	    runResetVelocityShader(env, currTexInfo.random, simSys.kT, textures.velFB);
	}
	
	// For Verlet, we keep the random number and position buffers together
	// But there are times when we want to swap one and not the other
	function advanceRandomNumbers() {
	    // Cycle the random numbers
	    runUpdateRandomShader(env, currTexInfo.random, nextTexInfo.randomFB);
	    const tempTex = currTexInfo.random;
	    currTexInfo.random = nextTexInfo.random;
	    nextTexInfo.random = tempTex;

	    const tempFB = currTexInfo.randomFB;
	    currTexInfo.randomFB = nextTexInfo.randomFB;
	    nextTexInfo.randomFB = tempFB;
	}

	// Cycle positions
	function swapPositionBuffers() {
	    const tempTex = currTexInfo.pos;
	    currTexInfo.pos = nextTexInfo.pos;
	    nextTexInfo.pos = tempTex;
	    
	    const tempFB = currTexInfo.posFB;
	    currTexInfo.posFB = nextTexInfo.posFB;
	    nextTexInfo.posFB = tempFB;
	}

	
	// Insert a molecule into the system
	function insertMolecule(mol, name) {	    
	    // Get the newest positions
	    const currPos = extractData(env, currTexInfo.pos, textures.width, textures.height);

	    // Get the initial number of collisions in the hash table
	    const oldMaxCollisions = atomTexData.maxCollisions;

	    // Update the texture data
	    let insertPos = null;
	    if (Object.hasOwn(SIM_PARAMETERS, 'insertPos')) {
		insertPos = SIM_PARAMETERS.insertPos;
	    } else {
		insertPos = [-simSys.wall[0]*0.8+simSys.radius, -simSys.wall[1]*0.8+simSys.radius, 1.0];
	    }
	    const ok = atomTexData.insertMolecule(mol, currPos, insertPos, true);

	    if (ok) {
		// Has the maximum number of collisions changed?
		// If so we need to recompile the shaders using the hash table
		if (atomTexData.maxCollisions > oldMaxCollisions) {
		    console.log(`Maximum collisions in the hash table changed from ${oldMaxCollisions} to ${atomTexData.maxCollisions} on insertion of ${name}. Reinitializing hash table and shaders...`);
		    // We might as well regenerate the hash table while we're at it
		    // Maybe we can improve maxCollisions
		    atomTexData.makeHashTable();
		    
		    // Reinitialize the shader code with potentially different hash table dimensions and maxCollisions
		    shaders = new MDShaders(gl, atomTexData, atomMinDist, atomMaxDisplace);
		    // We have to remember to update this reference in the "env" variable
		    env.shaders = shaders;
		}
		
		// Reinitialize the textures and frame buffers
		textures = new MDTextures(gl, atomTexData, shaders, simSys);
		
		// We have to remember to update this reference in the "env" variable
		env.textures = textures;

		// We also have to update the references in the texture cycler
		currTexInfo = {posFB: textures.posFB0, pos: textures.pos0, randomFB: textures.randomFB0, random: textures.random0};
		nextTexInfo = {posFB: textures.posFB1, pos: textures.pos1, randomFB: textures.randomFB1, random: textures.random1};
		
		// Relax everything (also resets velocities)
		minimize(500);

		selectStatus.innerHTML = `<span class="selectInsert">Inserting <b>${name}</b><br>&nbsp;</span>`;
		selectStatusTime = performance.now();
	    } else {
		selectStatus.innerHTML = `<span class="selectDelete"><i>Could not insert! Too many atoms.</i></span>`;
		selectStatusTime = performance.now();
	    }
	}

	// Delete a molecule from the system
	function deleteSelection() {
	    // Check that there is a valid selected fragment
	    if (validSelection()) {		
		// Show the deletion in the selection status
		const selectStatus = document.getElementById('selectStatus');
		selectStatus.innerHTML = `<span class="selectDelete">Deleting <b>${mouseInfo.moleculeName}</b><br>&nbsp;</span>`;
		selectStatusTime = performance.now();
		
		// Get the newest positions
		const currPos = extractData(env, currTexInfo.pos, textures.width, textures.height);
		//console.log(`Deleting ${mouseInfo.moleculeName} (fragment ${mouseInfo.fragment}). Position: ${currPos[0]} ${currPos[1]} ${currPos[2]}`);

		// Update the texture data
		const selectTexColumn = 0; // select RGBA format: [fragment, monomer, material, radius]
		atomTexData.deleteMolecule(currPos, mouseInfo.fragment, selectTexColumn);
		
		// Reinitialize the textures and frame buffers
		textures = new MDTextures(gl, atomTexData, shaders, simSys);

		// We have to remember to update this reference in the "env" variable
		env.textures = textures;

		// We also have to update the references in the texture cycler
		currTexInfo = {posFB: textures.posFB0, pos: textures.pos0, randomFB: textures.randomFB0, random: textures.random0};
		nextTexInfo = {posFB: textures.posFB1, pos: textures.pos1, randomFB: textures.randomFB1, random: textures.random1};
	
		// Relax everything
		minimize();
	    }
	}

	
	
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	// Prepare for the simulation loop
	let runSteps = 20;
	const analysisFrames = 20;
	let selectStatusTime = performance.now();
	let sumSimTime = 0.0;

	// 2D canvas for molecule identification
	const selectContext2d = document.getElementById('selectCanvas').getContext('2d');
	const insertContext2d = document.getElementById('insertCanvas').getContext('2d');
	updateInsertPanel();

	///////////////////////////////////////////////////////////////////////////////
	// Simulation loop
	function run() {
	    // Run multiple steps before rendering
	    const frameStartTime = performance.now();
	    if (currState == SimState.minimizing) {
		// Run minimization
		minimize(runSteps);
		step+=runSteps;
	    } else if (currState == SimState.running) {
		// Velocity Verlet for Langevin dynamics from
		// Izaguirre et al. (2001) Langevin stabilization of molecular dynamics. J. Chem. Phys. 114(5), 2090–2098. DOI: 10.1063/1.1332996
		for (let i = 0; i < runSteps; i++) {
		    if (thermostatOn) {
			// First half-kick
			// We correctly use same random numbers as the previous second half-kick
			runUpdateVelocityShader(env, shaders.firstKick, shaders.firstKickLocs, textures.vel,
						textures.force, currTexInfo.random, textures.halfVelFB);
		    } else {
			runUpdateVelocityShader(env, shaders.verletKick, shaders.verletKickLocs, textures.vel,
						textures.force, currTexInfo.random, textures.halfVelFB);
		    }
		    
		    // Drift
		    runUpdatePositionShader(env, currTexInfo.pos, textures.halfVel, nextTexInfo.posFB);
		    
    		    // Calculate the new forces
    		    // F(t + dt) = f(r(t+dt))
		    // true: fast version for dynamics
   		    runForceShader(env, nextTexInfo.pos, textures.forceFB, true);

		    // Second half-kick
		    if (thermostatOn) {
			// The second half-kick uses the t+dt random numbers
    			runUpdateRandomShader(env, currTexInfo.random, nextTexInfo.randomFB);
			runUpdateVelocityShader(env, shaders.secondKick, shaders.secondKickLocs, textures.halfVel,
						textures.force, nextTexInfo.random, textures.velFB);
		    } else {
			runUpdateVelocityShader(env, shaders.verletKick, shaders.verletKickLocs, textures.halfVel,
						textures.force, nextTexInfo.random, textures.velFB);
		    }
		    
		    // Cycle between textures to read to and textures to write to
		    {
      			const tempInfo = currTexInfo;
			currTexInfo = nextTexInfo;
      			nextTexInfo = tempInfo;
		    }
		    step++;
		    dynamicsStep++;
		} // Run multiple step loop
	    }
	    // Sum the amount of time just running the simulation loop
	    sumSimTime += performance.now() - frameStartTime;

	    // At the beginning of the simulation, reset the velocities after a relaxation 
	    if (frame == analysisFrames-1) {
		advanceRandomNumbers();
		runResetVelocityShader(env, currTexInfo.random, simSys.kT, textures.velFB);
	    } else if (frame % analysisFrames == analysisFrames-1) {
		// Calculate performance and energies
		// Set the temperature from the input form
		updateThermostatTemperature();

		// Calculate total performance
		const t_ps = simSys.timestep*dynamicsStep*0.001;
		const t = performance.now();
		const totalTime = t - lastTime;
		const totalSteps = step - lastStep;
		const stepTime = sumSimTime/totalSteps;
		const runPerf = 86.4*simSys.timestep/stepTime; // ns/day
		const totalPerf = 86.4*simSys.timestep/(totalTime/totalSteps); // ns/day
		const totalPerfPs = 86.4*simSys.timestep/(totalTime/totalSteps)/1.44 // ps/min
		const framesPerSecond = 1000.0*analysisFrames/totalTime;
		lastTime = t;
		lastStep = step;
		sumSimTime = 0.0;

		// Set the performance text
		let statusText = "";
		if (totalSteps > 0) {
		    statusText += `Time per simulation step: ${stepTime.toFixed(5)} ms`;
		    statusText += `<br>Total performance: ${totalPerfPs.toFixed(2)} ps/min`;
		    statusText += `<br>Total performance: ${totalPerf.toFixed(2)} ns/day`;
		    statusText += `<br>Simulation performance: ${runPerf.toFixed(2)} ns/day`;
		}
		statusText += `<br>Framerate: ${framesPerSecond.toFixed(2)} 1/s`;
		statusText += `<br>Steps per frame: ${runSteps}`;
		const status = document.getElementById("status");
		status.innerHTML = statusText;

		// Get the kinetic and potential energies from the textures
		// This comes after the swap, so we use currTexInfo
		const currPos = extractData(env, currTexInfo.pos, textures.width, textures.height);
		const currVel = extractData(env, textures.vel, textures.width, textures.height);
		const currForce = extractData(env, textures.force, textures.width, textures.height);
		const currNonbond = extractData(env, textures.nonbond, textures.width, textures.height);
		const [enerKin, temperKin] = calcKineticEnergy(currPos, currVel, currNonbond, simSys);
		const enerPot = calcPotentialEnergy(currPos, currForce);
		const enerTotal = enerKin + enerPot;
		
		// Write the status table
		let prec = 2;
		//document.getElementById("tabState").innerHTML = `${SimStateName[currState]}`;
		document.getElementById("tabTime").innerHTML = `${t_ps.toFixed(prec)}`;
		document.getElementById("tabEnergy").innerHTML = `${enerTotal.toFixed(prec)}`;
		document.getElementById("tabTemperature").innerHTML = `${temperKin.toFixed(prec)}`;
		
		// Adjust runSteps to improve performance
		if (currState != SimState.none) {
		    if (framesPerSecond < 10.0) {
			runSteps = Math.ceil(runSteps*framesPerSecond/10.0);
		    } else if (framesPerSecond > 25.0) {
			runSteps = Math.ceil(runSteps*framesPerSecond/25.0);
		    }
		    // Put some kind of sane limits on runSteps
		    if (runSteps < 5) {
			runSteps = 5;
		    } else if (runSteps > 500) {
			runSteps = 500;
		    }
		}
	    }

	    // Draw the system
	    // Is the mouse in the 3d canvas?
	    // Update the world position of the mouse (at the selected atom's z value)
	    if (mouseInfo.x > 0 && mouseInfo.y > 0) {
		mouseInfo.toWorld(gl, renderInfo.viewProjectionMatrix);
		if (!mouseInfo.down) {
		    mousePoint();
		}
	    }
		
	    // Set the selection status
	    const selectStatus = document.getElementById('selectStatus');
	    if (mouseInfo.dragging) {
		selectStatusTime = performance.now();
		selectStatus.innerHTML = `<span class="selectDrag">Dragging ${mouseInfo.moleculeName}</span><br>&nbsp;`;
	    } else if (mouseInfo.fragment >= 0 && mouseInfo.material != MaterialEnum.background) {
		selectStatusTime = performance.now();
		// selectStatus.innerHTML = `<span class="selectSelect">Selecting fragment ${mouseInfo.fragment}</span><br>&nbsp;`;
		selectStatus.innerHTML = `<span class="selectSelect">Selecting ${mouseInfo.moleculeName}</span><br>&nbsp;`;
	    } else {
		// No selection action
		// Leave the last update for a little bit
		if (performance.now() - selectStatusTime > 1000.0) {
		    selectStatusTime = performance.now();
		    //selectStatus.innerHTML = `To move molecule, drag with mouse. Press "F" or "B" to move toward/away. Press "Del" or drag off the viewport to delete.`;
		    selectStatus.innerHTML = 'Drag to move. Drag off or "del" to delete. F/B: move away/toward';
		}		
	    }
	    
	    // Set mouse status info
	    //const mouseStatus1 = document.getElementById("mouseStatus1");
	    //mouseStatus1.innerHTML = `<br>Mouse: ${mouseInfo.x} ${mouseInfo.y}<br>Atom: ${mouseInfo.atom}<br>Atom pos.: ${mouseInfo.atomX.toFixed(2)} ${mouseInfo.atomY.toFixed(2)} ${mouseInfo.atomZ.toFixed(2)}<br>World: ${mouseInfo.worldX.toFixed(2)} ${mouseInfo.worldY.toFixed(2)} ${mouseInfo.worldZ.toFixed(2)}<br>Fragment: ${mouseInfo.fragment}<br>Residue: ${mouseInfo.resName}<br>Monomer: ${mouseInfo.monomer}<br>Material: ${mouseInfo.material}<br>Type: ${mouseInfo.type}<br>Mouse scale: ${mouseInfo.scale.toFixed(2)}`;
	    
	    // Make the selection pulse with time (using frame makes it stutter)
	    const selectScale = 0.85 + 0.6*Math.cos(2.0*Math.PI*(performance.now()-cycleStartTime)/800.0);
	    // Highlight the whole fragment
	    //const selectInfo = {id: mouseInfo.fragment, mask: [1, 0, 0, 0], atom: mouseInfo.atom, scale: selectScale, hideMaterial: hideMaterial};
	    // Highlight just the monomer
	    const selectInfo = {id: mouseInfo.monomer, mask: [0, 1, 0, 0], atom: mouseInfo.atom, scale: selectScale, hideMaterial: hideMaterial};
	    // Draw the scene
	    // We have already cycled from nextTexInfo to currTexInfo, so use currTexInfo	    
	    drawAtoms(env, shaders.drawSpheres, shaders.drawSpheresLocs, currTexInfo.pos, selectInfo, renderInfo, null);
	    //drawAtoms(env, shaders.drawPositions, shaders.drawPositionsLocs, currTexInfo.pos, selectInfo, renderInfo, null);
	    frame++;

	    // Loop continuously 
	    if (!stopSimulation) requestAnimationFrame(run);
	}
	
	///////////////////////////////////////////////////////////////////////////////
	// Begin time loop
	const cycleStartTime = performance.now();
	let lastTime = cycleStartTime;
	let lastStep = 0;
	let stopSimulation = false;
	requestAnimationFrame(run);
    }
    
    main();
})();
