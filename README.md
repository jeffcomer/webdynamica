WebDynamica
===========

WebDynamica is a browser-based interactive molecular dynamics program using WebGL.

The goal of WebDynamica is to help educate the public and students about molecular interactions, while being useful to experts in molecular simulation for prototyping and constructing initial conditions.

It allows the user to manipulate molecules using the mouse or touchscreen, and runs in any web browser including Firefox, Chrome, and Safari, taking advantage of WebGL and graphics hardware to rapidly perform molecular calculations on devices from mobile phones to high-end desktops. 

Molecules can be moved along the *x*- and *y*-axes by dragging, along the *z*-axis by pressing the "F" and "B" keys. The selected molecule can be deleted by dragging it out of the view window or by pressing "Delete". A wide variety of other molecules can be added to the system using the Insertion tool. 

The current demo begins with benzene and benzoate and benzamidinium ions on top of a graphene sheet. These molecules are represented with the [CHARMM General Force Field](https://doi.org/10.1002/jcc.21367). WebDynamica uses a JSON-like structure to represent molecular structures and associated interatomic interactions. The node.js script `charmm_to_javascript.js` is can produce a JavaScript file in the necessary format from the set of CHARMM input files as would be used with the [NAMD](https://www.ks.uiuc.edu/Research/namd/) molecular dynamics engine. The needed files are:

- CHARMM-format .psf file
- .pdb file (Protein Data Bank format) file containing the coordinates
- .xsc file (NAMD extended system text file: step size_x 0 0 0 size_y 0 0 0 size_z). Non-orthogonal boxes are not supported by WebDynamica
- A list the necessary CHARMM-format parameter files

The usage of `charmm_to_javascript.js` is the following:

`node charmm_to_javascript.js psfFile pdbFile xscFile prmFile0 [prmFile1]... moleculeName [fragments|monomers] outFile`

- The presence of "fragments" means that fragment numbers will be assigned based on the bond network. Note that searching the bond network can be slow.
- The presence of "monomers" means that fragment numbers will be assigned based the segment and resid in the psf.
- With neither fragments nor monomers, no fragment data will be written
- Fragment and monomer data is used by WebDynamica for selection, moving, and deletion
    
For example, to generate the initial structure of the current demo, the following command can be used:

`node charmm_to_javascript.js examples/graph30_neg-pos.psf examples/graph30_neg-pos_pos1.pdb examples/graph30.xsc examples/par_all36_cgenff.prm examples/toppar_water_ions.str graph30_neg-pos fragments examples/graph30_neg-pos.webdyn.js`

WebDynamica makes use of utility functions in `m4.js` and `webgl-utils.js`, which are from [WebGL Fundamentals](https://webglfundamentals.org/). The redistribution conditions for these two files are given in these files. The author of WebDynamica is not affiliated, associated, authorized, endorsed by, or in any way officially connected with the author of `m4.js` and`webgl-utils.js`, GFXFundamentals, or WebGL Fundamentals. The inclusion of these two files and the mention of these names does not imply endorsement by GFXFundamentals or WebGL Fundamentals.
