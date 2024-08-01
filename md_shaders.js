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

////////////////////////////////////////////////////////////////////////////////
// SHADERS
//
// References for Velocity Verlet and Brooks-Brünger-Karplus Langevin dynamics schemes
// Izaguirre et al. (2001) Langevin stabilization of molecular dynamics. J. Chem. Phys. 114(5), 2090–2098. DOI: 10.1063/1.1332996
// Phillips et al. (2005) Scalable Molecular Dynamics with NAMD. J. Comput. Chem. 26(16), 1781–1802. DOI: 10.1002/jcc.20289
// Brooks et al. (1984) Stochastic boundary conditions for molecular dynamics simulations of ST2 water. Chem. Phys. Lett. 105(5), 495–500. DOI: 10.1016/0009-2614(84)80098-6
// I'm purposely not indenting here to make editing the shaders easier
class MDShaders {
constructor(gl, atomTexData, atomMinDist, atomMaxDisplace, forceSplit = 500.0) {
this.SPHERE_VERTICES = 39; // Our simple half-sphere model
this.ATOM_MIN_DIST = atomMinDist.toFixed(1);
// Typically, the maximum displacement should be a few times bigger than
// units 'sqrt(k*370*K/(hydrogen*u))' 'Å/fs' = 0.01747031
this.ATOM_MAX_DISPLACE = atomMaxDisplace.toFixed(1);
this.FORCE_SPLIT = forceSplit.toFixed(1);
// Force and energy should be in different units, but it works well numerically
this.ENERGY_SPLIT = forceSplit.toFixed(1); 

console.log('MDShaders texture size:',atomTexData.width, atomTexData.height);
// Texture dimension strings
this.ATOM_TEX_WIDTH = atomTexData.width.toFixed(1);
this.ATOM_TEX_HEIGHT = atomTexData.height.toFixed(1);
this.BOND_TEX_WIDTH = atomTexData.bondDim.width.toFixed(1);
this.BOND_TEX_HEIGHT = atomTexData.bondDim.height.toFixed(1);
this.BOND_TERMS_PER_ATOM = atomTexData.bondDim.termsPerAtom.toFixed(1);
this.ANGLE_TEX_WIDTH = atomTexData.angleDim.width.toFixed(1);
this.ANGLE_TEX_HEIGHT = atomTexData.angleDim.height.toFixed(1);
this.ANGLE_TERMS_PER_ATOM = atomTexData.angleDim.termsPerAtom.toFixed(1);
this.DIHEDRAL_TEX_WIDTH = atomTexData.dihedralDim.width.toFixed(1);
this.DIHEDRAL_TEX_HEIGHT = atomTexData.dihedralDim.height.toFixed(1);
this.DIHEDRAL_TERMS_PER_ATOM = atomTexData.dihedralDim.termsPerAtom.toFixed(1);
this.EXCLUDE_TEX_WIDTH = atomTexData.excludeDim.width.toFixed(1);
this.EXCLUDE_TEX_HEIGHT = atomTexData.excludeDim.height.toFixed(1);
this.EXCLUDE_TERMS_PER_ATOM = atomTexData.excludeDim.termsPerAtom.toFixed(1);

// Linear congruential generator parameters from
// Pierre L'Ecuyer (1999) Mathematics of Computation, 68(225):249–260
// https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00996-5/
// We choose with small enough moduli (m) and multipliers (a) that 
// exact integer math with 32-bit floats is possible
this.RAND_M0 = 32749;
this.RAND_A0 = 209;
this.RAND_M1 = 16381;
this.RAND_A1 = 665;
this.RAND_M2 = 16381;
this.RAND_A2 = 572;
this.RAND_M3 = 8191;
this.RAND_A3 = 1716;

// Render the atoms using a half-sphere model with 13 triangles
this.drawSpheresVS = `
#define SPHERE_VERTICES ${this.SPHERE_VERTICES}
#define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
#define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

attribute float id;
// Atom position texture
uniform sampler2D posTex;
uniform sampler2D colorTex;
uniform sampler2D selectTex;
uniform float selectId;
uniform vec4 selectMask;
uniform float selectScale;
uniform float hideMaterial;

// View geometry
uniform mat4 matrix;
uniform vec3 cameraPos;
uniform vec3 lightPos;
varying vec3 v_surfaceToView;
varying vec3 v_surfaceToLight;
varying vec3 v_normal;
varying vec3 v_color;
varying float v_background;
varying float v_atomIndex;

vec2 atomTextureCoords(float atomIndex) {
  float texX = (mod(atomIndex, ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_WIDTH;
  float texY = (floor(atomIndex/ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_HEIGHT;
  return vec2(texX, texY);
}


void main() {
  float vertices = float(SPHERE_VERTICES);

  vec3 d[SPHERE_VERTICES];
  d[0]  = vec3(0.000, -0.500, 0.866);
  d[1]  = vec3(0.433, 0.250, 0.866);
  d[2]  = vec3(-0.433, 0.250, 0.866);
  d[3]  = vec3(-0.000, 1.000, 0.000);
  d[4]  = vec3(-0.433, 0.250, 0.866);
  d[5]  = vec3(0.433, 0.250, 0.866);
  d[6]  = vec3(-0.000, 1.000, 0.000);
  d[7]  = vec3(0.433, 0.250, 0.866);
  d[8]  = vec3(0.643, 0.766, 0.000);
  d[9]  = vec3(0.985, 0.174, 0.000);
  d[10] = vec3(0.643, 0.766, 0.000);
  d[11] = vec3(0.433, 0.250, 0.866);
  d[12] = vec3(0.866, -0.500, 0.000);
  d[13] = vec3(0.985, 0.174, 0.000);
  d[14] = vec3(0.433, 0.250, 0.866);
  d[15] = vec3(0.866, -0.500, 0.000);
  d[16] = vec3(0.433, 0.250, 0.866);
  d[17] = vec3(0.000, -0.500, 0.866);
  d[18] = vec3(0.342, -0.940, 0.000);
  d[19] = vec3(0.866, -0.500, 0.000);
  d[20] = vec3(0.000, -0.500, 0.866);
  d[21] = vec3(-0.342, -0.940, 0.000);
  d[22] = vec3(0.342, -0.940, 0.000);
  d[23] = vec3(0.000, -0.500, 0.866);
  d[24] = vec3(-0.866, -0.500, 0.000);
  d[25] = vec3(-0.342, -0.940, 0.000);
  d[26] = vec3(0.000, -0.500, 0.866);
  d[27] = vec3(-0.866, -0.500, 0.000);
  d[28] = vec3(0.000, -0.500, 0.866);
  d[29] = vec3(-0.433, 0.250, 0.866);
  d[30] = vec3(-0.985, 0.174, 0.000);
  d[31] = vec3(-0.866, -0.500, 0.000);
  d[32] = vec3(-0.433, 0.250, 0.866);
  d[33] = vec3(-0.643, 0.766, 0.000);
  d[34] = vec3(-0.985, 0.174, 0.000);
  d[35] = vec3(-0.433, 0.250, 0.866);
  d[36] = vec3(-0.000, 1.000, 0.000);
  d[37] = vec3(-0.643, 0.766, 0.000);
  d[38] = vec3(-0.433, 0.250, 0.866);

  // Which vertex are we on?
  float sel = mod(id,vertices);
  vec3 z = vec3(0);
  vec3 vert = ((sel==0.0)?d[0]:z) + ((sel==1.0)?d[1]:z) + ((sel==2.0)?d[2]:z) + ((sel==3.0)?d[3]:z) + ((sel==4.0)?d[4]:z) + ((sel==5.0)?d[5]:z) + ((sel==6.0)?d[6]:z) + ((sel==7.0)?d[7]:z) + ((sel==8.0)?d[8]:z) + ((sel==9.0)?d[9]:z) + ((sel==10.0)?d[10]:z) + ((sel==11.0)?d[11]:z) + ((sel==12.0)?d[12]:z) + ((sel==13.0)?d[13]:z) + ((sel==14.0)?d[14]:z) + ((sel==15.0)?d[15]:z) + ((sel==16.0)?d[16]:z) + ((sel==17.0)?d[17]:z) + ((sel==18.0)?d[18]:z) + ((sel==19.0)?d[19]:z) + ((sel==20.0)?d[20]:z) + ((sel==21.0)?d[21]:z) + ((sel==22.0)?d[22]:z) + ((sel==23.0)?d[23]:z) + ((sel==24.0)?d[24]:z) + ((sel==25.0)?d[25]:z) + ((sel==26.0)?d[26]:z) + ((sel==27.0)?d[27]:z) + ((sel==28.0)?d[28]:z) + ((sel==29.0)?d[29]:z) + ((sel==30.0)?d[30]:z) + ((sel==31.0)?d[31]:z) + ((sel==32.0)?d[32]:z) + ((sel==33.0)?d[33]:z) + ((sel==34.0)?d[34]:z) + ((sel==35.0)?d[35]:z) + ((sel==36.0)?d[36]:z) + ((sel==37.0)?d[37]:z) + ((sel==38.0)?d[38]:z);

  // Get the position from the texture
  v_atomIndex = floor(id/vertices);
  vec2 texCoord = atomTextureCoords(v_atomIndex);
  vec4 position = texture2D(posTex, texCoord);

  // Is this atom selected?
  vec4 select = texture2D(selectTex, texCoord);
  // Use the mask to choose what component to use
  float ourId = dot(select, selectMask);
  // Background material is has select.z == 2.0
  v_background = (select.z == 2.0) ? 1.0 : 0.0;
  // Hide specific materials
  float show = (select.z == hideMaterial) ? 0.0 : 1.0;
  
  // Is this the selected atom?
  float selectFactor = (ourId == selectId) ? selectScale*(1.0-v_background) : 0.0;

  // Get the color and radius from the color texture
  vec4 colorRadius = texture2D(colorTex, texCoord);
  //vec4 selectColor = vec4(0, 0, selectFactor*(colorRadius.x+colorRadius.y+colorRadius.z), 0);
  //vec4 selectColor = vec4(selectFactor, selectFactor, (selectFactor>0.0)?1.0:0.0, 1);
  vec4 selectColor = vec4(selectFactor*vec3(0.3, 0.3, 1), 1);
  vec4 color = vec4(colorRadius.xyz, 1.0) + selectColor;
  // position.w is 1 if the atom is active, 0 if it is inactive
  float radius = 0.6*colorRadius.w*position.w*show; 

  // Get the surface position and normal.
  v_normal = vert; // Don't bother normalizing, since we will have to normalize after interpolation anyway.
  vec3 surfacePos = position.xyz + radius*vert;

  float bkg = 0.65*v_background;
  float gray = (color.x + color.y + color.z)/3.0;

  v_color = (1.0-bkg)*color.xyz + bkg*gray;
  v_surfaceToLight = lightPos - surfacePos;
  v_surfaceToView = cameraPos - surfacePos;
  gl_Position = matrix * vec4(surfacePos, 1.0);
}

`;

// Render the atoms using the full lighting model
this.drawSpheresFS = `
     precision highp float;

     uniform vec3 lightPos;
     uniform vec3 cameraPos;
     uniform float shininess;
     uniform float ambient;
     uniform float fadeDistance;
     varying vec3 v_normal;
     varying vec3 v_color;
     varying vec3 v_surfaceToView;
     varying vec3 v_surfaceToLight;
     varying float v_background;

     void main() {
       // Normalize the interpolated varying normal
       vec3 normal = normalize(v_normal);

       //float fadeFactor = fadeDistance*fadeDistance/dot(v_surfaceToLight,v_surfaceToLight);
       float fadeFactor = fadeDistance/length(v_surfaceToLight);
       vec3 surfaceToLightDirection = normalize(v_surfaceToLight);
       vec3 surfaceToViewDirection = normalize(v_surfaceToView);
       vec3 halfVector = normalize(surfaceToLightDirection + surfaceToViewDirection);

       float lightComp = dot(normal, surfaceToLightDirection);
       float lightDiffuse = lightComp*(1.0 - 0.5*v_background);

       float specular = 0.0;
       float specComp = dot(normal, halfVector);
       if (specComp > 0.0) {
	   specular = pow(specComp, shininess);
       }

       // Set the color
       gl_FragColor = vec4(v_color, 1);
       // Multiply just the color portion (not the alpha) by the light
       gl_FragColor.rgb *= fadeFactor*(lightDiffuse + ambient*(1.0 + 3.0*v_background));
       gl_FragColor.rgb += fadeFactor*specular*(1.0 - 0.5*v_background);
     }
`;


// Render atoms the color texture with no lighting model
this.drawFlatFS = `
     precision highp float;
     varying vec3 v_color;

     void main() {
        // Nothing fancy, just set the color
       gl_FragColor = vec4(v_color, 1);
     }
`;


// Render atoms with the positions placed in the x, y, z channels
// and the atom index plus 1 in the w channel
// This draws without any lighting model (for picking) and
this.drawPositionsFS = `
     precision highp float;
     uniform sampler2D posTex;
     varying float v_atomIndex;
     #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
     #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

     vec2 atomTextureCoords(float atomIndex) {
       float texX = (mod(atomIndex, ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_WIDTH;
       float texY = (floor(atomIndex/ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_HEIGHT;
       return vec2(texX, texY);
     }

     void main() {
       // Get the z coordinate from the position texture
       vec2 texCoord = atomTextureCoords(v_atomIndex);
       vec4 position = texture2D(posTex, texCoord);
       // Add one to the atom index since zero is no atom
       gl_FragColor = vec4(position.xyz, v_atomIndex + 1.0);
     }
`;


// This is a dummy vertex shader used by all the physics shaders
this.dummyVS = `
  attribute vec4 dummy;
  void main() {
    gl_Position = dummy;
  }
`;


// Cycle the linear congruential generators
this.randomFS = `
   precision highp float;
   uniform sampler2D randomTex;

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}
  // A linear congruential generator
  // We choose with small enough moduli (m) and multipliers (a) that 
  // exact integer math with 32-bit floats is possible
   #define RAND_M vec4(${this.RAND_M0}, ${this.RAND_M1}, ${this.RAND_M2}, ${this.RAND_M3})
   #define RAND_A vec4(${this.RAND_A0}, ${this.RAND_A1}, ${this.RAND_A2}, ${this.RAND_A3})

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // The position of this atom
     vec4 v = texture2D(randomTex, texCoord);
     gl_FragColor = mod(RAND_A*v, RAND_M);
   }
`;


// Reinitialize the velocities to a Maxwell-Boltzmann distribution
this.initVelFS = `
   precision highp float;
   uniform sampler2D randomTex;
   uniform sampler2D nonbondTex;
   uniform float kT;

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

   #define RAND_M vec4(${this.RAND_M0}, ${this.RAND_M1}, ${this.RAND_M2}, ${this.RAND_M3})
   #define PI 3.14159265358979323846

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // Nonbonded parameters for this atom (same dimensions as position)
     // {x: mass, y: charge, z: epsilon_LJ, w: Rmin_LJ}
     vec4 nonbondPar = texture2D(nonbondTex, texCoord);
     float mass = (nonbondPar.x==0.0)?1.0:nonbondPar.x;

     // The standard deviation of the velocity
     float vel0 = sqrt(kT/mass);

     // The random value 
     // Normalize the integers produced by the linear congruential generator 
     // (1 <= x <= m-1) => (1/m <= x' <= 1-1/m)
     vec4 random = texture2D(randomTex, texCoord)/RAND_M;
     // Box-Müller transform
     vec2 rad = sqrt(-2.0*log(random.xy));
     vec2 phi = 2.0*PI*random.zw;
     // We have a fourth normal deviate, but we don't use it
     // Different pattern to keep this uncorrelated with the Langevin thermostat
     vec3 gauss = vec3(rad.x*cos(phi.y), rad.x*sin(phi.y), rad.y*sin(phi.x));

     gl_FragColor = vec4(vel0*gauss, 1);
   }
`;


// Displace a selected group of atoms in a position texture
this.displaceFS = `
   precision highp float;
   uniform sampler2D posTex;
   uniform sampler2D selectTex;
   uniform vec4 selectMask;
   uniform float selectId;
   uniform vec3 displace;
   uniform vec3 box;

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

   // Wrap with the nearest image convention
   vec3 wrap(vec3 diffVec, vec3 box) {
     return diffVec - box*floor(diffVec/box + 0.5);
   }

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // Get the position and whether the atom is active
     vec4 posActive = texture2D(posTex, texCoord);

     // Is this atom selected?
     vec4 select = texture2D(selectTex, texCoord);
     float ourId = dot(select, selectMask);
     float selectFactor = (ourId == selectId) ? 1.0 : 0.0;

     gl_FragColor = vec4(wrap(posActive.xyz + selectFactor*displace, box), posActive.w);
     //gl_FragColor = vec4(posActive.xyz + selectFactor*displace, posActive.w);
   }
`;

    
// Displace a selected group of atoms so that a chosen
// atom is at a given position
this.displaceRelativeFS = `
   precision highp float;
   uniform sampler2D posTex;
   uniform sampler2D selectTex;
   uniform float selectId;
   uniform vec4 selectMask;
   uniform float pickedAtomIndex;
   uniform vec3 desiredPos;
   uniform vec3 axisMask;
   uniform vec3 box;

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

   // Wrap with the nearest image convention
   vec3 wrap(vec3 diffVec, vec3 box) {
     return diffVec - box*floor(diffVec/box + 0.5);
   }

   // Convert from the 1-dimensional atom index to the 2D atom texture coordinates
   vec2 atomTextureCoords(float atomIndex) {
     float texX = (mod(atomIndex, ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_WIDTH;
     float texY = (floor(atomIndex/ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_HEIGHT;
     return vec2(texX, texY);
   }

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);
     // Get the position and whether the atom is active
     vec4 posActive = texture2D(posTex, texCoord);

     // Get the position of the picked atom
     vec2 pickedTexCoord = atomTextureCoords(pickedAtomIndex);
     vec3 pickedPos = texture2D(posTex, pickedTexCoord).xyz;
     // This will displace the selected atoms so that the picked atom is at mouse's world position
     vec3 displace = axisMask*(desiredPos - pickedPos); 

     // Is this atom selected?
     vec4 select = texture2D(selectTex, texCoord);
     float ourId = dot(select, selectMask);
     float selectFactor = (ourId == selectId) ? 1.0 : 0.0;

     vec3 newPos = wrap(posActive.xyz + selectFactor*displace, box);
     gl_FragColor = vec4(newPos, posActive.w);
   }
`;


// Half-kick of velocity-Verlet time integration for NVE dynamics
// This shader is used for both the first and second half-kicks
// Izaguirre et al. (2001) DOI: 10.1063/1.1332996
this.verletHalfKickFS = `
   precision highp float;
   uniform sampler2D velTex;
   uniform sampler2D forceTex;
   uniform sampler2D nonbondTex;
   uniform float dt; // time step (in simulation units)

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // The current velocity
     vec3 vel = texture2D(velTex, texCoord).xyz;
     // The current force
     vec3 force = texture2D(forceTex, texCoord).xyz;

     // Nonbonded parameters for this atom (same dimensions as position)
     // {x: mass, y: charge, z: epsilon_LJ, w: Rmin_LJ}
     vec4 nonbondPar = texture2D(nonbondTex, texCoord);
     float mass = (nonbondPar.x==0.0)?1.0:nonbondPar.x;

     // Verlet time integration
     // Move the velocity forward one half step: v(t + dt/2)
     vec3 halfVel = vel + 0.5*dt/mass*force;
     gl_FragColor = vec4(halfVel, 1.0);
   }
`;


// Drift step of velocity-Verlet time integration
// This can be used for both NVE and Langevin thermostat (NVT)
// Equation 22 of Izaguirre et al. (2001) DOI: 10.1063/1.1332996
this.verletDriftFS = `
   precision highp float;
   uniform sampler2D velTex;
   uniform sampler2D posTex;
   uniform vec3 box;
   uniform float dt; // time step (in simulation units)

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}
   #define ATOM_MAX_DISPLACE ${this.ATOM_MAX_DISPLACE}

   // Wrap with the nearest image convention
   vec3 wrap(vec3 diffVec, vec3 box) {
     return diffVec - box*floor(diffVec/box + 0.5);
   }

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // The current position of this atom
     vec4 posActive = texture2D(posTex, texCoord);
     vec3 currPos = posActive.xyz;

     // The half velocity
     vec3 halfVel = texture2D(velTex, texCoord).xyz;

     // Verlet time integration
     // Get the next position
     vec3 nextPos = wrap(currPos + clamp(halfVel*dt, -ATOM_MAX_DISPLACE, ATOM_MAX_DISPLACE), box);
     //vec3 nextPos = wrap(currPos + halfVel*dt, box);
     //vec3 nextPos = currPos + halfVel*dt;
     gl_FragColor = vec4(nextPos, posActive.w);
   }
`;


// First half-kick of velocity-Verlet time integration for Langevin dynamics
// Equation 21 of Izaguirre et al. (2001) DOI: 10.1063/1.1332996
this.langevinFirstKickFS = `
   precision highp float;
   uniform sampler2D velTex;
   uniform sampler2D forceTex;
   uniform sampler2D nonbondTex;
   uniform float dt; // time step (in simulation units)

   // Langevin thermostat
   uniform sampler2D randomTex;
   uniform float langevinRatio; // 0.5*langevinGamma*dt (unitless)
   uniform float langevinRandom; // 2.0*langevinGamma*kT/dt (units of m*x^2/t^4)

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}
   #define RAND_M vec4(${this.RAND_M0}, ${this.RAND_M1}, ${this.RAND_M2}, ${this.RAND_M3})
   #define PI 3.14159265358979323846

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // The current velocity
     vec3 vel = texture2D(velTex, texCoord).xyz;
     // The current force
     vec3 force = texture2D(forceTex, texCoord).xyz;

     // The random value used for the random force
     // Normalize the integers produced by the linear congruential generator 
     // (1 <= x <= m-1) => (1/m <= x' <= 1-1/m)
     vec4 random = texture2D(randomTex, texCoord)/RAND_M;
     // Box-Müller transform
     vec2 rad = sqrt(-2.0*log(random.xy));
     vec2 phi = 2.0*PI*random.zw;
     // We have a fourth normal deviate, but we don't use it
     vec3 gauss = vec3(rad.x*cos(phi.x), rad.x*sin(phi.x), rad.y*cos(phi.y));

     // Nonbonded parameters for this atom (same dimensions as position)
     // {x: mass, y: charge, z: epsilon_LJ, w: Rmin_LJ}
     vec4 nonbondPar = texture2D(nonbondTex, texCoord);
     float mass = (nonbondPar.x==0.0)?1.0:nonbondPar.x;

     // Verlet time integration
     // Apply the Brooks-Brünger-Karplus algorithm for Langevin dynamics
     // This is the random force
     vec3 forceRandom = gauss*sqrt(mass*langevinRandom);

     // Velocity at the half step: v(t + dt/2)
     vec3 halfVel = (1.0 - langevinRatio)*vel + 0.5*dt/mass*(force + forceRandom);
     gl_FragColor = vec4(halfVel, 1.0);
   }
`;


// Second half-kick of velocity-Verlet time integration for Langevin dynamics
// Equation 23 of Izaguirre et al. (2001) DOI: 10.1063/1.1332996
this.langevinSecondKickFS = `
   precision highp float;
   uniform sampler2D velTex;
   uniform sampler2D forceTex;
   uniform sampler2D nonbondTex;
   uniform float dt; // time step (in simulation units)

   // Langevin thermostat
   uniform sampler2D randomTex;
   uniform float langevinRatio; // 0.5*langevinGamma*dt (unitless)
   uniform float langevinRandom; // 2.0*langevinGamma*kT/dt (units of m*x^2/t^4)

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}
   #define RAND_M vec4(${this.RAND_M0}, ${this.RAND_M1}, ${this.RAND_M2}, ${this.RAND_M3})
   #define PI 3.14159265358979323846

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // The current velocity
     vec3 vel = texture2D(velTex, texCoord).xyz;
     // The current force
     vec3 force = texture2D(forceTex, texCoord).xyz;

     // The random value used for the random force
     // Normalize the integers produced by the linear congruential generator 
     // (1 <= x <= m-1) => (1/m <= x' <= 1-1/m)
     vec4 random = texture2D(randomTex, texCoord)/RAND_M;
     // Box-Müller transform
     vec2 rad = sqrt(-2.0*log(random.xy));
     vec2 phi = 2.0*PI*random.zw;
     // We have a fourth normal deviate, but we don't use it
     vec3 gauss = vec3(rad.x*cos(phi.x), rad.x*sin(phi.x), rad.y*cos(phi.y));

     // Nonbonded parameters for this atom (same dimensions as position)
     // {x: mass, y: charge, z: epsilon_LJ, w: Rmin_LJ}
     vec4 nonbondPar = texture2D(nonbondTex, texCoord);
     float mass = (nonbondPar.x==0.0) ? 1.0 : nonbondPar.x;

     // Verlet time integration
     // Apply the Brooks-Brünger-Karplus algorithm for Langevin dynamics
     // This is the random force
     vec3 forceRandom = gauss*sqrt(mass*langevinRandom);

     // Velocity at the half step: v(t + dt/2)
     vec3 nextVel = (vel + 0.5*dt/mass*(force + forceRandom))/(1.0 + langevinRatio);
     gl_FragColor = vec4(nextVel, 1.0);
   }
`;


// Standard Verlet time integration
// x(t + dt) = 2*x(t) - x(t-dt) +  dt^2*force/mass
// This has poorer numerical properties than velocity-Verlet
// due to loss of precision in [x(t) - x(t-dt)]
this.verletFS = `
   precision highp float;
   uniform sampler2D prevPosTex;
   uniform sampler2D currPosTex;
   uniform sampler2D currForceTex;
   uniform sampler2D nonbondTex;
   uniform vec3 box;
   uniform float dt2; // squared time step (in simulation units)

   // Langevin thermostat parameters (not used in this shader)
   uniform sampler2D randomTex;
   uniform float langevinRatio;
   uniform float langevinRandom;

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}
   #define ATOM_MAX_DISPLACE ${this.ATOM_MAX_DISPLACE}

   // Wrap with the nearest image convention
   vec3 wrap(vec3 diffVec, vec3 box) {
     return diffVec - box*floor(diffVec/box + 0.5);
   }

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // The current position of this atom
     vec4 posActive = texture2D(currPosTex, texCoord);
     vec3 currPos = posActive.xyz;

     // The previous position 
     vec3 prevPos = texture2D(prevPosTex, texCoord).xyz;
     // The current force
     vec3 force = texture2D(currForceTex, texCoord).xyz;

     // Nonbonded parameters for this atom (same dimensions as position)
     // {x: mass, y: charge, z: epsilon_LJ, w: Rmin_LJ}
     vec4 nonbondPar = texture2D(nonbondTex, texCoord);
     float mass = (nonbondPar.x==0.0) ? 1.0 : nonbondPar.x;

     // Verlet time integration
     vec3 del = wrap(currPos - prevPos, box);
     vec3 displace = del + dt2*force/mass;
     vec3 nextPos = wrap(currPos + clamp(displace, -ATOM_MAX_DISPLACE, ATOM_MAX_DISPLACE), box);
     gl_FragColor = vec4(nextPos, posActive.w);
   }
`;


// Standard Verlet time integration of Langevin dynamics
// x(t + dt) = 2*x(t) - x(t-dt) +  dt^2*force/mass
// This has poorer numerical properties than velocity-Verlet
// due to loss of precision in [x(t) - x(t-dt)]
// Equation 16 of Phillips et al. (2005) DOI: 10.1002/jcc.20289
this.langevinFS = `
   precision highp float;
   uniform sampler2D prevPosTex;
   uniform sampler2D currPosTex;
   uniform sampler2D currForceTex;
   uniform sampler2D nonbondTex;
   uniform vec3 box;
   uniform float dt2; // squared time step (in simulation units)

   // Langevin thermostat
   uniform sampler2D randomTex;
   uniform float langevinRatio; // 0.5*langevinGamma*dt (unitless)
   uniform float langevinRandom; // 2.0*langevinGamma*kT/dt (units of m*x^2/t^4)

   #define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
   #define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}
   #define ATOM_MAX_DISPLACE ${this.ATOM_MAX_DISPLACE}
   #define RAND_M vec4(${this.RAND_M0}, ${this.RAND_M1}, ${this.RAND_M2}, ${this.RAND_M3})
   #define PI 3.14159265358979323846

   // Wrap with the nearest image convention
   vec3 wrap(vec3 diffVec, vec3 box) {
     return diffVec - box*floor(diffVec/box + 0.5);
   }

   void main() {
     // Get the texture coordinates for this atom in the position texture
     vec2 texCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

     // The current position of this atom
     vec4 posActive = texture2D(currPosTex, texCoord);
     vec3 currPos = posActive.xyz;

     // The previous position 
     vec3 prevPos = texture2D(prevPosTex, texCoord).xyz;
     // The current force
     vec3 force = texture2D(currForceTex, texCoord).xyz;

     // The random value used for the random force
     // Normalize the integers produced by the linear congruential generator 
     // (1 <= x <= m-1) => (1/m <= x' <= 1-1/m)
     vec4 random = texture2D(randomTex, texCoord)/RAND_M;
     // Box-Müller transform
     vec2 rad = sqrt(-2.0*log(random.xy));
     vec2 phi = 2.0*PI*random.zw;
     // We have a fourth normal deviate, but we don't use it
     vec3 gauss = vec3(rad.x*cos(phi.x), rad.x*sin(phi.x), rad.y*cos(phi.y));

     // Nonbonded parameters for this atom (same dimensions as position)
     // {x: mass, y: charge, z: epsilon_LJ, w: Rmin_LJ}
     vec4 nonbondPar = texture2D(nonbondTex, texCoord);
     float mass = nonbondPar.x;

     // Verlet time integration
     // Apply the Brooks-Brünger-Karplus algorithm for Langevin dynamics
     // This is the random force
     vec3 forceRandom = gauss*sqrt(mass*langevinRandom);
     vec3 del = wrap(currPos - prevPos, box);
     vec3 displace = ((1.0 - langevinRatio)*del + dt2*(force + forceRandom)/mass)/(1.0 + langevinRatio);
     vec3 nextPos = wrap(currPos + clamp(displace, -ATOM_MAX_DISPLACE, ATOM_MAX_DISPLACE), box);

     gl_FragColor = vec4(nextPos, posActive.w);
   }
`;


////////////////////////////////////////////////////////////////////////////////
// Provides all the force functions, which can be used in multiple shaders
this.forceEnergyPreamble = `
uniform sampler2D posTex;
uniform sampler2D nonbondTex;
uniform sampler2D bondTex;
uniform sampler2D angleTex;
uniform sampler2D dihedralTex;
uniform sampler2D excludeTex;
uniform sampler2D restrainTex;
uniform vec3 box;
uniform vec3 wall;
uniform float wallSpring;
uniform float coulomb;

#define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
#define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}
#define ATOM_MAX_DISPLACE ${this.ATOM_MAX_DISPLACE}
#define ATOM_MIN_DIST ${this.ATOM_MIN_DIST}
#define FORCE_SPLIT ${this.FORCE_SPLIT}
#define ENERGY_SPLIT ${this.ENERGY_SPLIT}

#define BOND_TEX_WIDTH ${this.BOND_TEX_WIDTH}
#define BOND_TEX_HEIGHT ${this.BOND_TEX_HEIGHT}
#define BOND_TERMS_PER_ATOM ${this.BOND_TERMS_PER_ATOM}
#define ANGLE_TEX_WIDTH ${this.ANGLE_TEX_WIDTH}
#define ANGLE_TEX_HEIGHT ${this.ANGLE_TEX_HEIGHT}
#define ANGLE_TERMS_PER_ATOM ${this.ANGLE_TERMS_PER_ATOM}
#define DIHEDRAL_TEX_WIDTH ${this.DIHEDRAL_TEX_WIDTH}
#define DIHEDRAL_TEX_HEIGHT ${this.DIHEDRAL_TEX_HEIGHT}
#define DIHEDRAL_TERMS_PER_ATOM ${this.DIHEDRAL_TERMS_PER_ATOM}
#define EXCLUDE_TEX_WIDTH ${this.EXCLUDE_TEX_WIDTH}
#define EXCLUDE_TEX_HEIGHT ${this.EXCLUDE_TEX_HEIGHT}
#define EXCLUDE_TERMS_PER_ATOM ${this.EXCLUDE_TERMS_PER_ATOM}

// Wrap with the nearest image convention
vec3 wrap(vec3 diffVec, vec3 box) {
  return diffVec - box*floor(diffVec/box + 0.5);
}

// Convert from the 1-dimensional atom index to the 2D atom texture coordinates
vec2 atomTextureCoords(float atomIndex) {
  float texX = (mod(atomIndex, ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_WIDTH;
  float texY = (floor(atomIndex/ATOM_TEX_WIDTH) + 0.5)/ATOM_TEX_HEIGHT;
  return vec2(texX, texY);
}


////////////////////////////////////////////////////////////////
// Forces for Lennard-Jones, Coulomb, and exclusions
// Here we handle the exclusions with a conditional
//
// The nonbonded parameter texture (nonbondTex) has the format:
// {x: mass, y: charge, z: epsilon_LJ, w: RMin_LJ}
// The exclusion texture (excludeTex) has the format:
// {x: partner atom index, y: unused, z: epsilon_LJ, w: RMin_LJ}
//
// epsilon_LJ = 0.0,  RMin_LJ = -1.0 is a real exclusion (LJ & Coulomb removed)
// epsilon_LJ = 0.0,  RMin_LJ = 0.0 is an inactive exclusion (does nothing)
// Other values represent special LJ (1-4 or bespoke); Coulomb should be normal 
vec4 nonbondForceExclude(vec4 posActive0, vec2 posTexCoord) {
  // Position of this atom
  vec3 pos0 = posActive0.xyz;
  float atomIndex0 = floor(gl_FragCoord.x) + floor(gl_FragCoord.y)*ATOM_TEX_WIDTH;

  // Nonbonded parameters for this atom (same dimensions as position)
  // {x: mass, y: charge, z: epsilon_LJ, w: RMin_LJ}
  vec4 nonbondPar0 = texture2D(nonbondTex, posTexCoord);

  float energy = 0.0;
  //float energyCoulomb = 0.0;
  vec3 force = vec3(0);
  //vec3 forceCoulomb = vec3(0);
     
  ////////////////////////////////////////////////////////////////
  // All pairs Coulomb and Lennard-Jones force calculation
  // We loop over all pixels in the position texture
  float excludeStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*EXCLUDE_TEX_WIDTH);
  for (float ix = 0.0; ix < ATOM_TEX_WIDTH; ix+=1.0) {
    for (float iy = 0.0; iy < ATOM_TEX_HEIGHT; iy+=1.0) {
      // Get the texture coordinates of the other atom
      float texX = (ix + 0.5) / ATOM_TEX_WIDTH;
      float texY = (iy + 0.5) / ATOM_TEX_HEIGHT;
      float atomIndex1 = ix + iy*ATOM_TEX_WIDTH;

      // Get the position of the other atom
      vec4 posActive1 = texture2D(posTex, vec2(texX, texY));
      vec3 pos1 = posActive1.xyz;

      // Get the nonbonded parameters of the other atom
      vec4 nonbondPar1 = texture2D(nonbondTex, vec2(texX, texY));
      // Lorentz-Berthelot mixing rules
      float epsMix = sqrt(nonbondPar0.z*nonbondPar1.z);
      float RMix = 0.5*(nonbondPar0.w + nonbondPar1.w);

      // Check for exclusions or special Lennard-Jones parameters
      float isExcluded = 0.0;
      for (float ei = 0.0; ei < EXCLUDE_TERMS_PER_ATOM; ei+=1.0) {
        float excludeTexX = excludeStart + ei/EXCLUDE_TEX_WIDTH;
        vec4 exclude = texture2D(excludeTex, vec2(excludeTexX, posTexCoord.y));
        float atomIndexExclude = exclude.x;
        float epsSpecial = exclude.z; 
        float RSpecial = exclude.w;

        // Check if this is a real exclusion (inactive exclusions are all zeros)
        // that affects these two atoms
        if (!(epsSpecial == 0.0 && RSpecial == 0.0) && atomIndexExclude == atomIndex1) {
          // Set special Lennard-Jones parameters
          // The conditional above ignores inactive exclusions (all zeros)
          epsMix = epsSpecial;
          RMix = RSpecial;
          // The signature for a real exclusion is eps=0.0, RMin=-1.0
          isExcluded = (epsSpecial == 0.0 && RSpecial == -1.0) ? 1.0 : isExcluded;
        }
      }

      // Don't include excluded atoms or self-interaction or inactive atoms
      if (isExcluded == 0.0 && atomIndex1 != atomIndex0 && posActive0.w != 0.0 && posActive1.w != 0.0) {
        // Get the separation vector and distance
        vec3 del = wrap(pos1 - pos0, box);
        float dist = length(del);

	float dist2 = dot(del,del);
	float dist3 = dist*dist2;
	float ratio = pow(RMix*RMix/dist2, 3.0);

	// Calculate the Coulomb and Lennard-Jones forces
	force += 12.0*epsMix*ratio*(1.0 - ratio)/dist2 * del;
        force += -coulomb*nonbondPar0.y*nonbondPar1.y/dist3 * del;
        // The 0.5 is because this energy is shared between two atoms
        energy += 0.5*(epsMix*ratio*(ratio - 2.0) + coulomb*nonbondPar0.y*nonbondPar1.y/dist);
      }
    }
  }
  return vec4(force, energy);
}


////////////////////////////////////////////////////////////////
// Forces for Lennard-Jones, Coulomb, and exclusions
// Here, excluded interactions are subtracted off in the second stage
// To avoid loss of precision due to large excluded forces,
// the forces are split into forceNeg, forceMed, forcePos
//
// The nonbonded parameter texture (nonbondTex) has the format:
// {x: mass, y: charge, z: epsilon_LJ, w: RMin_LJ}
// The exclusion texture (excludeTex) has the format:
// {x: partner atom index, y: unused, z: epsilon_LJ, w: RMin_LJ}
//
// epsilon_LJ = 0.0,  RMin_LJ = -1.0 is a real exclusion (LJ & Coulomb removed)
// epsilon_LJ = 0.0,  RMin_LJ = 0.0 is an inactive exclusion (does nothing)
// Other values represent special LJ (1-4 or bespoke); Coulomb should be normal 
vec4 nonbondForce(vec4 posActive0, vec2 posTexCoord) {
  // Position of this atom
  vec3 pos0 = posActive0.xyz;
  float atomIndex0 = floor(gl_FragCoord.x) + floor(gl_FragCoord.y)*ATOM_TEX_WIDTH;

  // Nonbonded parameters for this atom (same dimensions as position)
  // {x: mass, y: charge, z: epsilon_LJ, w: RMin_LJ}
  vec4 nonbondPar0 = texture2D(nonbondTex, posTexCoord);

  // Summing these separately typically reduces numerical error in the energy
  float energy12 = 0.0;
  float energy12Big = 0.0;
  float energy6 = 0.0;
  float energyCoulomb = 0.0;

  // Splitting the forces in into different sums, we can the error in the force by a factor of 2!
  vec3 forceNeg = vec3(0);
  vec3 forceMed = vec3(0);
  vec3 forcePos = vec3(0);
  vec3 forceCoulomb = vec3(0);
     
  ////////////////////////////////////////////////////////////////
  // All pairs Coulomb and Lennard-Jones force calculation
  // We loop over all pixels in the position texture
  // We do this blindly, ignoring exclusions (including special 1-4 and bespoke)
  // These forces are later subtracted below in these cases
  // Self-interactions are zero due to atomIndex comparison
  for (float ix = 0.0; ix < ATOM_TEX_WIDTH; ix+=1.0) {
    for (float iy = 0.0; iy < ATOM_TEX_HEIGHT; iy+=1.0) {
      // Get the texture coordinates of the other atom
      float texX = (ix + 0.5) / ATOM_TEX_WIDTH;
      float texY = (iy + 0.5) / ATOM_TEX_HEIGHT;
      float atomIndex1 = ix + iy*ATOM_TEX_WIDTH;

      // Get the position of the other atom
      vec4 posActive1 = texture2D(posTex, vec2(texX, texY));
      vec3 pos1 = posActive1.xyz;

      // Get the nonbonded parameters of the other atom
      vec4 nonbondPar1 = texture2D(nonbondTex, vec2(texX, texY));
      // Lorentz-Berthelot mixing rules
      float epsMix = sqrt(nonbondPar0.z*nonbondPar1.z);
      float RMix = 0.5*(nonbondPar0.w + nonbondPar1.w);

      // Get the separation vector and distance
      vec3 del = wrap(pos1 - pos0, box);
      float dist = length(del);

      // The conditional is to avoid infinities from self-interaction or inactive atoms
      if (atomIndex1 != atomIndex0 && posActive0.w != 0.0 && posActive1.w != 0.0) {
	float dist2 = dot(del,del);
	float dist3 = dist*dist2;
	float ratio = pow(RMix*RMix/dist2, 3.0);

	// Calculate the Coulomb and Lennard-Jones forces
	vec3 f = 12.0*epsMix*ratio*(1.0 - ratio)/dist2 * del;
        forceCoulomb -= coulomb*nonbondPar0.y*nonbondPar1.y/dist3 * del;

        // Splitting into different components significantly reduces numerical error
        if (f.x < -FORCE_SPLIT) {
           forceNeg.x += f.x;
        } else if (f.x > FORCE_SPLIT) {
           forcePos.x += f.x;
        } else {
           forceMed.x += f.x;
        }
        if (f.y < -FORCE_SPLIT) {
           forceNeg.y += f.y;
        } else if (f.y > FORCE_SPLIT) {
           forcePos.y += f.y;
        } else {
           forceMed.y += f.y;
        }
        if (f.z < -FORCE_SPLIT) {
           forceNeg.z += f.z;
        } else if (f.z > FORCE_SPLIT) {
           forcePos.z += f.z;
        } else {
           forceMed.z += f.z;
        }

	// We share this energy equally among the atoms, hence the 1/2
        // Splitting the energy into 12, 6, and Coulomb components improves numerics
        float e12 = 0.5*epsMix*ratio*ratio;
        if (e12 > ENERGY_SPLIT) {
           energy12Big += e12;
        } else {
           energy12 += e12;
        }
        energy6 += -epsMix*ratio;
        energyCoulomb += 0.5*coulomb*nonbondPar0.y*nonbondPar1.y/dist;
      }
    }
  }

  ////////////////////////////////////////////////////////////////
  // Exclusion calculation
  // This includes special 1-4 and bespoke LJ parameters
  float excludeStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*EXCLUDE_TEX_WIDTH);
  for (float ix = 0.0; ix < EXCLUDE_TERMS_PER_ATOM; ix+=1.0) {
    float excludeTexX = excludeStart + ix/EXCLUDE_TEX_WIDTH;
    vec4 exclude = texture2D(excludeTex, vec2(excludeTexX, posTexCoord.y));
    float atomIndex1 = exclude.x;
    float epsSpecial = exclude.z; 
    float RSpecial = exclude.w;

    // Get the position of the partner atom
    vec2 posTexCoord1 = atomTextureCoords(atomIndex1);
    vec4 posActive1 = texture2D(posTex, posTexCoord1);
    vec3 pos1 = posActive1.xyz;

    // The signature for a real exclusion is eps=0.0, RMin=-1.0
    float realExclude = (epsSpecial == 0.0 && RSpecial == -1.0) ? 1.0 : 0.0;

    // We need to recalculate the standard LJ so we can subtract it
    // Get the nonbonded parameters of the other atom
    vec4 nonbondPar1 = texture2D(nonbondTex, posTexCoord1);
    // Lorentz-Berthelot mixing rules
    float epsMix = sqrt(nonbondPar0.z*nonbondPar1.z);
    float RMix = 0.5*(nonbondPar0.w + nonbondPar1.w);
         
    // Get the separation vector and distance
    vec3 del = wrap(pos1 - pos0, box);
    float dist = length(del);

    // The conditional is to avoid infinities from self-interaction or inactive atoms
    // RSpecial checks if this an active exclusion, or just all zeros?
    if (RSpecial != 0.0 && posActive0.w != 0.0 && posActive1.w != 0.0) {
      float dist2 = dot(del,del);
      float dist3 = dist*dist2;
      float ratio = pow(RMix*RMix/dist2, 3.0);

      // Calculate the Coulomb and Lennard-Jones forces to subtract
      // The conditional is to avoid the self-interaction (or maybe clashes)
      vec3 fc = -12.0*epsMix*ratio*(1.0 - ratio)/dist2 * del;
      // Splitting into different components significantly reduces numerical error
      if (fc.x < -FORCE_SPLIT) {
          forceNeg.x += fc.x;
      } else if (fc.x > FORCE_SPLIT) {
          forcePos.x += fc.x;
      } else {
          forceMed.x += fc.x;
      }
      if (fc.y < -FORCE_SPLIT) {
         forceNeg.y += fc.y;
      } else if (fc.y > FORCE_SPLIT) {
         forcePos.y += fc.y;
      } else {
         forceMed.y += fc.y;
      }
      if (fc.z < -FORCE_SPLIT) {
         forceNeg.z += fc.z;
      } else if (fc.z > FORCE_SPLIT) {
         forcePos.z += fc.z;
      } else {
         forceMed.z += fc.z;
      }

      // If there is a bespoke or special LJ interaction apply it
      // For an empty exclusion, the force isn't modified since epsSpecial is zero
      float ratioSpecial = pow(RSpecial*RSpecial/dist2, 3.0);
      // If this is a real exclusion, we remove the Coulomb part too
      forceCoulomb += realExclude * coulomb*nonbondPar0.y*nonbondPar1.y/dist3 * del;

      // We remove the Lennard-Jones forces for this exclusion
      vec3 fSpecial = 12.0*epsSpecial*ratioSpecial*(1.0 - ratioSpecial)/dist2 * del;
      forceMed += fSpecial;

      // We share this energy equally among the atoms, hence the 1/2
      // Splitting the energy into 12, 6, and Coulomb components improves numerics
      float e12 = 0.5*epsMix*ratio*ratio;
      if (e12 > ENERGY_SPLIT) {
         energy12Big -= e12;
      } else {
         energy12 -= e12;
      }
      energy6 -= -epsMix*ratio;
      energyCoulomb -= realExclude*0.5*coulomb*nonbondPar0.y*nonbondPar1.y/dist;

      // Special energy
      energy12 += 0.5*epsSpecial*ratioSpecial*ratioSpecial;
      energy6 += -epsSpecial*ratioSpecial;
    }
  }

  return vec4(forcePos + forceNeg + forceMed + forceCoulomb, energy12 + energy6 + energyCoulomb);
}


////////////////////////////////////////////////////////////////
// Forces for Lennard-Jones, Coulomb, and exclusions
// Here, excluded interactions are subtracted off in the second stage
// Nothing is done to accelerate performance
//
// The nonbonded parameter texture (nonbondTex) has the format:
// {x: mass, y: charge, z: epsilon_LJ, w: RMin_LJ}
// The exclusion texture (excludeTex) has the format:
// {x: partner atom index, y: unused, z: epsilon_LJ, w: RMin_LJ}
//
// epsilon_LJ = 0.0,  RMin_LJ = -1.0 is a real exclusion (LJ & Coulomb removed)
// epsilon_LJ = 0.0,  RMin_LJ = 0.0 is an inactive exclusion (does nothing)
// Other values represent special LJ (1-4 or bespoke); Coulomb should be normal 
vec4 nonbondForceFast(vec4 posActive0, vec2 posTexCoord) {
  // Position of this atom
  vec3 pos0 = posActive0.xyz;
  float atomIndex0 = floor(gl_FragCoord.x) + floor(gl_FragCoord.y)*ATOM_TEX_WIDTH;

  // Nonbonded parameters for this atom (same dimensions as position)
  // {x: mass, y: charge, z: epsilon_LJ, w: RMin_LJ}
  vec4 nonbondPar0 = texture2D(nonbondTex, posTexCoord);

  // Summing these separately typically reduces numerical error in the energy
  float energy12 = 0.0;
  float energy6 = 0.0;
  float energyCoulomb = 0.0;

  // Splitting the forces in into different sums, we can the error in the force by a factor of 2!
  vec3 forceLJ = vec3(0);
  vec3 forceCoulomb = vec3(0);
     
  ////////////////////////////////////////////////////////////////
  // All pairs Coulomb and Lennard-Jones force calculation
  // We loop over all pixels in the position texture
  // We do this blindly, ignoring exclusions (including special 1-4 and bespoke)
  // These forces are later subtracted below in these cases
  // Self-interactions are zero due to atomIndex comparison
  for (float ix = 0.0; ix < ATOM_TEX_WIDTH; ix+=1.0) {
    for (float iy = 0.0; iy < ATOM_TEX_HEIGHT; iy+=1.0) {
      // Get the texture coordinates of the other atom
      float texX = (ix + 0.5) / ATOM_TEX_WIDTH;
      float texY = (iy + 0.5) / ATOM_TEX_HEIGHT;
      float atomIndex1 = ix + iy*ATOM_TEX_WIDTH;

      // Get the position of the other atom
      vec4 posActive1 = texture2D(posTex, vec2(texX, texY));
      vec3 pos1 = posActive1.xyz;

      // Get the nonbonded parameters of the other atom
      vec4 nonbondPar1 = texture2D(nonbondTex, vec2(texX, texY));
      // Lorentz-Berthelot mixing rules
      float epsMix = sqrt(nonbondPar0.z*nonbondPar1.z);
      float RMix = 0.5*(nonbondPar0.w + nonbondPar1.w);

      // Get the separation vector and distance
      vec3 del = wrap(pos1 - pos0, box);
      float dist = length(del);

      // The conditional is to avoid infinities from self-interaction or inactive atoms
      if (atomIndex1 != atomIndex0 && posActive0.w != 0.0 && posActive1.w != 0.0) {
	float dist2 = dot(del,del);
	float dist3 = dist*dist2;
	float ratio = pow(RMix*RMix/dist2, 3.0);

	// Calculate the Coulomb and Lennard-Jones forces
	vec3 f = 12.0*epsMix*ratio*(1.0 - ratio)/dist2 * del;
        forceCoulomb -= coulomb*nonbondPar0.y*nonbondPar1.y/dist3 * del;

        // Splitting into different components significantly reduces numerical error
        forceLJ += f;

	// We share this energy equally among the atoms, hence the 1/2
        // Splitting the energy into 12, 6, and Coulomb components improves numerics
        float e12 = 0.5*epsMix*ratio*ratio;
        energy12 += e12;
        energy6 += -epsMix*ratio;
        energyCoulomb += 0.5*coulomb*nonbondPar0.y*nonbondPar1.y/dist;
      }
    }
  }

  ////////////////////////////////////////////////////////////////
  // Exclusion calculation
  // This includes special 1-4 and bespoke LJ parameters
  float excludeStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*EXCLUDE_TEX_WIDTH);
  for (float ix = 0.0; ix < EXCLUDE_TERMS_PER_ATOM; ix+=1.0) {
    float excludeTexX = excludeStart + ix/EXCLUDE_TEX_WIDTH;
    vec4 exclude = texture2D(excludeTex, vec2(excludeTexX, posTexCoord.y));
    float atomIndex1 = exclude.x;
    float epsSpecial = exclude.z; 
    float RSpecial = exclude.w;

    // Get the position of the partner atom
    vec2 posTexCoord1 = atomTextureCoords(atomIndex1);
    vec4 posActive1 = texture2D(posTex, posTexCoord1);
    vec3 pos1 = posActive1.xyz;

    // The signature for a real exclusion is eps=0.0, RMin=-1.0
    float realExclude = (epsSpecial == 0.0 && RSpecial == -1.0) ? 1.0 : 0.0;

    // We need to recalculate the standard LJ so we can subtract it
    // Get the nonbonded parameters of the other atom
    vec4 nonbondPar1 = texture2D(nonbondTex, posTexCoord1);
    // Lorentz-Berthelot mixing rules
    float epsMix = sqrt(nonbondPar0.z*nonbondPar1.z);
    float RMix = 0.5*(nonbondPar0.w + nonbondPar1.w);
         
    // Get the separation vector and distance
    vec3 del = wrap(pos1 - pos0, box);
    float dist = length(del);

    // The conditional is to avoid infinities from self-interaction or inactive atoms
    // RSpecial checks if this an active exclusion, or just all zeros?
    if (RSpecial != 0.0 && posActive0.w != 0.0 && posActive1.w != 0.0) {
      float dist2 = dot(del,del);
      float dist3 = dist*dist2;
      float ratio = pow(RMix*RMix/dist2, 3.0);

      // Calculate the Coulomb and Lennard-Jones forces to subtract
      // The conditional is to avoid the self-interaction (or maybe clashes)
      vec3 fc = -12.0*epsMix*ratio*(1.0 - ratio)/dist2 * del;
      forceLJ += fc;

      // If there is a bespoke or special LJ interaction apply it
      // For an empty exclusion, the force isn't modified since epsSpecial is zero
      float ratioSpecial = pow(RSpecial*RSpecial/dist2, 3.0);
      // If this is a real exclusion, we remove the Coulomb part too
      forceCoulomb += realExclude * coulomb*nonbondPar0.y*nonbondPar1.y/dist3 * del;

      // We remove the Lennard-Jones forces for this exclusion
      vec3 fSpecial = 12.0*epsSpecial*ratioSpecial*(1.0 - ratioSpecial)/dist2 * del;
      forceLJ += fSpecial;

      // We share this energy equally among the atoms, hence the 1/2
      // Splitting the energy into 12, 6, and Coulomb components improves numerics
      float e12 = 0.5*epsMix*ratio*ratio;
      energy12 -= e12;
      energy6 -= -epsMix*ratio;
      energyCoulomb -= realExclude*0.5*coulomb*nonbondPar0.y*nonbondPar1.y/dist;

      // Special energy
      energy12 += 0.5*epsSpecial*ratioSpecial*ratioSpecial;
      energy6 += -epsSpecial*ratioSpecial;
    }
  }

  return vec4(forceLJ + forceCoulomb, energy12 + energy6 + energyCoulomb);
}


////////////////////////////////////////////////////////////////
// Calculate the bond force
// The bond texture (bondTex) has the format:
// {x: partner atom index, y: unused, z: bond energy constant, w: bond length}
//
// Each atom has multiple bonds, multiplying texture width by BOND_TERMS_PER_ATOM
// Textures are accessed on the domain 0 to 1
// The formula for the texture x coordinate of the bond, 
// bondStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*BOND_TEX_WIDTH),
// can be understood by the following diagram.
//
// Example: With 3 bond terms per atom
// posTex  |  o  |  o  |  o  |  o  | 
// bondTex |o| | |o| | |o| | |o| | |
vec4 bondForce(vec3 pos0, vec2 posTexCoord) {
  vec3 force = vec3(0.0);
  float energy = 0.0;

  float bondStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*BOND_TEX_WIDTH);
  for (float ix = 0.0; ix < BOND_TERMS_PER_ATOM; ix+=1.0) {
    // Get the coordinates of the bond parameter
    float bondTexX = bondStart + ix/BOND_TEX_WIDTH;
    vec4 bond = texture2D(bondTex, vec2(bondTexX, posTexCoord.y));
    float atomIndex1 = bond.x;
    float parK = bond.z;
    float parLen = bond.w;

    // Get the position of the partner atom
    vec2 posTexCoord1 = atomTextureCoords(atomIndex1);
    vec3 pos1 = texture2D(posTex, posTexCoord1).xyz;

    // Calculate the bond force
    vec3 del = wrap(pos1 - pos0, box);
    float dist = (parK==0.0)?1.0:length(del); // Avoid division by zero for inactive bonds
    // Inactive bonds have parK==0.0, so they add nothing to the force
    if (dist > ATOM_MIN_DIST) {
      float delDist = dist - parLen;
      force += 2.0*parK*delDist * del/dist;
      // Shared among the two atoms, so a factor of 1/2
      energy += 0.5*parK*delDist*delDist;
    }
  }

  return vec4(force, energy);
}

////////////////////////////////////////////////////////////////
// Calculate the angle force
// The angle texture (angleTex) has the format:
// - Two pixels for each angle term
// - {x: atom index 1, y: atom index 2, z: unused, w: unused}
// - {x: role (edge or central), y: energy constant, z: equil. angle, w: unused}
// - role 0 means the angle is ourAtom-atomIndex1-atomIndex2
// - role 1 means the angle is atomIndex1-ourAtom-atomIndex2
//
// Each atom has multiple angles, multiplying texture width by 2*ANGLE_TERMS_PER_ATOM
// Textures are accessed on the domain 0 to 1
// The formula for the texture x coordinate of the angle, 
// angleStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*ANGLE_TEX_WIDTH);
// can be understood by the following diagram.
//
// Example: With 3 angles (each with 2 pixels) per atom
// ANGLE_TEX_WIDTH = termsPerAtom * pixels * ATOM_TEX_WIDTH = 24
// posTex   |     o     |     o     |     o     |     o     | 
// angleTex |0| |0| |0| |0| |0| |0| |0| |0| |0| |0| |0| |0| | 
vec4 angleForce(vec3 pos0, vec2 posTexCoord) {
  vec3 force = vec3(0.0);
  float energy = 0.0;

  float angleStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*ANGLE_TEX_WIDTH);
  for (float ix = 0.0; ix < ANGLE_TERMS_PER_ATOM; ix+=1.0) {
    // Get the coordinates of the angle parameter
    float angleTexX0 = angleStart + 2.0*ix/ANGLE_TEX_WIDTH;
    vec4 angleAtoms = texture2D(angleTex, vec2(angleTexX0, posTexCoord.y));
    float angleTexX1 = angleTexX0 + 1.0/ANGLE_TEX_WIDTH;
    vec4 anglePar = texture2D(angleTex, vec2(angleTexX1, posTexCoord.y));

    // Get the parameters
    // Role: zero for edge atom, one for central atom
    float role = anglePar.x;
    float parK = anglePar.y;
    float parTheta = anglePar.z; // in radians

    // Get the atom positions
    float atomIndex1 = angleAtoms.x;
    vec2 posTexCoord1 = atomTextureCoords(atomIndex1);
    vec3 pos1 = texture2D(posTex, posTexCoord1).xyz;
    float atomIndex2 = angleAtoms.y;
    vec2 posTexCoord2 = atomTextureCoords(atomIndex2);
    vec3 pos2 = texture2D(posTex, posTexCoord2).xyz;

    // The angle geometry
    vec3 r1 = (role==0.0)?pos0:pos1;
    vec3 r2 = (role==0.0)?pos1:pos0;
    vec3 r3 = pos2;
    vec3 a = wrap(r1 - r2, box);
    vec3 b = wrap(r3 - r2, box);

    // Reciprocal magnitudes       
    // Avoid division by zero for inactive angles
    float a_mag = length(a);
    float b_mag = length(b);

    if (parK != 0.0 && a_mag > ATOM_MIN_DIST && b_mag > ATOM_MIN_DIST) {
      // Calculate the force 
      float a_rmag = 1.0/a_mag;
      float b_rmag = 1.0/b_mag;
      vec3 a_unit = a*a_rmag;
      vec3 b_unit = b*b_rmag;
      float cosTheta = dot(a_unit, b_unit);
      // The conditional can avoid some numerical problems
      float sinTheta = (abs(cosTheta) >= 1.0) ? 0.0 : sqrt(1.0 - cosTheta*cosTheta);
      float delTheta = acos(cosTheta) - parTheta;
      // Straight angles with very small sinTheta need this approximation for stability
      float factor = (sinTheta > 1e-5) ? (2.0*parK*delTheta/sinTheta) : sign(delTheta)*2.0*parK;

      // Calculate the angle force
      // Role 0 (leg atom) has a different formula than role 1 (central atom)
      if (role == 0.0) {
	force += factor * a_rmag * (b_unit - cosTheta*a_unit);
      } else {
	force += factor*(cosTheta*(a_unit*a_rmag + b_unit*b_rmag) - (a + b)*(a_rmag*b_rmag));
      }
      // Shared among the three atoms, so a factor of 1/3
      energy += parK*delTheta*delTheta/3.0;
    }
  }

  return vec4(force, energy);
}


////////////////////////////////////////////////////////////////
// Calculate the dihedral force
// The dihedral texture (dihedralTex) has the format:
// - Two pixels for each dihedral term
// - {x: atom index 1, y: atom index 2, z: atom index 3, w: unused}
// - {x: role (edge or central), y: energy constant, z: multiplicity, w: equil. angle}
// - role 0 means the dihedral is ourAtom-atomIndex1-atomIndex2-atomIndex3
// - role 1 means the dihedral is atomIndex1-ourAtom-atomIndex2-atomIndex3
//
// Each atom has multiple dihedrals, multiplying texture width by 2*DIHEDRAL_TERMS_PER_ATOM
// Textures are accessed on the domain 0 to 1
// The formula for the texture x coordinate of the dihedral, 
// dihedralStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*DIHEDRAL_TEX_WIDTH);
// can be understood by the following diagram.
//
// Example: With 3 dihedrals (each with 2 pixels) per atom
// DIHEDRAL_TEX_WIDTH = termsPerAtom * pixels * ATOM_TEX_WIDTH = 24
// posTex      |     o     |     o     |     o     |     o     | 
// dihedralTex |0| |0| |0| |0| |0| |0| |0| |0| |0| |0| |0| |0| | 
vec4 dihedralForce(vec3 pos0, vec2 posTexCoord) {
  vec3 force = vec3(0.0);
  float energy = 0.0;

  float dihedralStart = posTexCoord.x - 1.0/(2.0*ATOM_TEX_WIDTH) + 1.0/(2.0*DIHEDRAL_TEX_WIDTH);
  for (float ix = 0.0; ix < DIHEDRAL_TERMS_PER_ATOM; ix+=1.0) {
    // Get the coordinates of the dihedral parameter
    float dihedralTexX0 = dihedralStart + 2.0*ix/DIHEDRAL_TEX_WIDTH;
    vec4 dihedralAtoms = texture2D(dihedralTex, vec2(dihedralTexX0, posTexCoord.y));
    float dihedralTexX1 = dihedralTexX0 + 1.0/DIHEDRAL_TEX_WIDTH;
    vec4 dihedralPar = texture2D(dihedralTex, vec2(dihedralTexX1, posTexCoord.y));

    // Get the parameters
    // Role: zero for edge atom, one for central atom
    float role = dihedralPar.x;
    float parK = dihedralPar.y;
    float parMulti = dihedralPar.z;
    float parDelta = dihedralPar.w; // in radians

    // Get the atom positions
    float atomIndex1 = dihedralAtoms.x;
    vec2 posTexCoord1 = atomTextureCoords(atomIndex1);
    vec3 pos1 = texture2D(posTex, posTexCoord1).xyz;
    float atomIndex2 = dihedralAtoms.y;
    vec2 posTexCoord2 = atomTextureCoords(atomIndex2);
    vec3 pos2 = texture2D(posTex, posTexCoord2).xyz;
    float atomIndex3 = dihedralAtoms.z;
    vec2 posTexCoord3 = atomTextureCoords(atomIndex3);
    vec3 pos3 = texture2D(posTex, posTexCoord3).xyz;

    // The dihedral geometry
    vec3 r1 = (role==0.0)?pos0:pos1;
    vec3 r2 = (role==0.0)?pos1:pos0;
    vec3 r3 = pos2;
    vec3 r4 = pos3;
    vec3 a = wrap(r2 - r1, box);
    vec3 b = wrap(r3 - r2, box);
    vec3 c = wrap(r4 - r3, box);
    vec3 r13 = wrap(r3 - r1, box);
       
    // Calculate the needed cross products
    vec3 M0 = cross(a, b);
    vec3 N0 = cross(b, c);
    float M0_mag = length(M0);
    float N0_mag = length(N0);

    // Avoid division by zero for inactive dihedrals and colinear systems
    if (parK != 0.0 &&  M0_mag > ATOM_MIN_DIST && N0_mag > ATOM_MIN_DIST) {
      float M0_rmag = 1.0/M0_mag;
      float N0_rmag = 1.0/N0_mag;
      vec3 M = M0_rmag*M0;
      vec3 N = N0_rmag*N0;
      vec3 M_cross_b = cross(M, b);
      vec3 N_cross_b = cross(N, b);
      vec3 M_cross_c = cross(M, c);
      vec3 N_cross_c = cross(N, c);
      vec3 M_cross_r13 = cross(M, r13);
      vec3 N_cross_r13 = cross(N, r13);

      // Calculate the dihedral angle
      float cosPhi = dot(M, N);
      // The conditional is to avoid erroneous phi values (pi/4) on some renderers
      float sinPhi = (abs(cosPhi) >= 1.0) ? 0.0 : sqrt(1.0 - cosPhi*cosPhi);
      float phi = -atan(sinPhi, cosPhi);
      // Impropers have parMulti = 0
      // Surprisingly, there is little change between regular and improper dihedrals, 
      // since most of the work goes into applying the chain rule to the angle.
      float factor = (parMulti==0.0) ? -2.0*parK*(phi - parDelta) : parMulti*parK*sin(parMulti*phi - parDelta);

      // The derivative of the potential energy is
      // D_xi(V) = -parK*sin(n*phi - delta)*D_xi(phi)
      // There are two ways to derive D_xi(phi):
      // (1) phi = arccos( dot(M,N) / (|M|*|N|) )
      // (2) phi = arcsin( |cross(M,N)| / (|M|*|N|) )
      // The first results in a 1/sin(phi)
      // The first becomes numerically unstable for small sin(phi)
      // We have derived the result for the second case, but
      // the first works out fine if we just avoid zero, it seems
      
      if (sinPhi != 0.0) {
	// Solved with SymPy
	force += (role==0.0) ? factor/sinPhi*M0_rmag*(N_cross_b - cosPhi*M_cross_b) : factor/sinPhi*(N0_rmag*M_cross_c + cosPhi*(M0_rmag*M_cross_r13 - N0_rmag*N_cross_c) - M0_rmag*N_cross_r13);
      }
      // We share the energy among 4 atoms, hence the factor of 1/4
      energy += (parMulti==0.0) ? 0.25*parK*(phi-parDelta)*(phi-parDelta) : 0.25*parK*(1.0 + cos(parMulti*phi - parDelta));
    } // End active if
  } // End for over dihedrals (active and inactive)
  
  return vec4(force, energy);
}


// Calculate the restraint force
// The restraint texture (restrainTex) has the format:
// {x: x_reference, y: y_ref, z: z_ref, w: force constant}
vec4 restraintForce(vec3 pos0, vec2 posTexCoord) {
 // The original restrained position is orig.xyz
  // The spring constant is orig.w 
  vec4 orig = texture2D(restrainTex, posTexCoord);
  vec3 del = wrap(orig.xyz - pos0, box);
  vec3 restraintForce = orig.w*del;
  float restraintEnergy = 0.5*orig.w*dot(del,del);
  return vec4(restraintForce, restraintEnergy);
}

// Calculate the force of flat-bottomed harmonic walls
// Restrained atoms are ignored
// The restraint texture (restrainTex) has the format:
// {x: x_reference, y: y_ref, z: z_ref, w: force constant}
vec4 wallForce(vec3 pos0, vec2 posTexCoord) {
  // Restrained atoms don't feel the wall force
  // The spring constant is orig.w 
  vec4 orig = texture2D(restrainTex, posTexCoord);
  float unrestrained = (orig.w > 0.0) ? 0.0 : 1.0;

  vec3 delWall = max(vec3(0.0), pos0 - wall) + min(vec3(0.0), pos0 + wall);
  vec3 wallForce = -unrestrained*wallSpring*delWall;
  float wallEnergy =  0.5*unrestrained*wallSpring*dot(delWall,delWall);
  return vec4(wallForce, wallEnergy);    
}

// Sum the force and energy of all components
vec4 computeForceEnergy(vec4 posActive0, vec2 posTexCoord) {
   vec3 pos0 = posActive0.xyz;

  // Initialize the force to zero
  vec4 forceEnergy = vec4(0.0);

  // Lennard-Jones, Coulomb, exclusions
  //forceEnergy += nonbondForceExclude(posActive0, posTexCoord); 
  forceEnergy += nonbondForce(posActive0, posTexCoord); 
  //forceEnergy += nonbondForceFast(posActive0, posTexCoord); 
  // Bonds
  forceEnergy += bondForce(pos0, posTexCoord);
  // Angles
  forceEnergy += angleForce(pos0, posTexCoord);
  // Dihedral force
  forceEnergy += dihedralForce(pos0, posTexCoord);

  // Restraints
  // Atomic restraints
  forceEnergy += restraintForce(pos0, posTexCoord);
  // Wall force
  forceEnergy += wallForce(pos0, posTexCoord);

  return forceEnergy;
}
`;


////////////////////////////////////////////////////////////////////////////////
// Calculate the force and potential energy using functions in the preamble
this.forceEnergyFS = `
precision highp float;
${this.forceEnergyPreamble}

void main() {
  // Get the texture coordinates for this atom in the position texture
  vec2 posTexCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);
  // The position of this atom
  vec4 posActive0 = texture2D(posTex, posTexCoord);
  vec3 pos0 = posActive0.xyz;

  // Compute the force and energy
  vec4 forceEnergy = computeForceEnergy(posActive0, posTexCoord);

  // Return the force and energy
  gl_FragColor = forceEnergy;
}
`;

// Move down the potential gradient for a steepest descent energy minimizer
this.moveDownhillFS = `
precision highp float;
uniform sampler2D posTex;
uniform sampler2D forceTex;
uniform float stepSize;
uniform vec3 box;

#define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
#define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

// Wrap with the nearest image convention
vec3 wrap(vec3 diffVec, vec3 box) {
  return diffVec - box*floor(diffVec/box + 0.5);
}

void main() {
  // Get the texture coordinates for this atom in the position texture
  vec2 posTexCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);
  // The position of this atom
  vec4 posActive = texture2D(posTex, posTexCoord);

  // Get the force on this atom
  vec3 force = texture2D(forceTex, posTexCoord).xyz;

  // Shift downhill
  float forceSq = dot(force, force);
  vec3 dir = (forceSq == 0.0) ? vec3(0) : normalize(force);
  vec3 nextPos = wrap(posActive.xyz + stepSize*dir, box);

  // Return the new position
  gl_FragColor = vec4(nextPos, posActive.w);
}
`;

// Just copy one texture to another
this.copyFS = `
precision highp float;
uniform sampler2D posTex;

#define ATOM_TEX_WIDTH ${this.ATOM_TEX_WIDTH}
#define ATOM_TEX_HEIGHT ${this.ATOM_TEX_HEIGHT}

void main() {
  // Get the texture coordinates for this atom in the position texture
  vec2 posTexCoord = gl_FragCoord.xy / vec2(ATOM_TEX_WIDTH, ATOM_TEX_HEIGHT);

  // Get the position of this atom
  vec4 posActive = texture2D(posTex, posTexCoord);

  // Return the same position
  gl_FragColor = posActive;
}
`;


////////////////////////////////////////////////////////////////////////////////
// Extract floating point values from a floating-point texture
// The destination texture has 4 times the source texture width since
// the floating point values are written as bytes 
this.extractFS = `
precision highp float;
uniform sampler2D dataTex;
uniform vec2 destTexDim;

// The 5 functions below are from StackOverflow user Adrian Seeley
// https://stackoverflow.com/questions/17981163/webgl-read-pixels-from-floating-point-render-target
float shift_right (float v, float amt) { 
  v = floor(v) + 0.5; 
  return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
  return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
  return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
  from = floor(from + 0.5); to = floor(to + 0.5); 
  return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
  if (val == 0.0) return vec4(0, 0, 0, 0); 
  float sign = val > 0.0 ? 0.0 : 1.0; 
  val = abs(val); 
  float exponent = floor(log2(val)); 
  float biased_exponent = exponent + 127.0; 
  float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
  float t = biased_exponent / 2.0; 
  float last_bit_of_biased_exponent = fract(t) * 2.0; 
  float remaining_bits_of_biased_exponent = floor(t); 
  float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
  float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
  float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
  float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
  return vec4(byte4, byte3, byte2, byte1); 
}

void main() {
  // gl_FragCoord is in pixel coordinates of the destination
  // but the textures are referenced on 0 to 1
  // For this reason, we don't actually need the source dimensions
  vec2 texCoord = gl_FragCoord.xy/destTexDim;
  vec4 value = texture2D(dataTex, texCoord);
  float r = floor(mod(gl_FragCoord.x, 4.0));
  float c = (r==0.0)?value.x:((r==1.0)?value.y:((r==2.0)?value.z:value.w));
  gl_FragColor = encode_float(c);
}
`;

// Compile the programs and get variable locations
this.preparePrograms(gl);    
} // end of constructor


    // Begin normal identation
    preparePrograms(gl) {
	/////////////////////////////////////////////////////////
	// Create the shader programs
	/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////
	// Drawing programs
	/////////////////////////////////////////////////////////

	// Draw the atoms as spheres with shading
	this.drawSpheres = webglUtils.createProgramFromSources(
	    gl, [this.drawSpheresVS, this.drawSpheresFS]);
	this.drawSpheresLocs = {
	    id: gl.getAttribLocation(this.drawSpheres, 'id'),
	    posTex: gl.getUniformLocation(this.drawSpheres, 'posTex'),
	    colorTex: gl.getUniformLocation(this.drawSpheres, 'colorTex'),
	    selectTex: gl.getUniformLocation(this.drawSpheres, 'selectTex'),
	    selectId: gl.getUniformLocation(this.drawSpheres, 'selectId'),
	    selectMask: gl.getUniformLocation(this.drawSpheres, 'selectMask'),
	    selectScale: gl.getUniformLocation(this.drawSpheres, 'selectScale'),
	    hideMaterial: gl.getUniformLocation(this.drawSpheres, 'hideMaterial'),
	    matrix: gl.getUniformLocation(this.drawSpheres, 'matrix'),
	    shininess: gl.getUniformLocation(this.drawSpheres, 'shininess'),
	    ambient: gl.getUniformLocation(this.drawSpheres, 'ambient'),
	    fadeDistance: gl.getUniformLocation(this.drawSpheres, 'fadeDistance'),
	    cameraPos: gl.getUniformLocation(this.drawSpheres, 'cameraPos'),
	    lightPos: gl.getUniformLocation(this.drawSpheres, 'lightPos'),
	};

	// Draw the atoms with flat shading colored by position
	// (x, y, z, atomIndex+1)
	this.drawPositions = webglUtils.createProgramFromSources(
	    gl, [this.drawSpheresVS, this.drawPositionsFS]);
	this.drawPositionsLocs = {
	    id: gl.getAttribLocation(this.drawPositions, 'id'),
	    posTex: gl.getUniformLocation(this.drawPositions, 'posTex'),
	    colorTex: gl.getUniformLocation(this.drawPositions, 'colorTex'),
	    selectTex: gl.getUniformLocation(this.drawPositions, 'selectTex'),
	    selectId: gl.getUniformLocation(this.drawPositions, 'selectId'),
	    selectMask: gl.getUniformLocation(this.drawPositions, 'selectMask'),
	    hideMaterial: gl.getUniformLocation(this.drawPositions, 'hideMaterial'),
	    selectScale: gl.getUniformLocation(this.drawPositions, 'selectScale'),
	    matrix: gl.getUniformLocation(this.drawPositions, 'matrix'),
	    shininess: gl.getUniformLocation(this.drawPositions, 'shininess'),
	    ambient: gl.getUniformLocation(this.drawPositions, 'ambient'),
	    fadeDistance: gl.getUniformLocation(this.drawPositions, 'fadeDistance'),
	    cameraPos: gl.getUniformLocation(this.drawPositions, 'cameraPos'),
	    lightPos: gl.getUniformLocation(this.drawPositions, 'lightPos'),
	};


	/////////////////////////////////////////////////////////
	// Molecular dynamics programs
	//
	// References for Velocity Verlet and Brooks-Brünger-Karplus Langevin dynamics schemes
	// Phillips et al. (2005) Scalable Molecular Dynamics with NAMD. J. Comput. Chem. 26(16), 1781–1802. DOI: 10.1002/jcc.20289
	// Izaguirre et al. (2001) Langevin stabilization of molecular dynamics. J. Chem. Phys. 114(5), 2090–2098. DOI: 10.1063/1.1332996
	// Brooks et al. (1984) Stochastic boundary conditions for molecular dynamics simulations of ST2 water. Chem. Phys. Lett. 105(5), 495–500. DOI: 10.1016/0009-2614(84)80098-6
	//
	// Velocity Verlet drift stage
	this.drift = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.verletDriftFS]);
	this.driftLocs = {
	    dummy: gl.getAttribLocation(this.drift, 'dummy'),
	    velTex: gl.getUniformLocation(this.drift, 'velTex'),
	    posTex: gl.getUniformLocation(this.drift, 'posTex'),
	    box: gl.getUniformLocation(this.drift, 'box'),
	    dt: gl.getUniformLocation(this.drift, 'dt'),
	};

	// Velocity Verlet half-kick stage
	this.verletKick = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.verletHalfKickFS]);
	this.verletKickLocs = {
	    dummy: gl.getAttribLocation(this.verletKick, 'dummy'),
	    velTex: gl.getUniformLocation(this.verletKick, 'velTex'),
	    forceTex: gl.getUniformLocation(this.verletKick, 'forceTex'),
	    nonbondTex: gl.getUniformLocation(this.verletKick, 'nonbondTex'),
	    dt: gl.getUniformLocation(this.verletKick, 'dt'),
	};

	// Langevin velocity Verlet first half-kick stage
	this.firstKick = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.langevinFirstKickFS]);
	this.firstKickLocs = {
	    dummy: gl.getAttribLocation(this.firstKick, 'dummy'),
	    velTex: gl.getUniformLocation(this.firstKick, 'velTex'),
	    forceTex: gl.getUniformLocation(this.firstKick, 'forceTex'),
	    nonbondTex: gl.getUniformLocation(this.firstKick, 'nonbondTex'),
	    dt: gl.getUniformLocation(this.firstKick, 'dt'),
	    // Langevin thermostat variables
	    randomTex: gl.getUniformLocation(this.firstKick, 'randomTex'),
	    langevinRatio: gl.getUniformLocation(this.firstKick, 'langevinRatio'),
	    langevinRandom: gl.getUniformLocation(this.firstKick, 'langevinRandom'),
	};

	// Langevin velocity Verlet second half-kick stage
	this.secondKick = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.langevinSecondKickFS]);
	this.secondKickLocs = {
	    dummy: gl.getAttribLocation(this.secondKick, 'dummy'),
	    velTex: gl.getUniformLocation(this.secondKick, 'velTex'),
	    forceTex: gl.getUniformLocation(this.secondKick, 'forceTex'),
	    nonbondTex: gl.getUniformLocation(this.secondKick, 'nonbondTex'),
	    dt: gl.getUniformLocation(this.secondKick, 'dt'),
	    // Langevin thermostat variables
	    randomTex: gl.getUniformLocation(this.secondKick, 'randomTex'),
	    langevinRatio: gl.getUniformLocation(this.secondKick, 'langevinRatio'),
	    langevinRandom: gl.getUniformLocation(this.secondKick, 'langevinRandom'),
	};

	// Program for calculating forces from atom positions
	this.force = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.forceEnergyFS]);
	this.forceLocs = {
	    dummy: gl.getAttribLocation(this.force, 'dummy'),
	    bondTex: gl.getUniformLocation(this.force, 'bondTex'),
	    angleTex: gl.getUniformLocation(this.force, 'angleTex'),
	    dihedralTex: gl.getUniformLocation(this.force, 'dihedralTex'),
	    excludeTex: gl.getUniformLocation(this.force, 'excludeTex'),
	    nonbondTex: gl.getUniformLocation(this.force, 'nonbondTex'),
	    restrainTex: gl.getUniformLocation(this.force, 'restrainTex'),
	    posTex: gl.getUniformLocation(this.force, 'posTex'),
	    box: gl.getUniformLocation(this.force, 'box'),
	    wall: gl.getUniformLocation(this.force, 'wall'),
	    wallSpring: gl.getUniformLocation(this.force, 'wallSpring'),
	    coulomb: gl.getUniformLocation(this.force, 'coulomb'),
	};

	// Program for moving downhill with the force
	this.downhill = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.moveDownhillFS]);
	this.downhillLocs = {
	    dummy: gl.getAttribLocation(this.downhill, 'dummy'),
	    posTex: gl.getUniformLocation(this.downhill, 'posTex'),
	    forceTex: gl.getUniformLocation(this.downhill, 'forceTex'),
	    stepSize: gl.getUniformLocation(this.downhill, 'stepSize'),
	    box: gl.getUniformLocation(this.downhill, 'box'),
	};

	// Just copy the texture to a new texture
	this.copy = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.copyFS]);
	this.copyLocs = {
	    dummy: gl.getAttribLocation(this.copy, 'dummy'),
	    posTex: gl.getUniformLocation(this.copy, 'posTex'),
	};
	
	// Program for displacing atoms
	this.displace = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.displaceFS]);
	this.displaceLocs = {
	    dummy: gl.getAttribLocation(this.displace, 'dummy'),
	    posTex: gl.getUniformLocation(this.displace, 'posTex'),
	    selectTex: gl.getUniformLocation(this.displace, 'selectTex'),
	    selectId: gl.getUniformLocation(this.displace, 'selectId'),
	    selectMask: gl.getUniformLocation(this.displace, 'selectMask'),
	    displace: gl.getUniformLocation(this.displace, 'displace'),
	    box: gl.getUniformLocation(this.displace, 'box'),
	};

	// Program for displacing atoms relative to another atom's position
	// All atoms of the selection will be moved so that the picked atom is at desiredPos
	// The axis mask, for example, (1, 1, 0), allows the displacement to only happen
	// along chosen axes (x and y in the example)
	this.displaceRelative = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.displaceRelativeFS]);
	this.displaceRelativeLocs = {
	    dummy: gl.getAttribLocation(this.displaceRelative, 'dummy'),
	    posTex: gl.getUniformLocation(this.displaceRelative, 'posTex'),
	    selectTex: gl.getUniformLocation(this.displaceRelative, 'selectTex'),
	    selectId: gl.getUniformLocation(this.displaceRelative, 'selectId'),
	    selectMask: gl.getUniformLocation(this.displaceRelative, 'selectMask'),
	    pickedAtomIndex: gl.getUniformLocation(this.displaceRelative, 'pickedAtomIndex'),
	    desiredPos: gl.getUniformLocation(this.displaceRelative, 'desiredPos'),
	    axisMask: gl.getUniformLocation(this.displaceRelative, 'axisMask'),
	    box: gl.getUniformLocation(this.displaceRelative, 'box')
	};

	// Program for resetting the velocities
	this.initVel = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.initVelFS]);
	this.initVelLocs = {
	    dummy: gl.getAttribLocation(this.initVel, 'dummy'),
	    randomTex: gl.getUniformLocation(this.initVel, 'randomTex'),
	    nonbondTex: gl.getUniformLocation(this.initVel, 'nonbondTex'),
	    kT: gl.getUniformLocation(this.initVel, 'kT'),
	};

	// Program for advancing the random numbers
	this.random = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.randomFS]);
	this.randomLocs = {
	    dummy: gl.getAttribLocation(this.random, 'dummy'),
	    randomTex: gl.getUniformLocation(this.random, 'randomTex'),
	};

	// Program for extracting data from textures
	this.extract = webglUtils.createProgramFromSources(gl, [this.dummyVS, this.extractFS]);
	this.extractLocs = {
	    dummy: gl.getAttribLocation(this.extract, 'dummy'),
	    dataTex: gl.getUniformLocation(this.extract, 'dataTex'),
	    destTexDim: gl.getUniformLocation(this.extract, 'destTexDim'),
	};
    }

    
} // end of MDShaders class
