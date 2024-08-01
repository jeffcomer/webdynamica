#!/bin/bash

tab="    "
moleculeList=(../../Molecules/mol_{water,benzene,benzoate,benzamidine,NAFT,methane,octane,dodecane,glucose,N2,O2,argon,graphene,sodium_ion,chloride_ion,potassium,calcium,magnesium,bisulfate,cellulose_fragment,lignin_fragment,quartz,SiO2}.webdyn.js ../../Molecules/mol_zwit_aa_{A,E,F,G,K,R,S,T,V,W,Y}.webdyn.js ../../Molecules/mol_zwit_PRO.webdyn.js)

# Insert the molecules if they exist
if [[ -f ${moleculeList[0]} ]]; then
    # Make a copy of index.html
    cp index.html index.html.old
    input=index.html.old
    output=index.html
    
    # Write everything up to LOAD_MOLECULES
    awk '{print $0} /<!--BEGIN_LOAD_MOLECULES-->/ {exit};' $input > $output

    nameStr="const moleculeNames = ["
    for f in ${moleculeList[@]}
    do
	fn=$( basename $f )
	pre=${fn%.webdyn.js}
	name=${pre#mol_}
	
	echo "${tab}<script>" >> $output
	cat $f >> $output
	echo "" >> $output
	echo "${tab}</script>" >> $output
	nameStr="${nameStr}'${name}',"
    done

     # Write everything after END_LOAD_MOLECULES
    awk 'BEGIN {ok=0}; /<!--END_LOAD_MOLECULES-->/ {ok=1}; ok {print $0};' $input >> $output
    nameStr="${nameStr}];"
    echo $nameStr
else
    echo "File ${moleculeList[0]} does not exist. Not modifying index.html"
fi

# iPhone 15 has a 6.1-inch screen with a screen size (resolution): 1179px × 2556px , 393px × 852px viewport 1, and a CSS Pixel Ratio of 3. 
google-chrome --window-size=393,852 --new-window index.html &


