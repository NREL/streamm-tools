#!/bin/bash


source mol_gen_setup.sh 


# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

python ${TOOLS_PATH}/da_builder/donoracceptorsystems.py -h

# to see a list of the commands 

# Generate a gaussian input file for thiophene by running     

python ${TOOLS_PATH}/da_builder/donoracceptorsystems.py  "thiophene ( R_hydrogen  ) " -b ${BB_PATH}  -r 1 

# this creates the file
echo " The structure file "
echo "    mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.xyz "
echo "  has been generated "
more mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.xyz

# which is just a file containing the cartesian coordinates of a
# thiophene molecule.  You can view with your favorite viewer. The -r option is set to 1 to generate a single molecule rather than an oligomer 
 
python ${TOOLS_PATH}/da_builder/donoracceptorsystems.py  "thiophene ( R_hydrogen  ) " -b ${BB_PATH}  -r 1 

# creates the file
#    mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.com    
# which is the gaussian input file based on "donoracceptor.com.template"
#    mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.pbs    
# which is the pbs script file based on "donoracceptor.pbs.template"
#    mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.json   
# which is the json file contain the molecular information 
#    mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.meta   
#    mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.r1.com 
#    mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.r2.com

# These templates can be modified for your HPC system or what ever properties of the molecule you wish to calculate 