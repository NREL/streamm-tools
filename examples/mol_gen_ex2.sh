#!/bin/bash

# Example 2
# Generate a gaussian input file for oligo methylthiophene 

python ${TOOLS_PATH}/da_builder/donoracceptorsystems.py  "thiophene ( R_methane  ) " -b ${BB_PATH}  -r 5

# will generate a rigioregular oligomer of methylthiol n=1-5 

python ${TOOLS_PATH}/da_builder/donoracceptorsystems.py  "thiophene ( R_methane  ) " -b ${BB_PATH}  -r 5  -p "0"

# sets the inter-ring dihedral angle to zero making the sulfurs of thiophene in the cis configuration 

python ${TOOLS_PATH}/da_builder/donoracceptorsystems.py  "thiophene ( R_methane ) " -b ${BB_PATH}  -r 5  -p "180 0 "

# sets the inter-ring dihedral angle to alternate between sero and 180 making the sulfurs of thiophene in the trans configuration 
