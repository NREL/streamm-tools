#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build C60
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

# Generate structure files for C60
donoracceptorsystems.py  "C60" -b BuildingBlocks-release  -r 1 

# Generate topology  files for C60
xyz2data.py --in_itp conj.itp --in_xyz mols/C60/acc1_C60_n1.xyz --out_data  mols/C60/acc1_C60_n1.data




