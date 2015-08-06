#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build oligothiophene
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

# Generate structure files for oligothiophene

donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release  -r 5

donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release    -r 5 -p "0"

donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release   -r 5  -p "180 0 "

# Generate topology  files for oligothiophene
xyz2data.py --in_itp conj.itp --in_xyz mols/thiophene/acc1_thiophene_n5.xyz --out_data  mols/thiophene/acc1_thiophene_n5.data




