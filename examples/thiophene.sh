#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build thiophene
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

# Generate structure files for thiophene
donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release  -r 1 

# Generate topology  files for thiophene
xyz2data.py --in_itp conj.itp --in_xyz mols/thiophene/acc1_thiophene_n1.xyz --out_data  mols/thiophene/acc1_thiophene_n1.data




