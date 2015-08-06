#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build P3HT
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

# Generate structure files for P3HT

donoracceptorsystems.py  "thiophene  ( R_hexane  )" -b  BuildingBlocks-release -r 5 -p "140 40"

# Generate topology  files for P3HT
xyz2data.py  --in_itp conj.itp     --in_xyz  mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.xyz     --out_data  mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data





