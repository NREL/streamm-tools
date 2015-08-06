#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build oligomethylthiophene
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

# Generate structure files for oligomethylthiophene

    echo " 
     ========================================================================
         Generating: 
              mols/thiophene_R_methane_/acc1_thiophene_R_methane__n5.xyz
     ========================================================================
     "
     donoracceptorsystems.py  "thiophene  ( R_methane  )" -b  BuildingBlocks-release -r 5 -p "180 0 "

# Generate topology  files for oligomethylthiophene
    echo " 
     ========================================================================
         Generating: 
              mols/thiophene_R_methane_/acc1_thiophene_R_methane__n5.data
     ========================================================================
     "
     xyz2data.py  --in_itp conj.itp --in_xyz  mols/thiophene_R_methane_/acc1_thiophene_R_methane__n5.xyz --out_data  mols/thiophene_R_methane_/acc1_thiophene_R_methane__n5.data




