#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build benzene
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

# Generate structure files for benzene
    echo " 
     ========================================================================
         Generating: 
              mols/benzene/acc1_benzene_n1.xyz
     ========================================================================
     "
     donoracceptorsystems.py  "benzene" -b BuildingBlocks-release  -r 1 

# Generate topology  files for benzene
    echo " 
     ========================================================================
         Generating: 
              mols/benzene/acc1_benzene_n1.xyz
     ========================================================================
     "
     xyz2data.py --in_itp conj.itp --in_xyz mols/benzene/acc1_benzene_n1.xyz --out_data  mols/benzene/acc1_benzene_n1.data




