#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build oligothiophene
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

exitcode=$?
if [ $exitcode == 0  ]
then
    # Generate structure files for oligothiophene
    echo " 
     ========================================================================
         Generating: 
              mols/thiophene/acc1_thiophene_n5.xyz
     ========================================================================
     "
     donoracceptorsystems.py  "thiophene" -b BuildingBlocks -r 5

     donoracceptorsystems.py  "thiophene" -b BuildingBlocks -r 5 -p "0"

     donoracceptorsystems.py  "thiophene" -b BuildingBlocks -r 5  -p "180 0 "

# Generate topology  files for oligothiophene
    echo " 
     ========================================================================
         Generating: 
              mols/thiophene/acc1_thiophene_n5.data
     ========================================================================
     "
     xyz2data.py --in_itp conj.itp --in_xyz mols/thiophene/acc1_thiophene_n5.xyz --out_data  mols/thiophene/acc1_thiophene_n5.data
     echo " 
     ========================================================================
         Generating: 
              mols/thiophene/acc1_thiophene_n5.gro
              mols/thiophene/acc1_thiophene_n5.top
              mols/thiophene/acc1_thiophene_n5.itp
     ========================================================================
     "
     xyz2gromacs.py --in_itp conj.itp --in_xyz  mols/thiophene/acc1_thiophene_n5.xyz --out_gro mols/thiophene/acc1_thiophene_n5.gro   --out_top mols/thiophene/acc1_thiophene_n5.top   --out_itp  acc1_thiophene_n5.itp
     mv acc1_thiophene_n5.itp  mols/thiophene/
    echo " 
     ========================================================================
         Generateration of input files for oligo-thiophene finished  
     ========================================================================
     "
fi

exit 



