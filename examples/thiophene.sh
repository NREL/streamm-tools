#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build thiophene
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

exitcode=$?
if [ $exitcode == 0  ]
then
    # Generate structure files for thiophene
    echo " 
     ========================================================================
         Generating: 
              mols/thiophene/acc1_thiophene_n1.xyz
     ========================================================================
     "
     donoracceptorsystems.py  "thiophene" -b BuildingBlocks  -r 1 

# Generate topology  files for thiophene
    echo " 
     ========================================================================
         Generating: 
              mols/thiophene/acc1_thiophene_n1.xyz
     ========================================================================
     "
     xyz2data.py --in_itp conj.itp --in_xyz mols/thiophene/acc1_thiophene_n1.xyz --out_data  mols/thiophene/acc1_thiophene_n1.data
     echo " 
     ========================================================================
         Generating: 
              mols/thiophene/acc1_thiophene_n1.gro
              mols/thiophene/acc1_thiophene_n1.top
              mols/thiophene/acc1_thiophene_n1.itp
     ========================================================================
     "
     xyz2gromacs.py --in_itp conj.itp --in_xyz  mols/thiophene/acc1_thiophene_n1.xyz --out_gro mols/thiophene/acc1_thiophene_n1.gro   --out_top mols/thiophene/acc1_thiophene_n1.top   --out_itp  acc1_thiophene_n1.itp
     mv acc1_thiophene_n1.itp  mols/thiophene/
    echo " 
     ========================================================================
         Generateration of input files for thiophene finished  
     ========================================================================
     "
fi 

exit 
