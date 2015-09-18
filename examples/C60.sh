#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build C60
./examples_setup.sh 

# The backbone of the molecule can be composed of  acceptors, spacers or
# donors, and can be  decorated with molecules from the
# functional_groups directory  and terminated with molecules from the
# terminals directory. Run 

exitcode=$?
if [ $exitcode == 0  ]
then
    # Generate structure files for C60
    echo " 
     ========================================================================
         Generating: 
              mols/C60/acc1_C60_n1.xyz
     ========================================================================
     "
     donoracceptorsystems.py  "C60" -b BuildingBlocks  -r 1 

# Generate topology  files for C60
    echo " 
     ========================================================================
         Generating: 
              mols/C60/acc1_C60_n1.data
     ========================================================================
     "
     xyz2data.py --in_itp conj.itp --in_xyz mols/C60/acc1_C60_n1.xyz --out_data  mols/C60/acc1_C60_n1.data
     echo " 
     ========================================================================
         Generating: 
              mols/C60/acc1_C60_n1.gro
              mols/C60/acc1_C60_n1.top
              mols/C60/acc1_C60_n1.itp
     ========================================================================
     "
     xyz2gromacs.py --in_itp conj.itp --in_xyz  mols/C60/acc1_C60_n1.xyz --out_gro mols/C60/acc1_C60_n1.gro   --out_top mols/C60/acc1_C60_n1.top   --out_itp  acc1_C60_n1.itp
     mv acc1_C60_n1.itp  mols/C60/
    echo " 
     ========================================================================
         Generateration of input files for C60 finished  
     ========================================================================
     "
fi

exit 
