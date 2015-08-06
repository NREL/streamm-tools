#!/bin/bash

# Examples of replicating P3HT, C50 and benzene in molecular dynamics simulation box 
./examples_setup.sh 

./P3HT.sh 
./C60.sh 
./benzene.sh 

    echo " 
     ========================================================================
         Generating: 
              mols/P3HTC60benzene_mix1.data
              mols/P3HTC60benzene_mix1.xyz
     ========================================================================
     "
    ./replicate_multimol.py  --mol1_data mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data  --mol1_n 10 --mol2_data mols/C60/acc1_C60_n1.data --mol2_n 10 --mol3_data mols/benzene/acc1_benzene_n1.data   --out_data  mols/P3HTC60benzene_mix1.data    --out_xyz mols/P3HTC60benzene_mix1.xyz  --mol3_n 50 -v  
