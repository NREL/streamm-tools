#!/bin/bash

# Examples of replicating P3HT, C50 and benzene in molecular dynamics simulation box 
./examples_setup.sh 

exitcode=$?
if [ $exitcode == 0  ]
then
    
    echo " 
     ========================================================================
         Generating: 
              thiophene.cply
              R_hexane.cply 

     with new cply format 
     ========================================================================
     "
    python ../functions/cply2cply.py  -l -b -r -f   --in_cply BuildingBlocks/donors/thiophene.cply  -o thiophene 
    python ../functions/cply2cply.py  -l -b -r -f   --in_cply BuildingBlocks/functional_groups/R_hexane.cply -o R_hexane 
    sed "s/rg/term/g" R_hexane.cply > R_hexane_v2.cply
    
    # Edit charges now 

    echo " 
     ========================================================================
         Generating: 
              thiophene_R_hexane_v2_n5
     ========================================================================
     "
    
    python ../functions/catcply.py  "thiophene  ( R_hexane_v2  )"  -r 5 

    python ../functions/cply2cply.py  -l -b -r -f   --in_cply BuildingBlocks/fullerene/C60.cply  -o C60
    
    python ../functions/cply2film.py  --ranmol   "thiophene_R_hexane_v2_n5  10  C60 10 "  -v  -o P3HT_x10_C60_x10 
    
    python ../functions/cply2data.py  --in_cply P3HT_x10_C60_x10.cply  --in_itp biaryl_v1.itp  --out_data P3HT_x10_C60_x10.data 


fi

exit 0 