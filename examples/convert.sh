#!/bin/bash

# Convert .com gaussain files into file types for other software 
MOLDIR=mols/thiophene_R_methane_/
MOLID=acc1_thiophene_R_methane__n1

python com2xyz.py --in_com ${MOLDIR}${MOLID}.com --out_xyz ${MOLDIR}${MOLID}.xyz
python com2gro.py --in_com ${MOLDIR}${MOLID}.com --out_gro ${MOLDIR}${MOLID}.gro 
python com2gromacs.py --in_com ${MOLDIR}${MOLID}.com  --out_gro ${MOLDIR}${MOLID}.gro --in_itp conj.itp --out_itp ${MOLID}.itp --out_top ${MOLDIR}${MOLID}.top 
mv ${MOLID}.itp ${MOLDIR}

python com2data.py --in_com ${MOLDIR}${MOLID}.com  --in_itp conj.itp --out_data ${MOLDIR}${MOLID}.data

# python replicate_data.py  --in_data ${MOLDIR}${MOLID}.data --out_data 


