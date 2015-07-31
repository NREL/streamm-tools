#!/bin/bash
#
# Convert structure files into input files for molecular dynamics simulations 
#
#
# Check the correct path is set
#
./mol_md_setup.sh 
#
# Choose a structure generatated by mol_gen_ex3 (P3HT)
#
MOLDIR=mols/thiophene_R_hexane_/
MOLID=acc1_thiophene_R_hexane__n5
#
# Check to make sure the files were generated 
#
if [ ! -s  ${MOLDIR}${MOLID}.com ]; then 
    echo " "
    echo "---------------------------------------------------------------"
    echo " Example .com file not found:"
    echo "    "${MOLDIR}${MOLID}.com
    echo " rerun mol_gen_ex3.sh  "
    echo "---------------------------------------------------------------"
    echo " "
    exit
fi 
#
# Transplate the .com tile to .xyz file to view with your favorate viewer
#
python com2xyz.py --in_com ${MOLDIR}${MOLID}.com --out_xyz ${MOLDIR}${MOLID}.xyz
#
# Transplate the .com tile to .gro file 
#
python com2gro.py --in_com ${MOLDIR}${MOLID}.com --out_gro ${MOLDIR}${MOLID}.gro 
#
# Transplate the .com tile to .gro file, .itp file and .top file for use GROMACS
#
python com2gromacs.py --in_com ${MOLDIR}${MOLID}.com  --in_itp conj.itp --out_gro ${MOLDIR}${MOLID}.gro --out_itp ${MOLID}.itp --out_top ${MOLDIR}${MOLID}.top 
mv ${MOLID}.itp ${MOLDIR}
#
# Transplate the .com tile to .data file 
#
python com2data.py --in_com ${MOLDIR}${MOLID}.com  --in_itp conj.itp --out_data ${MOLDIR}${MOLID}.data
#
# Replicate the .data file to produce a larger simulations cell
#
# Replicate to random points with random rotations
python replicate_data.py  --in_data ${MOLDIR}${MOLID}.data --out_data ${MOLDIR}${MOLID}x10.data  --out_xyz t1.xyz  --mol_n 10
# Replicate on a grid
python replicate_data.py  --in_data ${MOLDIR}${MOLID}.data --out_data ${MOLDIR}${MOLID}x10_g.data  --out_xyz t1.xyz  --mol_n 10 --grid
#

