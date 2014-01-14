Running torsion


Example: 

cp ${TOOL_PATH}/AtomicPy/torsion/*template ./

TOOL_PATH=/home/${USER}/streamm/tools
DAT=D100


# Make input files for torsional calculation based on input files listed in
# $DAT file 

python ${TOOL_PATH}/AtomicPy/torsion/mk_qm.py  --cluster  "peregrine" --submit --dih_temp acc1_dih.com.template   $DAT

# *dih.list can be modified for each dihedral angle to loop / fix /relax  and
# mk_qm can be reran. But the qm_dih lines in $DAT will need to be delete before
# mk_qm is reran. 

python ${TOOL_PATH}/AtomicPy/torsion/calc_qm.py  --cluster  "peregrine" --submit   $DAT
python ${TOOL_PATH}/AtomicPy/torsion/data_qm.py $DAT









