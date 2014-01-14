Running torsion


Example: 

cp ${TOOL_PATH}/AtomicPy/torsion/*template ./

TOOL_PATH=/home/${USER}/streamm/tools
DAT=D100n2-5.dat 



# Make input files for torsional calculation based on input files listed in
# $DAT file 

python ${TOOL_PATH}/AtomicPy/torsion/mk_qm.py  --cluster  "peregrine" --submit --dih_temp acc1_dih.com.template   $DAT

# *dih.list can be modified for each dihedral angle to loop / fix /relax  and
# mk_qm can be reran. But the qm_dih lines in $DAT will need to be delete before
# mk_qm is reran. 


# Run calculations 
python ${TOOL_PATH}/AtomicPy/torsion/calc_qm.py  --cluster  "peregrine" --submit   $DAT

# Compile results into data files 
python ${TOOL_PATH}/AtomicPy/torsion/data_qm.py $DAT

# On local machine the data files can be retrieved using 
TOOL_PATH=/Users/${USER}/Scripts/streamm/tools
DAT=D100n2-5.dat


python ${TOOL_PATH}/AtomicPy/torsion/get_data.py  --cluster  "peregrine"   --info_dir /scratch/tkemper/opv_gen_torsion2  --info_file  $DAT
python ${TOOL_PATH}/AtomicPy/torsion/plot_data.py $DAT
python ${TOOL_PATH}/AtomicPy/torsion/vmd_data.py $DAT






 









