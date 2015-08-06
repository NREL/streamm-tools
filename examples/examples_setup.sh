#!/bin/bash


# Check dependencies 
echo "================================================================================"
echo "                 Checking  PYTHONPATH and dependent repos                    "
echo " "

command -v donoracceptorsystems.py >/dev/null 2>&1 || { echo >&2 "I require donoracceptorsystems.py but it's not installed.  Aborting."; exit 1; }

pull_BuildingBlocks.py


exe_file="replicate_data.py      replicate_multimol.py  xyz2data.py            xyz2gromacs.py"
for EXE_FILE in $exe_file
do
    if [ -s $EXE_FILE ]; then 
	echo " Executable file " $EXE_FILE " has been found "
	chmod 755 $EXE_FILE 
    else
	echo " Executable file " $EXE_FILE " has not been found "
	exit 1
    fi

done


template_files="donoracceptor.com.template  donoracceptor.pbs.template   donoracceptor.pbs.template   donoracceptor.com.template.r1    donoracceptor.com.template.r2"

for TEMPLATE_FILE in $template_files
do
    if [ -s $TEMPLATE_FILE ]; then 
	echo " Template file " $TEMPLATE_FILE " has been found "
    else
	echo " Template file " $TEMPLATE_FILE " has not been found "
	exit 1
    fi

done

echo "mol_gen_setup.sh finished with out errors "
