#!/bin/bash


# Check dependencies 
echo "================================================================================"
#echo "                 Checking PYTHONPATH, PATH and TOOLS_PATH                        "
echo "                 Checking  OPV-PROJECT_PATH                        "
echo " "


#
# Check if TOOLS_PATH is set
#
if [ -n "${TOOLS_PATH}" ]; then
    echo " "
    echo "---------------------------------------------------------------"
    echo " TOOLS_PATH set to  " ${TOOLS_PATH}
    echo "---------------------------------------------------------------"
    echo " "

    if [ -s ${TOOLS_PATH}/src/structureContainer.py ]; then
	echo " "
	echo "---------------------------------------------------------------"
	echo " structureContainer.py found in:" 
	echo "      " ${TOOLS_PATH}/da_src/
	echo "---------------------------------------------------------------"
	echo " "
    else
	echo " "
	echo "---------------------------------------------------------------"
	echo " structureContainer.py not found  in:" 
	echo "      " ${TOOLS_PATH}/src/
	echo " rerun setup.sh in the opv-project directory "
	echo " or set TOOLS_PATH to location of the tools repo "
	echo " "
	echo " eg: TOOLS_PATH='path-to-tools'/tools"
	echo "     then "
	echo "     export TOOLS_PATH"
	echo "---------------------------------------------------------------"
	echo " "
	exit
    fi
else
    echo " "
    echo "---------------------------------------------------------------"
    echo "TOOLS_PATH is unset. Setting TOOLS_PATH to " '../' 
    echo "---------------------------------------------------------------"
    echo " "

    if [ -s ${TOOLS_PATH}/src/structureContainer.py ]; then
	echo " "
	echo "---------------------------------------------------------------"
	echo " structureContainer.py found in:" 
	echo "      " ${TOOLS_PATH}/src/
	echo "---------------------------------------------------------------"
	echo " "
    else
	echo " "
	echo "---------------------------------------------------------------"
	echo " structureContainer.py not found  in:" 
	echo "      " ${TOOLS_PATH}/src/
	echo " rerun setup.sh in the opv-project directory "
	echo " or set TOOLS_PATH to location of the tools repo "
	echo " "
	echo " eg: TOOLS_PATH='path-to-tools'/tools"
	echo "     then "
	echo "     export TOOLS_PATH"
	echo "---------------------------------------------------------------"
	echo " "
	exit
    fi
fi


echo "mol_md_setup.sh finished with out errors "
