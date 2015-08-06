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

    if [ -s ${TOOLS_PATH}/da_builder/donoracceptorsystems.py ]; then
	echo " "
	echo "---------------------------------------------------------------"
	echo " donoracceptorsystems.py found in:" 
	echo "      " ${TOOLS_PATH}/da_builder/
	echo "---------------------------------------------------------------"
	echo " "
    else
	echo " "
	echo "---------------------------------------------------------------"
	echo " donoracceptorsystems.py not found  in:" 
	echo "      " ${TOOLS_PATH}/da_builder/
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

    if [ -s ${TOOLS_PATH}/da_builder/donoracceptorsystems.py ]; then
	echo " "
	echo "---------------------------------------------------------------"
	echo " donoracceptorsystems.py found in:" 
	echo "      " ${TOOLS_PATH}/da_builder/
	echo "---------------------------------------------------------------"
	echo " "
    else
	echo " "
	echo "---------------------------------------------------------------"
	echo " donoracceptorsystems.py not found  in:" 
	echo "      " ${TOOLS_PATH}/da_builder/
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

#
# Check if the BuildingBlocks repo has been pulled
#
if [ -z  ${BB_PATH} ]; then
    BB_PATH=${TOOLS_PATH}/../BuildingBlocks-release  
    export BB_PATH
    echo " "
    echo "---------------------------------------------------------------"
    echo " BuildingBlocks repo path being set to:"
    echo "    " ${BB_PATH}
    echo "---------------------------------------------------------------"
    echo " "
else
    echo " "
    echo "---------------------------------------------------------------"
    echo " Using existing BuildingBlocks repo path  BB_PATH:"
    echo "    " ${BB_PATH}
    echo "---------------------------------------------------------------"
    echo " "

fi


if [ -d  ${BB_PATH}/donors -a  -d  ${BB_PATH}/acceptors -a -d  ${BB_PATH}/functional_groups -a -d  ${BB_PATH}/spacers -a -d  ${BB_PATH}/terminals   ]; then
    echo " "
    echo "---------------------------------------------------------------"
    echo " BuildingBlocks repo has been found in:"
    echo "    " ${BB_PATH}
    echo "---------------------------------------------------------------"
    echo " "
else
    echo " "
    echo "---------------------------------------------------------------"
    echo " BuildingBlocks repo has not been found in:"
    echo "    " ${BB_PATH}
    echo " or is missing a subdirectory "
    echo "    donors            spacers   acceptors         functional_groups terminals"
    echo " rerun setup.sh in the opv-project directory to pull this repo "
    echo " or set BB_PATH to it's location "
    echo " eg: BB_PATH='path-to-BuildingBlocks'/BuildingBlocks"
    echo "     then "
    echo "     export BB_PATH"
    echo "---------------------------------------------------------------"
    
fi



#
# Check if the templates are in the opv-project directory 
#
#if [ -z  ${TEMPLATE_PATH} ]; then
#    TEMPLATE_PATH=${TOOLS_PATH}/../
#    export TEMPLATE_PATH
#    echo " "
#    echo "---------------------------------------------------------------"
#    echo " Template path being set to:"
#    echo "    " ${TEMPLATE_PATH}
#    echo "---------------------------------------------------------------"
#    echo " "
#fi


#if [ -s ${TEMPLATE_PATH}/donoracceptor.com.template -a -s ${TEMPLATE_PATH}/donoracceptor.pbs.template  -a -s ${TEMPLATE_PATH}/donoracceptor.pbs.template  -a -s ${TEMPLATE_PATH}/donoracceptor.com.template.r1   -a -s ${TEMPLATE_PATH}/donoracceptor.com.template.r2 ]; then
#    echo " "
#    echo "---------------------------------------------------------------"
#    echo " TEMPLATEs have been found in:"
#    echo "    " ${TEMPLATE_PATH}
#    echo "---------------------------------------------------------------"
#    echo " "
#else
#    echo " "
#    echo "---------------------------------------------------------------"
#    echo " TEMPLATEs have not been found in:"
#    echo "    " ${TEMPLATE_PATH}
#    echo " rerun setup.sh in the opv-project directory to pull these templates  "
#    echo " or set TEMPLATE_PATH to it's location "
#    echo " eg: TEMPLATE_PATH='path-to-templates'/"
#    echo "     then "
#    echo "     export TEMPLATE_PATH"
#    echo "---------------------------------------------------------------"
#fi

# To create a gaussian input file you can add the gaussian templates to the current directory
#cp ${TEMPLATE_PATH}/donoracceptor.com.template ./
#cp ${TEMPLATE_PATH}/donoracceptor.com.template.r1 ./
#cp ${TEMPLATE_PATH}/donoracceptor.com.template.r2 ./
#cp ${TEMPLATE_PATH}/donoracceptor.pbs.template ./

echo "mol_gen_setup.sh finished with out errors "
