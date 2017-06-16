#!/bin/bash
echo " "

TOP_PATH=`pwd`
echo "Setting top path $TOP_PATH in setup"

#
# Setting paths for python tools
#
echo "================================================================================"
echo "                 Setting PYTHONPATH, PATH                                       "
echo " "

#
# Cleaning environment of previous STREAMM related paths
#
$TOP_PATH/src/cleanpath.py
echo " "
source clean-PYTHONPATH.sh
source clean-PATH.sh
echo ""
echo "Removing clean-*.sh scripts"
rm -rf clean-*PATH.sh

#
# Set PYTHONPATH, PATH  and TOOLS_PATH
#
PYTHONPATH=$TOP_PATH/src:$PYTHONPATH
PYTHONPATH=./:$PYTHONPATH
export PYTHONPATH

PATH=$TOP_PATH/src:$PATH
PATH=$TOP_PATH/da_builder:$PATH
PATH=./:$PATH
export PATH

echo " "
echo "Setting PYTHONPATH = $PYTHONPATH"
echo " "
echo "Setting PATH       = $PATH"
echo " "
echo "NOTE: be sure to have 'source'-ed this file to include paths in environment"
echo " "
echo "================================================================================"
