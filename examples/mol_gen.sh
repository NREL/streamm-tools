#!/bin/bash

# Examples of running the donoracceptorsystems.py script to build complex molecules 

# Files  and directories
# - donoracceptorsystems.py 
# Concatenates molecular cply files based on connectivity tags 
# 
# - BuildingBlocks
# 
#   Must contain the directories 
# 
#   -  acceptors   
#    -   donors            
#    - functional_groups spacers           
#    - terminals

source mol_gen_setup.sh 

echo "mol_gen_setup.sh finished with out errors "
echo "  Proceed to example 1 [y]yes [n]no"
read -e ans
# Example 1
if [ $ans == "y" ]; then
    source mol_gen_ex1.sh 
fi


echo "  Proceed to example 2 [y]yes [n]no"
read -e ans
# Example 2
if [ $ans == "y" ]; then
    source mol_gen_ex2.sh 
fi

echo "  Proceed to example 3 [y]yes [n]no"
read -e ans
# Example 2
if [ $ans == "y" ]; then
    source mol_gen_ex3.sh 
fi

echo "mol_gen finished"
exit
