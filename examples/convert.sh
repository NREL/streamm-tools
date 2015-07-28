#!/bin/bash

# Convert .com gaussain files into file types for other software 

python com2xyz.py --in_com mols/thiophene_R_methane_/acc1_thiophene_R_methane__n1.com --out_xyz t.xyz
python com2gro.py --in_com mols/thiophene_R_methane_/acc1_thiophene_R_methane__n1.com --out_gro t.gro 
