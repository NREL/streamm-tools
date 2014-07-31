#!/usr/bin/env python

import os

set_path="PYTHONPATH=../src/:../AtomicPy/src/"
os.system(set_path)
run_py="python ../scripts/rdf.py  -v -j acc1_D1_R2R200_A2_R3_n1_R41n1R41n1R40n1__n5.json   --in_lammpsxyz  n5.xyz   --fftype_i 'C!'   --fftype_j 'C!'  -o test_rdf "

os.system(run_py)
