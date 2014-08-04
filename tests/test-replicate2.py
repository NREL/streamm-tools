#!/usr/bin/env python

import os

set_path="PYTHONPATH=../src/:../AtomicPy/src/"
# os.system(set_path)
run_py="python ../scripts/replicate.py --gro SOL.gro --top SOL.top   --sol_gro SOL.gro   --sol_top SOL.top   --den_target 0.1  --atoms_target 200    --perc_sol 50 -v "
os.system(run_py)


