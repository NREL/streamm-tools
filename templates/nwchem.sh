#!/bin/bash

scratch_dir=<scratch_dir>
cd $scratch_dir

NPROC=<n_proc>
INPUT_FILE=<input_nw> 
mpirun -np $NPROC  nwchem   $INPUT_FILE    >& <output_log>
