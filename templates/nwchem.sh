#!/bin/bash

scratch_dir=<scratch>
cd $scratch_dir

NPROC=<nproc>
INPUT_FILE=<input_nw> 
mpirun -np $NPROC  nwchem   $INPUT_FILE    >& <tag>.log 
