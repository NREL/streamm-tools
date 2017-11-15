#!/bin/bash

scratch_dir=<scratch>
cd $scratch_dir

module purge;
module load python/2.7.8
module load mdanalysis/0.8.1
module load comp-intel/13.1.3      impi-intel/4.1.1-13.1.3
module load lammps/14Aug2013

LAMMPSEXE="lmp"
NPROC=<nproc>
echo " Running  $LAMMPSEXE   on  $NPROC processors "

mpirun -np $NPROC $LAMMPSEXE -in <input_in>  -log <tag>.log  >> <tag>.out

