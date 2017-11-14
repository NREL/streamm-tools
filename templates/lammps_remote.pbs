#!/bin/bash
#PBS -A <allocation>
#PBS -l walltime=<walltime>:00:00 # WALLTIME limit
#PBS -l nodes=<nodes>:ppn=<ppn> # Number of nodes, put 16 processes on each
#PBS -l pmem=1500
#PBS -M <e-mail>
#PBS -m e
#PBS -N <tag>
#PBS -q <queue>
#PBS -l feature=<feature>

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


