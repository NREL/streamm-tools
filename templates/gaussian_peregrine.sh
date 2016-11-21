#!/bin/bash

NPROC=<nproc>

scratch_dir=<scratch_dir>
cd $scratch_dir

SCRATCH=<scratch_dir>/GAUSSIAN-<tag>.$PBS_JOBID.src
SCRATCH2=/dev/shm

if [ -d $SCRATCH ]
then
   rm -rf $SCRATCH
fi
mkdir $SCRATCH


export GAUSS_SCRDIR=$SCRATCH2
ls -lah /dev/shm
rm -rf /dev/shm/*

echo "%NoSave" > header.text
echo "%RWF=$SCRATCH2/,1500MB,$SCRATCH/,-1" >> header.text

cp <input_com> body.text
cat header.text  body.text  >  <input_com> 

g09 < <input_com> >&  <tag>.log

formchk <tag>.chk
