#!/usr/bin/env python

import sys
import mpiNREL
p = mpiNREL.getMPIObject(False, localVerbose=False)

# MPI setup
rank = p.getRank()
size = p.getCommSize()

allPoints = range(10)
p.barrier()

myPoints = p.splitListOnProcs(allPoints) # Split elements on across processors

for iproc in range(size):
    if rank == iproc:
        print "myPoints = ", myPoints, " on proc ", rank
    p.barrier()
