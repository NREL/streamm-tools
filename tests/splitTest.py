#!/usr/bin/env python

import sys
import mpiNREL
p = mpiNREL.getMPIObject(False, localVerbose=False)
p.barrier()

# MPI setup
rank = p.getRank()
size = p.getCommSize()
p.barrier()

allPoints = range(10)
p.barrier()

myPoints = p.splitListOnProcs(allPoints) # Split elements on across processors
p.barrier()
    
if rank == 0:
    print "------------------------------------------------------------------"
    print "For processors = ", str(size) + "\n"    
p.barrier()

for iproc in range(size):
    if rank == iproc:
        print "myPoints = ", myPoints, " on proc ", rank
    p.barrier()

p.barrier()
