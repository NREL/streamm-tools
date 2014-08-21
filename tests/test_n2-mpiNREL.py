#!/usr/bin/env python

import mpiNREL
p = mpiNREL.getMPIObject(False, localVerbose=False)
# p = mpiNREL.getMPIObject(True, localVerbose=True)

# MPI setup
rank = p.getRank()
size = p.getCommSize()

x = 0              # Sets x on all processors
p.barrier()
if rank == 0:      # Resets x on proc 0
    x = 12.3456
p.barrier()

# Check x setting on each processor
for proc in range(size):
    if proc == rank:
        print "Before broadcast x = ", x, " on processor ", proc
    p.barrier()
p.barrier()

# Broadcast x to all processors (default root proc is '0')
xAll = p.bcast(x)
p.barrier()

# Check x setting on each processor
for proc in range(size):
    if proc == rank:
        print "After broadcast x = ", xAll, " on processor ", proc
    p.barrier()
p.barrier()


