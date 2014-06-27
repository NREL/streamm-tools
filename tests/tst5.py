#!/usr/bin/env python

import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from structureContainer import StructureContainer

try:
    import numpy as np
except:
    print "Error: numpy module not found, check PYTHONPATH for numpy install or module version"
    sys.exit(3)

try:
    import matplotlib.pyplot as plt
except:
    print "Error: matplotlib.pyplot module not found, check PYTHONPATH or install from source"
    sys.exit(3)

try:
    from mpl_toolkits.mplot3d import Axes3D
except:
    print "Error: mpl_toolkits.mplot3d module not found, check PYTHONPATH or install from source"
    sys.exit(3)


try:
    from runjobs import Geometry, IO
except:
    print "runjobs module not found"
    print "Try adding the following to your .bashrc file"
    print "   PYTHONPATH='path-to-runjobs.py':$PYTHONPATH'"
    print "   export PYTHONPATH"
    sys.exit(3)


try:
    import mpiNREL
except:
    print "mpiNREL Module not found"
    print "Try adding the following to your .bashrc file"
    print "   PYTHONPATH='path-to-runjobs.py':$PYTHONPATH'"
    print "   export PYTHONPATH"
    sys.exit(3)



# Get useful methods
g=Geometry()

# Get comm object
p = mpiNREL.getMPIObject()
rank = p.getRank()
size = p.getCommSize()


def setRandomGridPoints(nptsL, rshift):
    """
    Creates list of elements containing:
    -->   [ [global-index, x-pos, y-pos, z-pos], [ ] .... ]
    On a cubic grid with random displacements
    Index begins at '1'
    Box lengths are assumed to be from particles with ptcl_dist in between.
    End particles are next to box 'edges'
    
    lpts   -- number of points on each cubic side (list)
    rshift -- spatial distance to random shift a position
    """

    global g
    global rad_avg
    global ptcl_dist

    ptclCdist = (2.0*rad_avg) + ptcl_dist

    nx=nptsL[0]
    ny=nptsL[1]
    nz=nptsL[2]
    ptsList = []

    n = 1
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):

                if (x == 0) or (x == nx-1):
                    xNew = (x*ptclCdist)
                else:
                    xNew = (x*ptclCdist) + random.uniform(-rshift,rshift)

                # if (x == 0) or (x == nx-1):
                #    xNew = (x*ptclCdist)
                # else:
                #    xNew = (x*ptclCdist) + random.uniform(-rshift,rshift)

                xNew = ( (x*ptclCdist) + rad_avg) + random.uniform(-rshift,rshift)
                yNew = ( (y*ptclCdist) + rad_avg) + random.uniform(-rshift,rshift)
                zNew = ( (z*ptclCdist) + rad_avg) + random.uniform(-rshift,rshift)

                # xNew = (x*ptclCdist) + random.gauss(0.0,rshift)
                # yNew = (y*ptclCdist) + random.gauss(0.0,rshift)
                # zNew = (z*ptclCdist) + random.gauss(0.0,rshift)

                randomPt = [n, xNew, yNew, zNew]
                
                n = n+1
                ptsList.append(randomPt)

    return ptsList




def setSystemSizes(numQDs, spacing):

    Lx = (spacing*(numQDs[0]-1)) + spacing
    Ly = (spacing*(numQDs[1]-1)) + spacing
    Lz = (spacing*(numQDs[2]-1)) + spacing

    sizes = [Lx, Ly, Lz]
    return sizes




#############################################################################################
#
# Generate list of  [gid, x,y,z] points ---> allPoints list
# copied on all processors. # Global pt index included
# UNITS distance --> A [0.1 nm]
#
#############################################################################################

rad_avg = 18.5
rad_sig = 0.005
ptcl_dist = 10.0

# NOTE: this is an arbitrary choice
# Particles shift uniformly between -rshift,rshift
rshift = ptcl_dist/5.0
rshift = 0.0

cutoff_dist   = (2.00 * rad_avg) + ptcl_dist + rshift
cutoff_dist   = math.sqrt(2)*((2.00 * rad_avg) + ptcl_dist)
nQD_spacing   = (2.00 * rad_avg) + ptcl_dist
cutoffDistSqr = pow(cutoff_dist, 2.0)

# Makes number of points along side grid
nQDs = [20, 20, 20]
nQDs = [5, 5, 1]
nQDs = [5, 5, 5]
sysLs = setSystemSizes(nQDs, nQD_spacing)
# bcLs = sysLs
bcLs = [0.0, sysLs[1], sysLs[2]]
bcLs = [0.0, 0.0, 0.0]

if rank == 0:
    print "rad_avg = ", rad_avg
    print "rad_sig = ", rad_sig
    print "rshift  = ", rshift
    print " "
    print "cutoff_dist = ", cutoff_dist
    print " "
    print "bcLs  = ", bcLs
    print "sysLs = ", sysLs
    print " "


# Sync random stream for tests
random.seed(0)
p.barrier()

tmpPoints = []
tmpPoints = setRandomGridPoints(nQDs, rshift) # Generate global points list
allPoints = p.bcast(tmpPoints)                # To ensure random points are same across procs

numpoints = len(allPoints)                    # Total number of points (for diag)
if numpoints == 0:
    print "Number of points generated == 0"
    sys.exit(3)

myPoints  = p.splitListOnProcs(allPoints) # Split elements on across processors
#############################################################################################


nanoPtcls = ParticleContainer()

for ptPos in allPoints:

    if ptPos[0] & 1:
        pt = Particle(ptPos[1:], type="QDsmall", mass=1.0)
    else:
        pt = Particle(ptPos[1:], type="QDbig", mass=2.0)

    tagsD = {"molnum":1}
    pt.setTagsDict(tagsD)
    nanoPtcls.put(pt)



strucQD = StructureContainer(nanoPtcls)

boxSizes = [ [0.0, sysLs[0] ], [0.0, sysLs[1] ], [0.0, sysLs[2] ] ]
strucQD.setBoxLengths(boxSizes)

paramMap = {("QDsmall", "epsilon"):1.20, ("QDsmall", "sigma"):10.10,
            ("QDbig",   "epsilon"):1.11, ("QDbig",   "sigma"):12.38 }

if rank == 0:
    nanoPtcls.scatterPlot()
    strucQD.dumpLammpsInputFile("qd_lmp_in", paramMap)
