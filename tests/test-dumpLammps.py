#!/usr/bin/env python

import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from bonds import Bond
from bonds import BondContainer

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




def setRandomGridPoints(nptsL, rad_avg, ptcl_dist, rshift=0.0):
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

    #global rad_avg
    #global ptcl_dist

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

                xNew = ( (x*ptclCdist) + rad_avg) + random.uniform(-rshift,rshift)
                yNew = ( (y*ptclCdist) + rad_avg) + random.uniform(-rshift,rshift)
                zNew = ( (z*ptclCdist) + rad_avg) + random.uniform(-rshift,rshift)

                randomPt = [n, xNew, yNew, zNew]
                
                n = n+1
                ptsList.append(randomPt)

    return ptsList


def inBetween(val, bounds):
    """
    Test if val in between values in 'bounds'
    Not including bounds
    """

    if (val < bounds[1]) and (val > bounds[0]):
        return True
    else:
        return False


def setBC1InXFromPos(allPts, rshift, sysLs, rad_avg, ptcl_dist, rank):
    """
    Sets the boundary condition flag (BC) [-1,0,1] based
    on shift used to lay down particles. Assumes
    
    Args:
        allPts -- list of [index, x, y, z]
        rshift -- Gaussian random shift used to displace particles
    """

    xL = sysLs[0]  # _min/_max
    yL = sysLs[1]
    zL = sysLs[2]
    
    # Used to calculate (estimate) particles centers at left/right electrode
    leftBuffer  = rad_avg
    rightBuffer = rad_avg + ptcl_dist

    if (rshift == 0.0):
        winSize = 0.01
    else:
        winSize = 3.5*rshift

    # Range window at 'left  x' side
    # xRangeL = [leftBuffer-winSize, leftBuffer+winSize]
    xRangeL = [xL[0] + leftBuffer - winSize, xL[0] + leftBuffer + winSize]

    # Range window at 'right x' side
    # xRangeR = [sysLs[0] - rightBuffer - winSize, sysLs[0] - rightBuffer + winSize]
    xRangeR = [xL[1] - rightBuffer - winSize, xL[1] - rightBuffer + winSize]
        
    electrodeSizeL = 0
    electrodeSizeR = 0

    bcFlags = []

    for pt in allPts:

        gid  = pt[0]
        xpos = pt[1]

        if inBetween(xpos, xRangeL):
            bcFlags.append([gid, 1.0])
            electrodeSizeL += 1

        elif inBetween(xpos, xRangeR):
            bcFlags.append([gid, -1.0])
            electrodeSizeR += 1

        else:
            bcFlags.append([gid, 0.0])

    if rank == 0:
        print "xRangeL = ", xRangeL
        print "xRangeR = ", xRangeR
        print "Left  electrode contains ", electrodeSizeL, " nodes"
        print "Right electrode contains ", electrodeSizeR, " nodes"

    return bcFlags


def setSystemSizes(numQDs, spacing):

    Lx = (spacing*(numQDs[0]-1)) + spacing
    Ly = (spacing*(numQDs[1]-1)) + spacing
    Lz = (spacing*(numQDs[2]-1)) + spacing

    sizes = [[0.0, Lx], [0.0, Ly], [0.0, Lz]]
    return sizes


def isPtInCylinder(pt, axisPt, axisDir, radius=0.0):
    """
    Deterimine whether pt is inside a cylinder aligned with box
    axis. NOTE: must have Geometry object

    Args:
        pt      (list): point positions [x,y,z]
        axisPt  (list): location of 'projection' of cylinder axis [x,y], or [x,z]....
        axisDir (str):  axis alignment (eg "xy" "yz" "xz")
        radius  (float): size of cylinder
        
    Returns: true if inside cylinder
    """
    # Geometry methods
    global g

    # Select direction
    if axisDir == "xy":
        axisColumns = [0,1]
    elif axisDir == "yz":
        axisColumns = [1,2]
    elif axisDir == "xz":
        axisColumns = [0,2]
    else:
        print "isPtInCylinder: axisDir not recognized"
        sys.exit(3)

    # Pick out 'columns' specified in axisDir from input pt
    ptProj = [ pt[i] for i in axisColumns]

    dist = g.rdist(ptProj, axisPt)

    if (dist < radius):
        return True
    else:
        return False


def isPtInLayer(pt, layerPos, layerDir, width=0.0):
    """
    Deterimine whether pt is inside a layer, alternates for multiple layers
    NOTE: must have Geometry object

    Args:
        pt       (list): point positions [x,y,z]
        layerPos (list): pos of layer center (dir specified by layerDir)
        layerDir  (str): plane layer is in (eg "xy" "yz" "xz")
        width   (float): width of layer
        
    Returns:
           true if inside layer
    """
    
    # Geometry methods
    global g

    # Select direction
    if layerDir == "xy":
        layerColumns = 2
    elif layerDir == "yz":
        layerColumn = 0
    elif layerDir == "xz":
        layerColumn = 1
    else:
        print "isPtInLayer: axisDir not recognized"
        sys.exit(3)

    # Pick out coord dir specified in axisDir from input pt
    ptProj = pt[layerColumn]

    # Check if point in between edges of layers
    layerEdges = [layerPos-(width/2.0),layerPos+(width/2.0)]
    if inBetween(ptProj, layerEdges):
        return True
    else:
        return False


def setQDNeighbors(localList, globalList, cutoff_dist, bcLs):
    """
    Calculate N^2 distances on multiple processors and parse out neighbor graph info

    Args:
        localList   (list):  [ [x1, y1, z1], [x2, y2, z2], ... ] local to processor
        globalList  (list):  ' ' all global points for calculating dist-s to input local pts
        cutoff_dist (float): Distance at which to consider points neightbors
        bcLs        (list):  List of boundary condition flags
        
    Returns:
        allDist, allNeighbors, localTotalNeighs
    """

    # Local geometry object
    g = Geometry()

    allDist      = []
    allNeighbors = []
    allBonds = []
    localTotalNeighs = 0
    cutoffDistSqr = pow(cutoff_dist, 2.0)
    
    if rank == 0:
        print "Number of points on proc-0", len(localList)
    p.barrier()

    start_time = time.time()
    for n,pt_i in enumerate(localList): # Loop on subset of points per processor

        if rank == 0:
            if n % 20 == 0:
                percentDone = 100.0*float(n+1)/float(len(localList))
                print "Proc 0 ", percentDone, "-% done"

        gid_i = pt_i[0]                # Parse out global-ID from element
        pos_i = pt_i[1:]               # Parse out x,y,z from element
        inghList = []                  # Temp list for holder neighbors gids for pt_i

        for pt_j in globalList:         # Loop on entire set of points (owned on each proc)

            pos_j = pt_j[1:]                         # Parse out x,y,z from element

            # dist = g.rdist(pos_i, pos_j)           # Calculate raw distance between pt1--pt2
            # dist = g.rdistPBCbox(pos_i, pos_j, bcLs) # Calc distance between pt1--pt2 in PBC box

            # Quick bounds check to speed, NOTE: only for non-PBC (check)
            inBounds = g.rdistInBounds(pos_i, pos_j, cutoff_dist)
            if (not inBounds):
                continue

            distSqr = g.rdistPBCboxSqr(pos_i, pos_j, bcLs) # Calc distance between pt1--pt2 in PBC box
            gid_j = pt_j[0]                                # Parse out global-ID from element

            if ( (distSqr < cutoffDistSqr) and (gid_i != gid_j) ): # Distance check (exclude self)

                 # if not nanoBonds.hasBond([gid_i, gid_j]):     # Check if bond is unique
                 #   bnd = Bond(gid_i, gid_j, type="normBond") # Put into STREAMM object
                 #   nanoBonds.put(bnd)                        # add if not in container

                allBonds.append([gid_i, gid_j])               # Store bond info (for MD)
                inghList.append(gid_j)                        # Store indiv. neighbors
                allDist.append(distSqr)                       # Store distances for diagnostics

        numNeigh = len(inghList)                       # Find #-neighbors
        localTotalNeighs = localTotalNeighs + numNeigh # Total neigh count for final output
        neighInfoLine = [gid_i, numNeigh] + inghList   # Build line with gid, #, {neigh1, neigh2..}
        allNeighbors.append(neighInfoLine)             # Build neighbor info


    p.barrier()
    end_time = time.time()

    # Calculation times and close dist-s
    proc_calc_time = end_time - start_time
    tmp = p.allReduceSum(proc_calc_time)
    avg_calc_time = tmp/float(size)
    totalDist = p.allReduceSum(len(allDist))

    # Timing/status info
    statsName="tmp--" + str(rank) + ".txt"
    statsIO=IO(statsName, False)

    if rank == 0:
        statsIO.prtFile(" ")
        statsIO.prtFile("     Average calc-time/proc = " + str(avg_calc_time))
        statsIO.prtFile("Total close distances found = " + str(totalDist))
    p.barrier()

    numpoints = len(globalList) # Total number of points (for diag)
    for proc in range(size):
        if proc == rank:
            statsIO.prtFile(" ")
            statsIO.prtFile("Time of processor " + str(proc) + " " + str(end_time - start_time) + " seconds")
            statsIO.prtFile("Number of close distances found = " 
                            + str(len(allDist)) + " for " + str(numpoints)
                            + " points on proc " + str(proc))
        p.barrier()
    p.barrier()

    del(statsIO)

    # Merge files by hand
    if rank == 0:
        os.system("cat tmp--*.txt > stats.txt")
        os.system("rm -rf tmp--*.txt")
    p.barrier()


    # Gather and sort neighbors on proc-0
    if rank == 0:
        print "Begin gathering neighbor list"
    numTotalNeighbors = p.allReduceSum(localTotalNeighs)
    allNeighbors      = p.gatherList(allNeighbors)
    if rank == 0:
        print "Begin sorting neighbor list"
    allNeighbors.sort()

    allBonds = p.gatherList(allBonds)
    allBonds.sort()

    return allDist, allNeighbors, numTotalNeighbors, allBonds



# Get useful methods
g=Geometry()

# Get comm object
p = mpiNREL.getMPISerialObject()
rank = p.getRank()
size = p.getCommSize()

# Sync random stream for tests
random.seed(0)

#############################################################################################
#
# Generate list of  [gid, x,y,z] points ---> allPoints list
# copied on all processors. # Global pt index included
# UNITS distance --> A [0.1 nm]
#
#############################################################################################


rad_avg = 15.0     # Moving to use this to make
rad_sig = 0.025    # initial cubic grid thats used to pass to MD
ptcl_dist = 10.0   # 


# NOTE: this is an arbitrary choice
cutoff_dist   = math.sqrt(2)*((2.00 * rad_avg) + ptcl_dist)
nQD_spacing   = (2.00 * rad_avg) + ptcl_dist

# Makes number of points along side grid
nQDs = [5, 5, 5]
sysLs = setSystemSizes(nQDs, nQD_spacing)

# Set boundary condition list (used in PBC calc)
xL = sysLs[0]
yL = sysLs[1]
zL = sysLs[2]
bcLs = [0.0, yL[1]-yL[0], zL[1]-zL[0] ] # x is NOT PB

# Status info
if rank == 0:
    print "rad_avg     = ", rad_avg
    print "rad_sig     = ", rad_sig
    print "cutoff_dist = ", cutoff_dist
    print " "
    print "bcLs  = ", bcLs
    print "sysLs = ", sysLs
    print " "

#############################################################################################
# Generate initial points

p.barrier()

tmpPoints = setRandomGridPoints(nQDs, rad_avg, ptcl_dist) # Generate global pts lst
allPoints = p.bcast(tmpPoints)                            # ensures rand pts same across procs

numpoints = len(allPoints) # Total number of points (for diag)
if numpoints == 0:
    print "Number of points generated == 0"
    sys.exit(3)

myPoints  = p.splitListOnProcs(allPoints) # Split elements on across processors
p.barrier()
#############################################################################################


###################################################################
# Store particles in struc and find neighbors/bonds
# ptPos includes index eg. [1, 3.4, 1.2, -2.0]
    
nanoPtcls = ParticleContainer()
nanoBonds = BondContainer()

for ptPos in allPoints:

    axisLoc = [bcLs[1]/2.0, bcLs[2]/2.0]
    cylinderRadius = 120.0
    inside  = isPtInCylinder(ptPos[1:], axisLoc, "yz", cylinderRadius)
    
    if inside:
        pt = Particle(ptPos[1:], type="QDbig",   mass=2.0)
        tagsD = {"molnum":1, "QDSizeAvg":20.0}
    else:
        pt = Particle(ptPos[1:], type="QDsmall", mass=1.0)
        tagsD = {"molnum":1, "QDSizeAvg":15.0}


    pt.setTagsDict(tagsD)
    nanoPtcls.put(pt)

p.barrier()

# Get bond info
allDist, allNeighbors, numTotalNeighbors, allBonds = \
      setQDNeighbors(myPoints, allPoints, cutoff_dist, bcLs)

if rank == 0:
    print "Begin building bond container"

# Put bonds into container
for bond in allBonds:
    
    if not nanoBonds.hasBond(bond):                   # Check if bond is unique
        bnd = Bond(bond[0], bond[1], type="normBond") # Put into STREAMM object
        nanoBonds.put(bnd)                            # add if not in container

p.barrier()

if rank == 0:
    print "length of bond container = ", len(nanoBonds)
###############################################################################


############################################################################
# Write LAMMPS file with atoms/bonds

strucQD = StructureContainer(nanoPtcls, nanoBonds)
# strucQD = StructureContainer(nanoPtcls)

boxSizes = sysLs
strucQD.setBoxLengths(boxSizes)

small_rad_avg   = 15.0
big_bond_min    =  (2.0*      rad_avg) + ptcl_dist
small_sigma_min = ((2.0*small_rad_avg) + ptcl_dist) / pow(2, 1/6.)
big_sigma_min   = ((2.0*      rad_avg) + ptcl_dist) / pow(2, 1/6.)


ptclParamMap = {("QDsmall", "epsilon"):1.0, ("QDsmall", "sigma"):small_sigma_min,
                 ("QDbig",   "epsilon"):1.0, ("QDbig",  "sigma"):big_sigma_min}

# Just one bond type for now
bondParamMap = {("normBond", "Kenergy"):1.0, ("normBond", "r0"):big_bond_min}

if rank == 0:
    # nanoPtcls.scatterPlot()
    strucQD.dumpLammpsInputFile("qd.data", ptclParamMap, bondParamMap)

# For format of test check
os.system("cat qd.data")
############################################################################
