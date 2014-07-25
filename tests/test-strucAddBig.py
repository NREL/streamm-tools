#!/usr/bin/env python
"""
Class data structures for atomic data
"""
import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from bonds import Bond
from bonds import BondContainer

from structureContainer import StructureContainer


print "************************************************************************************"
print " This is a speed/memory test for adding StructureContainer objects together "
print "************************************************************************************ \n"

def runStrucAdd(numPtcls1, numPtcls2):
    """
    Creates structures with ptcls and bonds between consecutive IDs
    Adds and reports timing

    Args:
        numPtcls1 (int): # of ptcls in first  struc
        numPtcls2 (int): # of ptcls in second struc
    Returns:
           (initTime, addTime, compressTime)
    """

    atoms1 = ParticleContainer()
    bonds1 = BondContainer()
    atoms2 = ParticleContainer()
    bonds2 = BondContainer()

    start_time = time.time()
    for n in range(numPtcls1):
        pobj1 = Particle( [1.1, 1.1, 1.1], "Si", 2.0, float(n+1))
        atoms1.put(pobj1)

    for n in range(1,numPtcls1):
        bobj1 = Bond( n, n+1, 1.233, "hooke-1")
        bonds1.put(bobj1)

    for n in range(numPtcls2):
        pobj2 = Particle( [2.2, 2.2, 2.2], "C",  1.0, float(n+1))
        atoms2.put(pobj2)

    for n in range(1,numPtcls2):
        bobj2 = Bond( n, n+1, 1.233, "hooke-2")
        bonds2.put(bobj2)
    end_time = time.time()

    initTime = end_time-start_time

    polymer1 = StructureContainer(atoms1,bonds1)  # Complete structure 1
    polymer2 = StructureContainer(atoms2,bonds2)  # Complete structure 2
    del atoms1, atoms2
    print "\n Cleaning memory for initial objects \n" 

    start_time = time.time()
    polymer1 += polymer2
    end_time = time.time()
    addTime = end_time - start_time

    start_time = time.time()    
    polymer1.compressPtclIDs()
    end_time = time.time()
    compressTime = end_time - start_time

    return (numPtcls1, numPtcls2, initTime, addTime, compressTime)


def convertWithTolerance(x, sig):
    tolerance = pow(10,sig)
    tmp = x * tolerance
    tmp = round(tmp)
    tmp = tmp / tolerance
    return tmp


def printTiming(timeTuple, fileObj):
    """
    Take timing info and round within a certain tolerance.
    Tests will pass if within this tolerance

    Args:
        timeTuple (5-tuple): (#-ptcls1, #-ptcls2, initTime, addTime, compressTime)
        fileObj            : output file object
    """

    numPtcls1    = timeTuple[0]
    numPtcls2    = timeTuple[1]
    initTime     = timeTuple[2]
    addTime      = timeTuple[3]
    compressTime = timeTuple[4]

    fileObj.write("--------------------------------------------- \n")    
    fileObj.write("#-of-particles  = " + str(numPtcls1) + " " + str(numPtcls2) + "\n")
    fileObj.write("  Init time     = " + str(initTime) + "\n")
    fileObj.write("  Add  time     = " + str(addTime) + "\n")
    fileObj.write("  Compress time = " + str(compressTime) + "\n")
    fileObj.write("--------------------------------------------- \n")


timeTuple1 = runStrucAdd(100,  100)
timeTuple2 = runStrucAdd(1000, 1000)
timeTuple3 = runStrucAdd(5000, 5000)
timeTuple4 = runStrucAdd(5000, 100)

fobj = open("test-strucAddBig-timing.dat", 'w')
printTiming(timeTuple1, fobj)
printTiming(timeTuple2, fobj)
printTiming(timeTuple3, fobj)
printTiming(timeTuple4, fobj)
fobj.close()
