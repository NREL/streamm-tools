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
print " This shows operations on empty StructureContainer objects for 'robustness'"
print "************************************************************************************ \n"

atoms1   = ParticleContainer()
bonds1   = BondContainer()

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
b1 = Bond( 1, 2, 1.233, "hooke")
b2 = Bond( 2, 3, 0.500, "hooke")

atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)
bonds1.put(b1)
bonds1.put(b2)

atoms2   = ParticleContainer()
bonds2   = BondContainer()

polymer1 = StructureContainer(atoms1, bonds1)  # Non-empty structure 1
polymer2 = StructureContainer(atoms2, bonds2)  # Empty     structure 2

# ------------------------------------------------------------------------------------------------------------------------
 
print "Non-empty container", polymer1
print "Empty container",     polymer2

print "Testing empty dump"
polymer2.dump("empty")

print "Testing empty compressPtclIDs"
polymer2.compressPtclIDs()

print "Testing empty getSubStructure"
subStruc = polymer2.getSubStructure([])
print "subStruc = ", subStruc


print "Testing empty getPtclPositions"
posList = polymer2.getPtclPositions()
print "posList = ", posList

print "------------------------------------------------------------------------------------ \n"

print "Testing struc add --> non-empty (polymer1) += empty (polymer2)"
polymer1 += polymer2
print "polymer1 = ", polymer1
print "------------------------------------------------------------------------------------ \n"
print "polymer2 = ", polymer2

print "------------------------------------------------------------------------------------ \n"

print "Testing struc add --> empty (polymer2) += non-empty (polymer1)"
polymer2 += polymer1
print "polymer2 = ", polymer2
print "------------------------------------------------------------------------------------ \n"
print "polymer1 = ", polymer1


print "------------------------------------------------------------------------------------ \n"

print "Testing empty getSubStructure with non-zero id list (should return ERROR)"
polymer3 = StructureContainer()
subStruc = polymer3.getSubStructure([1,2])
print "subStruc = ", subStruc
