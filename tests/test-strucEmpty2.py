#!/usr/bin/env python
"""
Class data structures for atomic data
"""
import os, sys, math, random, time, copy

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

print "------------------------------------------------------------------------------------ \n"

print "Testing struc add --> empty (polymer2) += non-empty (polymer1)"
polymer2 += polymer1
print "polymer1 (non-empty)= ", polymer1
print "------------------------------------------------------------------------------------ \n"
print "polymer2 (empty)    = ", polymer2


polymer3 = copy.deepcopy(polymer1)

polymer1.ptclC[2].position = [11111.111, 22222.2222, 33333.333]



print "------------------------------------------------------------------------------------ \n"
print polymer1
print "------------------------------------------------------------------------------------ \n"
print polymer2
print "------------------------------------------------------------------------------------ \n"
print polymer3
print "------------------------------------------------------------------------------------ \n"
