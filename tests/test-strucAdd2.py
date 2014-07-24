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


print "***************************************************************************************"
print " This test shows how to add StructureContainer objects together "
print " A special test to double-check re-labeling ... add a big container to a 'small' one "
print "*************************************************************************************** \n"

p1 = Particle( [1.1, 1.1, 1.1], "Si", 2.0, 1.23)
p2 = Particle( [2.2, 2.2, 2.2], "C",  1.0, 2.34)
p3 = Particle( [3.3, 3.3, 3.3], "C",  1.0, 2.34)

b1 = Bond( 1, 2, 1.111, "hooke")
b2 = Bond( 2, 3, 2.222, "hooke")

atoms1   = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)

bonds1   = BondContainer()
bonds1.put(b1)
bonds1.put(b2)

polymer1 = StructureContainer(atoms1, bonds1)  # Complete structure 1 completely

p1other = Particle( [1.11, 1.11, 1.11], "C",  1.0, 2.34)
p2other = Particle( [2.22, 2.22, 2.22], "Ar", 2.0, 2.34)

b1other = Bond( 1, 2, 1.1, "hooke-2")    # Correct ptclIDs for second structure

atoms2 = ParticleContainer()
atoms2.put(p1other)
atoms2.put(p2other)

bonds2   = BondContainer()
bonds2.put(b1other)

polymer2 = StructureContainer(atoms2, bonds2)  # Complete structure 1 completely

del p1, p2, p3, p1other, p2other, b1, b2, b1other, atoms1, atoms2, bonds1, bonds2
print "\n Cleaning memory for initial objects \n" 

print "-------------------- Before adding --------------------"
print "polymer1 = ", polymer1
print "polymer2 = ", polymer2
print " "

print "-------------------- After adding --------------------"
# polymer1 += polymer2
polymer2 += polymer1
print "polymer2 = ", polymer2
