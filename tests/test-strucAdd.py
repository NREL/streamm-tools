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
print " This test shows how to add StructureContainer objects together "
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)

b1 = Bond( 1, 2, 1.233, "hooke")
b2 = Bond( 2, 3, 0.500, "hooke")

atoms1   = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)

bonds1   = BondContainer()
bonds1.put(b1)
bonds1.put(b2)

polymer1 = StructureContainer(atoms1, bonds1)  # Complete structure 1 completely

p1other = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)
p2other = Particle( [50.0, 0.3, -0.1], "Ar", 2.0, 2.34)

b1other = Bond( 1, 2, 1.233, "hooke")    # Correct ptclIDs for second structure

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
polymer1 += polymer2
print "polymer1 = ", polymer1

print "-------------------- Results (check above) --------------------"
print " 1---b1---2---b2---3 + 1---b1----2   should go to"
print " "
print " 1---b1---2---b2---3   4---b3----5 \n"

print " "
print "------ After adding polymer 2 should be unchanged--------------"
print "polymer2 = ", polymer2
