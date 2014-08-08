#!/usr/bin/env python
"""
Class data structures for atomic data
"""

import copy

# Old setup testing replace IDs
diagramBefore="""  
1--'b1'---2
 -    'a1'|
  -       'b2'
   'b4'   |
     -    |'a2'
       -  3---'b3'---4
"""

diagramAfter="""
5--'b1'---2
 -    'a1'|
  -       'b2'
   'b4'   |
     -    |'a2'
       -  3---'b3'---4
"""


import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from bonds import Bond
from bonds import BondContainer

from angles import Angle
from angles import AngleContainer

from structureContainer import StructureContainer

print "************************************************************************************"
print " This test shows how to set up Structure container with Particle/Bond-Containers \n"
print " Shows how ID changed in Structure propagate to values set in BondContainer for its "
print "    held particle ID values \n"
print " Illustrates how a substructure method could return subgroup"
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
p4 = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)
p5 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)

b1 = Bond( 5, 2, 1.233, "hooke")  # This setup mimics earlier state of test
b2 = Bond( 2, 3, 0.500, "hooke")  # Still matches with diagram
b3 = Bond( 3, 4, 2.301, "hooke")
b4 = Bond( 5, 3, 0.828, "hooke")

a1 = Angle(1, 2, 3, 1.111, "harmonic")
a2 = Angle(2, 3, 4, 2.222, "harmonic")

atoms1 = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)
atoms1.put(p4)
atoms1.put(p5)

del atoms1[1]  # This is to mimic state of container for an older test

bonds = BondContainer()
bonds.put(b1)
bonds.put(b2)
bonds.put(b3)
bonds.put(b4)

angles = AngleContainer()
angles.put(a1)
angles.put(a2)

del p1, p2, p3, p4, b1, b2, b3, b4, a1, a2
polymer1 = StructureContainer(atoms1, bonds, angles)
del atoms1, bonds, angles

print "********************************************************** \n"
print polymer1
print diagramAfter
print "********************************************************** \n"

print "-------------------------------------------------------------------------------- \n"

print "********************************************************** \n"
print "Testing polymer1.getSubStructure([5,2])"
print "   currently ID's are reassigned in substructure \n"
subpolymer = polymer1.getSubStructure([5,2])
print subpolymer
print diagramAfter
print "********************************************************** \n"

print "********************************************************** \n"
print "polymer1 Struture after returning sub-structure"
print polymer1
print "********************************************************** \n"

print "-------------------------------------------------------------------------------- \n"

print "********************************************************** \n"
print "Testing polymer1.getSubStructure([2,3,4])"
print "   currently ID's are reassigned in substructure \n"
subpolymer = polymer1.getSubStructure([2,3,4])
print subpolymer
print diagramAfter
print "********************************************************** \n"
