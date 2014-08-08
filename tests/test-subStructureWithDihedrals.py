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

from dihedrals import Dihedral
from dihedrals import DihedralContainer


from structureContainer import StructureContainer

print "************************************************************************************"
print " Illustrates how a substructure method returns subgroup with all 4 types"
print "************************************************************************************ \n"

p1 = Particle( [1.1, 1.1, 1.1], "Si", 2.0, 1.23)
p2 = Particle( [2.2, 2.2, 2.2], "Si",  1.0, 2.34)
p3 = Particle( [3.3, 3.3, 3.3], "Si",  1.0, 2.34)
p4 = Particle( [4.4, 4.4, 4.4], "Si",  1.0, 2.34)
p5 = Particle( [5.5, 5.5, 5.5], "C",  1.0, 2.34)   # 1
p6 = Particle( [6.6, 6.6, 6.6], "C",  1.0, 2.34)   # 2
p7 = Particle( [7.7, 7.7, 7.7], "C",  1.0, 2.34)   # 3
p8 = Particle( [8.8, 8.8, 8.8], "C",  1.0, 2.34)   # 4

b1 = Bond( 1, 2, 1.11, "hooke")
b2 = Bond( 2, 3, 2.22, "hooke")
b3 = Bond( 3, 4, 3.33, "hooke")
b4 = Bond( 4, 5, 4.44, "hooke")
b5 = Bond( 5, 6, 5.55, "hooke")
b6 = Bond( 6, 7, 6.66, "hooke")
b7 = Bond( 7, 8, 7.77, "hooke")

a1 = Angle(1, 2, 3, 1.11, "harmonic")
a2 = Angle(2, 3, 4, 2.22, "harmonic")
a3 = Angle(3, 4, 5, 3.33, "harmonic")
a4 = Angle(4, 5, 6, 4.44, "harmonic")
a5 = Angle(5, 6, 7, 5.55, "harmonic")
a6 = Angle(6, 7, 8, 6.66, "harmonic")


d1 = Dihedral(1, 2, 3, 4, 1.11, "stiff")
d2 = Dihedral(2, 3, 4, 5, 2.22, "stiff")
d3 = Dihedral(3, 4, 5, 6, 3.33, "stiff")
d4 = Dihedral(4, 5, 6, 7, 4.44, "stiff")
d5 = Dihedral(5, 6, 7, 8, 5.55, "stiff")


atoms1 = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)
atoms1.put(p4)
atoms1.put(p5)
atoms1.put(p6)
atoms1.put(p7)
atoms1.put(p8)

bonds = BondContainer()
bonds.put(b1)
bonds.put(b2)
bonds.put(b3)
bonds.put(b4)
bonds.put(b5)
bonds.put(b6)
bonds.put(b7)

angles = AngleContainer()
angles.put(a1)
angles.put(a2)
angles.put(a3)
angles.put(a4)
angles.put(a5)
angles.put(a6)

dihedrals = DihedralContainer()
dihedrals.put(d1)
dihedrals.put(d2)
dihedrals.put(d3)
dihedrals.put(d4)
dihedrals.put(d5)


del p1, p2, p3, p4, p5, p6, p7, p8
del b1, b2, b3, b4, b5, b6, b7
del a1, a2, a3, a4, a5, a6
del d1, d2, d3, d4, d5

polymer1 = StructureContainer(atoms1, bonds, angles, dihedrals)
del atoms1, bonds, angles, dihedrals


print "********************************************************** \n"
print "polymer1 Structure before returning sub-structure"
print polymer1
print "********************************************************** \n"

print "-------------------------------------------------------------------------------- \n"

print "********************************************************** \n"
print "Testing polymer1.getSubStructure([2,3,4])"
print "   currently ID's are reassigned in substructure \n"
subpolymer = polymer1.getSubStructure([2,3,4])
print subpolymer
print "********************************************************** \n"
