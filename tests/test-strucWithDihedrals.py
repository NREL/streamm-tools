#!/usr/bin/env python
"""
Class data structures for atomic data
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
print " This test shows structureContainer functionality with angles included"
print "************************************************************************************ \n"

p1 = Particle( [1.1, 1.1, 1.1], "Si", 2.0, 1.23)
p2 = Particle( [2.2, 2.2, 2.2], "Si",  1.0, 2.34)
p3 = Particle( [3.3, 3.3, 3.3], "Si",  1.0, 2.34)
p4 = Particle( [4.4, 4.4, 4.4], "Si",  1.0, 2.34)
#
p5 = Particle( [5.5, 5.5, 5.5], "C",  1.0, 2.34)   # 1
p6 = Particle( [6.6, 6.6, 6.6], "C",  1.0, 2.34)   # 2
p7 = Particle( [7.7, 7.7, 7.7], "C",  1.0, 2.34)   # 3
p8 = Particle( [8.8, 8.8, 8.8], "C",  1.0, 2.34)   # 4


b1 = Bond( 1, 2, 1.11, "hooke")
b2 = Bond( 2, 3, 2.22, "hooke")
b3 = Bond( 3, 4, 3.33, "hooke")
#
b4 = Bond( 1, 2, 4.44, "hooke")
b5 = Bond( 2, 3, 5.55, "hooke")
b6 = Bond( 3, 4, 6.66, "hooke")

a1 = Angle(1, 2, 3, 1.11, "harmonic")
#
a2 = Angle(1, 2, 3, 2.22, "harmonic")


d1 = Dihedral(1, 2, 3, 4, 1.11, "stiff")
#
d2 = Dihedral(1, 2, 3, 4, 2.22, "stiff")

atoms1 = ParticleContainer()
atoms2 = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)
atoms1.put(p4)
#
atoms2.put(p5)
atoms2.put(p6)
atoms2.put(p7)
atoms2.put(p8)

bonds1 = BondContainer()
bonds2 = BondContainer()
bonds1.put(b1)
bonds1.put(b2)
bonds1.put(b3)
#
bonds2.put(b4)
bonds2.put(b5)
bonds2.put(b6)

angles1 = AngleContainer()
angles2 = AngleContainer()
angles1.put(a1)
angles2.put(a2)

dih1 = DihedralContainer()
dih2 = DihedralContainer()
dih1.put(d1)
dih2.put(d2)

polymer1 = StructureContainer(atoms1, bonds1, angles1, dih1)
polymer2 = StructureContainer(atoms2, bonds2, angles2, dih2)


del p1, p2, p3, p4, p5, p6, p7, p8
del b1, b2, b3, b4, b5, b6, a1, a2, d1, d2
del atoms1, atoms2, bonds1, bonds2, angles1, angles2, dih1, dih2
print "\n Cleaning memory for initial objects \n" 

print "-------------------- Initial structures --------------------"
print "polymer1 = ", polymer1
print "polymer2 = ", polymer2
print " "

print "-------------- After adding  (polymer1 += polymer2) ----------------"
polymer2 += polymer1
print "polymer2 = ", polymer2
