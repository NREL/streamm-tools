#!/usr/bin/env python
"""
Class data structures for atomic data
"""
import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from bonds import Bond
from bonds import BondContainer

from dihedrals import Dihedral
from dihedrals import DihedralContainer


print "************************************************************************************"
print " This test shows various operators within Dihedral and DihedralContainer classes    "
print " Also shows memory management structure and access methods                          "
print "************************************************************************************ \n"

p1 = Particle( [1.1, 1.1, 1.1], "Si", 2.0, 1.23)
p2 = Particle( [2.2, 2.2, 2.2], "Si",  1.0, 2.34)
p3 = Particle( [3.3, 3.3, 3.3], "Si",  1.0, 2.34)
p4 = Particle( [4.4, 4.4, 4.4], "Si",  1.0, 2.34)
#
p5 = Particle( [5.5, 5.5, 5.5], "C",  1.0, 2.34)
p6 = Particle( [6.6, 6.6, 6.6], "C",  1.0, 2.34)
p7 = Particle( [7.7, 7.7, 7.7], "C",  1.0, 2.34)
p8 = Particle( [8.8, 8.8, 8.8], "C",  1.0, 2.34)

b1 = Bond( 1, 2, 1.11, "hooke")
b2 = Bond( 2, 3, 2.22, "hooke")
#
b3 = Bond( 1, 2, 3.33, "hooke")
b4 = Bond( 1, 2, 4.44, "hooke")

d1 = Dihedral(1, 2, 3, 4, 1.11, "stiff")
#
d2 = Dihedral(1, 2, 3, 4, 2.22, "stiff")

atoms1 = ParticleContainer()
atoms2 = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)
#
atoms2.put(p4)
atoms2.put(p5)
atoms2.put(p6)

bonds1 = BondContainer()
bonds2 = BondContainer()
bonds1.put(b1)
bonds1.put(b2)
bonds2.put(b3)
bonds2.put(b4)

dihedrals1 = DihedralContainer()
dihedrals2 = DihedralContainer()
dihedrals1.put(d1)
dihedrals2.put(d2)

del p1, p2, p3, p4, p5, p6, b1, b2, b3, b4, d1, d2
print "\n Cleaning memory for initial objects \n" 

print "dihedral1 container"
print dihedrals1

print " "
print "dihedral2 container"
print dihedrals2

print "Testing 'in' operator (1 in dihedrals1)"
if (1 in dihedrals1):
    print "dihedrals1 contains gid 1"
else:
    print "key not found in dihedrals1"

print "Testing 'in' operator (5 in dihedrals1)"
if (5 in dihedrals1):
    print "dihedrals1 contains gid 5"
else:
    print "key not found in dihedrals1"

print " "
dihedrals1 += dihedrals2
print "Will print the new dihedrals1 after adding dihedrals1 += dihedrals2"
print dihedrals1


print "Check for pre-existing dihedral"
if dihedrals1.hasDihedral([4,3,2,1]):
    print "dihedral 1--2--3--4 exists"
else:
    print "dihedral 1--2--3--4 does NOT exists"

print "Check for pre-existing dihedral"
if dihedrals1.hasDihedral([2,3,1,4]):
    print "dihedral 2--3--1--4 exists"
else:
    print "dihedral 2--3--1--4 does NOT exists"
