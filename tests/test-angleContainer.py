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


print "************************************************************************************"
print " This test shows various operators within Angle and AngleContainer classes          "
print " Also shows memory management structure and access methods                          "
print "************************************************************************************ \n"

p1 = Particle( [1.1, 1.1, 1.1], "Si", 2.0, 1.23)
p2 = Particle( [2.2, 2.2, 2.2], "Si",  1.0, 2.34)
p3 = Particle( [3.3, 3.3, 3.3], "Si",  1.0, 2.34)
#
p4 = Particle( [4.4, 4.4, 4.4], "C",  1.0, 2.34)  # ptcl 1
p5 = Particle( [5.5, 5.5, 5.5], "C",  1.0, 2.34)  # ptcl 2
p6 = Particle( [6.6, 6.6, 6.6], "C",  1.0, 2.34)  # ptcl 3

b1 = Bond( 1, 2, 1.11, "hooke")
b2 = Bond( 2, 3, 2.22, "hooke")
#
b3 = Bond( 1, 2, 3.33, "hooke")
b4 = Bond( 2, 3, 4.44, "hooke")

a1 = Angle(1, 2, 3, 1.11, "harmonic")
#
a2 = Angle(1, 2, 3, 2.22, "harmonic")

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

angles1 = AngleContainer()
angles2 = AngleContainer()
angles1.put(a1)
angles2.put(a2)

del p1, p2, p3, p4, p5, p6, b1, b2, b3, b4, a1, a2
print "\n Cleaning memory for initial objects \n" 

print "angles1 container"
print angles1

print " "
print "angles2 container"
print angles2

print "Testing 'in' operator (1 in angles1)"
if (1 in angles1):
    print "angles1 contains gid 1"
else:
    print "key not found in angles1"

print "Testing 'in' operator (5 in angles1)"
if (5 in angles1):
    print "angles1 contains gid 5"
else:
    print "key not found in angles1"

print " "
angles1 += angles2
print "Will print the new angles1 after adding angles1 += angles2"
print angles1


print "Check for pre-existing angle"
if angles1.hasAngle([3,2,1]):
    print "angle 1--2--3 exists"
else:
    print "angle 1--2--3 does NOT exists"

print "Check for pre-existing angle"
if angles1.hasAngle([2,3,1]):
    print "angle 2--3--1 exists"
else:
    print "angle 2--3--1 does NOT exists"    

