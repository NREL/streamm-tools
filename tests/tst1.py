#!/usr/bin/env python
"""
Class data structures for atomic data
"""
import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from bonds import Bond
from bonds import BondContainer

print "************************************************************************************"
print " This test shows various operators within Particle and ParticleContainer classes    "
print " Also shows memory management structure and access methods"
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
p4 = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)

atoms1 = ParticleContainer()
atoms2 = ParticleContainer()

atoms1.put(p1)
atoms1.put(p2)
atoms2.put(p3)
atoms2.put(p4)

del p1, p2, p3, p4
print "\n Cleaning memory for initial objects \n" 


# This is a shallow copy
print "x = atoms1[1] is returns x as an effective 'reference' \n"
x = atoms1[1]
print "x = ", x.__dict__, "\n"
x.position=[1.0, 1.0, 1.0]
print "after changing with x.position = [1.0, 1.0, 1.0]"
print "x = ", x.__dict__, "\n"

# This value has been changed by code above
print "atoms1 has been changed"
print atoms1

print "before, atoms1--> ", atoms1, "\n"
del atoms1[2]
print "after 'del atoms1[2]' atoms1 --> ", atoms1, "\n"

print "Testing 'in' operator (1 in atoms1)"
if (1 in atoms1):
    print "atoms1 contains gid 1"
else:
    print "key not found in atoms1"

print "Testing 'in' operator (5 in atoms1)"
if (5 in atoms1):
    print "atoms1 contains gid 5"
else:
    print "key not found in atoms1"

print " "
atoms1 += atoms2
print "Will print the new atoms1 after adding atoms1 += atoms2"
print atoms1

print " "
print "Here's an issue... globalID's are assigned by ParticleContainer class "
print "Setting Bond object with gid's of Particles needs to be done after PC class set "
print "or there could be conflicts \n"

print "Testing Bond and BondContainer"
b1 = Bond( 1, 2, 1.233, "hooke")
b2 = Bond( 3, 4, 0.5,   "hooke")

print "b1 = ", b1.__dict__
print "b2 = ", b2.__dict__

bonds = BondContainer()
bonds.put(b1)
bonds.put(b2)
del b1, b2
print "\n Cleaning memory for initial objects \n" 

print "Bonds: ", bonds
print " "

print "Testing 'in' operator (5 in bonds"
if 5 in bonds:
    print "ID 1 in b1"
else:
    print "nope"

