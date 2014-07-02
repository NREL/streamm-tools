#!/usr/bin/env python
"""
Class data structures for atomic data
"""
import os, sys, math, random, time

from bonds import Bond
from bonds import BondContainer

print "************************************************************************************"
print " This test shows various operators within Bond and BondContainer classes    "
print "************************************************************************************ \n"

print "Testing Bond and BondContainer"
b1 = Bond( 1, 2, 1.233, "hooke")
b2 = Bond( 3, 4, 0.5,   "hooke")
b3 = Bond( )

print "b1 = ", b1.__dict__
print "b2 = ", b2.__dict__
print "b3 = ", b3.__dict__

bonds = BondContainer()
bonds.put(b1)
bonds.put(b2)
bonds.put(b3)
del b1, b2, b3
print "\n Cleaning memory for initial objects \n" 

print "Bonds: ", bonds
print " "

print "Test iterator for bond container "
for gid, bondObj in bonds:
    print "bondID = ", gid, "bond object ", bondObj
print " "

searchID = 2
print "Testing 'in' operator for searchID ", searchID, " in bonds"
if searchID in bonds:
    print "ID ",searchID, " in bonds"
else:
    print "nope"


print " "
print "Here's an issue... globalID's are assigned by ParticleContainer class "
print "Setting Bond object with gid's of Particles needs to be done after PC class set "
print "or there could be conflicts \n"

