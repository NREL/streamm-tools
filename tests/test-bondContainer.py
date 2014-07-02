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
b4 = Bond( 3, 8, 0.5,   "stiff")

bonds = BondContainer()
bonds.put(b1)
bonds.put(b2)
bonds.put(b3)
bonds.put(b4)
del b1, b2, b3, b4
print "\n Cleaning memory for initial objects \n" 

print "Check for pre-existing bond"
if bonds.hasBond([2,1]):
    print "bond 1--2 exists"
else:
    print "bond 1--2 does NOT exists"

print "Check for pre-existing bond"
if bonds.hasBond([2,3]):
    print "bond 2--3 exists"
else:
    print "bond 2--3 does NOT exist"


print " "
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

print "-----------------------------------------------------"
print "bonds.getTypeInfoDict() ", bonds.getTypeInfoDict()
print "-----------------------------------------------------"

print " "
print "Here's an issue... globalID's are assigned by ParticleContainer class "
print "Setting Bond object with gid's of Particles needs to be done after PC class set "
print "or there could be conflicts \n"

