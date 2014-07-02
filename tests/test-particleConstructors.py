#!/usr/bin/env python
"""
Class data structures for atomic data
"""

import copy

import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

print "************************************************************************************"
print " This test illustrates the different ways Particle constructors can be used  "
print " by using default values etc "
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
p4 = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)


p5 = Particle( type="Si", charge=2.0, pos=[0.2, 1.3, 54.0], mass=100.0)
p6 = Particle( pos=[0.2, 1.3, 54.0], mass=100.0)
p7 = Particle()
p8 = Particle( [0.01, 2.4, 10.2], "Co", mass=134.0, charge=0.001 )
p9 = Particle( [0.2, 1.3, 54.0], mass=100.0)
p10 = Particle( [0.2, 1.3, 54.0] )

# The following fail with argument checking native python functionality
# p9  = Particle( [0.01, 2.4, 10.2], "Co", mass=134.0, 0.001 )
# p10 = Particle( [0.01, 2.4, 10.2], mass=134.0, "Co" )

print "p5  = ", p5.__dict__
print "p6  = ", p6.__dict__
print "p7  = ", p7.__dict__
print "p8  = ", p8.__dict__
print "p9  = ", p9.__dict__
print "p10 = ", p10.__dict__
