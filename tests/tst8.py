#!/usr/bin/env python
"""
Class data structures for atomic data
"""

import copy

import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from bonds import Bond
from bonds import BondContainer

print "************************************************************************************"
print " This test illustrates the setting a ParticleContainer with pre-existing IDs "
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
p4 = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)
p5 = Particle( [90.0, -0.0, -80.1], "Zr",  1.0, 4.0)

idList = [2, 3, 6, 10]
print "Using idList = ", idList, " to initialize particle container"
atoms1 = ParticleContainer(idList)
print "Initial atoms1 state from a pre-set ID list ", atoms1, "\n"

atoms1[2 ]=p1
atoms1[3 ]=p2
atoms1[6 ]=p3
atoms1[10]=p4
print "atoms1 state from a pre-set ID list populated atoms1[..] = ... ", atoms1, "\n"

atoms1.put(p5)
print "atoms1 state after a put(p5) call ", atoms1, "\n"
