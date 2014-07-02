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
print " This test illustrates the different ways ParticleContainer can set Particle "
print " object data"
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
p4 = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)

p5 = Particle( [9.0, -8.0, -80.1], "Ar",  1.0, 4.0)
p6 = Particle( [90.0, -0.0, -80.1], "Zr",  1.0, 4.0)

atoms1 = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)
atoms1.put(p4)

print "Initial atoms1 state from an empty container, populated with put()", atoms1, "\n"

tag = 4
print "Replacing ID tag", tag, "with p5 particle data \n"
atoms1[4]=p5
print "atoms1 state", atoms1, "\n"

tag = 5
print "Setting ID tag that does not exist", tag, "with p6 particle data ---> atoms1[5]=p6 \n"
atoms1[5]=p6
print "atoms1 state", atoms1, "\n"


