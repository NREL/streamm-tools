#!/usr/bin/env python
"""
Class data structures for atomic data
"""
import os, sys, math, random, time

from particles import Particle
from particles import ParticleContainer

from bonds import Bond
from bonds import BondContainer

from structureContainer import StructureContainer
from simulation import Simulation


print "************************************************************************************"
print " This test shows basics of the Simulation classes"
print " Highlights differences between setting reference to the structure"
print "  container and a deepcopy"
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)

atoms1   = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)

polymer1 = StructureContainer(atoms1,verbose=False)  # Complete structure 1 completely
polymer2 = StructureContainer(verbose=False)
polymer2 += polymer1

del p1, p2, p3, atoms1
print "\n Cleaning memory for initial objects \n" 

print "-------------------- Before adding --------------------"
print "polymer1 = ", polymer1
print "polymer2 = ", polymer2
print " "

simObj1 = Simulation("GaussianJuly2014")
simObj2 = Simulation("Lammps012038")

simObj1.setStructureContainer(polymer1)
simObj2.copyStructureContainerInto(polymer1)

del polymer1.ptclC[1]

print "-------------------- After changing polymer1 struc --------------------"
print "polymer1 = ", polymer1
simObj1.printStruc()
simObj2.printStruc()

#simC = SimulationContainer(verbose=False)
#simC[1] = (polymer1, simObj1, polymer2)
#simC[2] = (polymer1, simObj2, polymer2)
#del simObj1, simObj2
#print "\n Cleaning memory for simulation objects \n"
#print "simC = ", simC
#simTuple1 = simC[1]
#print "simTuple1 init struc  = ", simTuple1[0]
#print "simTuple1 simulation  = ", simTuple1[1].__dict__
#print "simTuple1 final struc = ", simTuple1[2]
