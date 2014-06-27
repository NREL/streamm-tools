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
print " This test illustrates the search capability for multiple tags                      "
print " This can be combined new class method that returns iterator over search results    "
print "************************************************************************************ \n"

p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
tagsD = {"molnum":1,"ringnum":4}
p1.setTagsDict(tagsD)

p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
tagsD = {"molnum":2,"ringnum":4}
p2.setTagsDict(tagsD)

p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
tagsD = {"molnum":1, "ringnum":4}
p3.setTagsDict(tagsD)

p4 = Particle( [0.0, 2.3, -20.1], "Si",  1.0, 2.34)
tagsD = {"molnum":2,"ringnum":4}
p4.setTagsDict(tagsD)

p5 = Particle( [1.0, 2.3, -20.1], "C",  1.0, 5.34)
tagsD = {"molnum":2,"ringnum":2}
p5.setTagsDict(tagsD)

p6 = Particle( [8.0, 2.3, -20.1], "Si",  1.0, 8.00)
tagsD = {"molnum":2,"ringnum":3}
p6.setTagsDict(tagsD)

p7 = Particle( [8.0, 2.3, -20.1], "O",  1.0, 8.00)
tagsD = {"molnum":2,"ringnum":3}
p7.setTagsDict(tagsD)

p8 = Particle( [8.0, 2.3, -20.1], "O",  1.0, 8.00)
tagsD = {"molnum":2,"ringnum":3}
p8.setTagsDict(tagsD)

p9 = Particle( [8.0, 2.3, -20.1], "H",  1.0, 8.00)
tagsD = {"molnum":2,"ringnum":3}
p9.setTagsDict(tagsD)

"""
print "p1 --> ", p1.__dict__
print "p1 has molnum=1 ",   p1.isTagEqualTo("molnum", 1)
print "p1 has molnum=100 ", p1.isTagEqualTo("molnum", 100)
print "p3 has molnum=2 ",   p3.isTagEqualTo("molnum", 2)
"""

atoms1 = ParticleContainer()
atoms1.put(p1)
atoms1.put(p2)
atoms1.put(p3)
atoms1.put(p4)
atoms1.put(p5)
atoms1.put(p6)
atoms1.put(p7)
atoms1.put(p8)
atoms1.put(p9)
del p1, p2, p3, p4, p5, p6, p7, p8, p9

print "atoms1 initially contains all of the following"
print atoms1
# This should fail
# print "Testing callable method ", atoms1()
print " "
# -----------------------------------------------------------


print "Print only atoms1 with search.... \n"
subList = atoms1.getParticlesWithTags({'molnum':2, 'type':"Si"})
for id, ptcl in atoms1(subList):
    print "id = ", id, " ptclObj = ", ptcl.__dict__
print "--------------------------------------------------------------------------------- \n"

print "Print only atoms1 with search.... \n"
subList = atoms1.getParticlesWithTags({'molnum':2, 'type':"Si",'ringnum':4})
for id, ptcl in atoms1(subList):
    print "id = ", id, " ptclObj = ", ptcl.__dict__
print "--------------------------------------------------------------------------------- \n"

print "Print only atoms1 with search.... \n"
subList = atoms1.getParticlesWithTags( {'type':"Si", 'ringnum':[3, 4]} )
for id, ptcl in atoms1(subList):
    print "id = ", id, " ptclObj = ", ptcl.__dict__
print "--------------------------------------------------------------------------------- \n"    

print "Print only atoms1 with search.... \n"
subList = atoms1.getParticlesWithTags( {'type':["Si","O"], 'ringnum':[2, 3, 4]} )
for id, ptcl in atoms1(subList):
    print "id = ", id, " ptclObj = ", ptcl.__dict__



