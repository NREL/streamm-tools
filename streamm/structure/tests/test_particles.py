# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0


'''
Unit tests for the particles module
'''

from __future__ import division, unicode_literals

import logging
logger = logging.getLogger(__name__)

import unittest

try:
    import streamm.structure.particles as particles
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    import particles
    
    
class TestParticle(unittest.TestCase):
    def setUp(self):
        self.part = particles.Particle("genpart")
        
    def test_checktype(self):
        self.assertEqual(self.part.type,"genpart")
        
    def test_str(self):
        self.assertEqual(str(self.part)," genpart")

    def tearDown(self):
        del self.part 
        self.part = None

class TestForceField(unittest.TestCase):
    def setUp(self):
        self.ff = particles.ForceField("CC*",'C24')
        self.ff.charge = 0.756
        self.ff.mass  = 12.0
        self.ff.lammps_index = 23
        self.ff.gromacs_index = 863
    
        
    def test_checktype(self):
        self.assertEqual(self.ff.type,"CC*")
        
    def test_checktlabel(self):
        self.assertEqual(self.ff.label,"C24")
        
    def test_str(self):
        self.assertEqual(str(self.ff)," C24 (CC*)")

    def test_charge(self):
        self.assertEqual(self.ff.charge,0.756)
        
    def test_mass(self):
        self.assertEqual(self.ff.mass,12.0)
        
    def test_lammps_index(self):
        self.assertEqual(self.ff.lammps_index,23)
        
    def test_gromacs_index(self):
        self.assertEqual(self.ff.gromacs_index,863)
        
    def tearDown(self):
        del self.ff 


if __name__ == '__main__':
    unittest.main()
