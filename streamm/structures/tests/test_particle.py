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
import os

import streamm.structures.particle as particle
    
    

from streamm_testutil import *

class TestParticle(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.part = particle.Particle(label="C137")
        
    def test_checktype(self):
        self.assertEqual(self.part.type,"atom")
        
    def test_str(self):
        self.assertEqual(str(self.part),"atom[None] C137 (C137)")
        
    def test_element_symbol(self):
        self.part.set_element(symbol='S')
        self.assertEqual(self.part.mass,32.065)
        self.assertEqual(self.part.bonded_radius,0.88)
        self.assertEqual(self.part.nonbonded_radius,1.8)

    def test_element_mass(self):
        self.part.set_element(mass=1.011)
        self.assertEqual(self.part.mass,1.00794)
        self.assertEqual(self.part.bonded_radius,0.53)
        self.assertEqual(self.part.nonbonded_radius,1.2)

    def test_element_number(self):
        self.part.set_element(number=8)
        self.assertEqual(self.part.mass,15.9994)
        self.assertEqual(self.part.bonded_radius,0.48)
        self.assertEqual(self.part.nonbonded_radius,1.52)
        

    @tearDown_streamm 
    def tearDown(self):
        del self.part 


class TestParticleAtom(unittest.TestCase):
    
    
    @setUp_streamm 
    def setUp(self):
        self.atom = particle.Particle(symbol='Au')
            
        self.atom.mol = 112
        self.atom.ring  = 3
        self.atom.residue  = 387
        self.atom.resname = "MET1"
        self.atom.qgroup = 23469

    def test_checktype(self):
        self.assertEqual(self.atom.type,"atom")

    def test_properties(self):
        self.assertEqual(self.atom.symbol,"Au")
        
        self.assertEqual(self.atom.mol,112)
        self.assertEqual(self.atom.ring ,3)
        self.assertEqual(self.atom.residue ,387)
        self.assertEqual(self.atom.resname,"MET1")
        self.assertEqual(self.atom.qgroup,23469)
        #
        self.assertEqual(self.atom.mass,196.966569)
        self.assertEqual(self.atom.bonded_radius,1.74)
        self.assertEqual(self.atom.nonbonded_radius,1.66)
        #
        # Test element setting from pymatgen
        self.assertEqual(self.atom.element.atomic_mass,196.966569)
        self.assertEqual(self.atom.element.covalent_radius,1.74)
        self.assertEqual(self.atom.element.vdw_radius,1.66)
    
    def test_str(self):
        self.assertEqual(str(self.atom),"atom[None] Au (Au)")

    @tearDown_streamm 
    def tearDown(self):
        del self.atom 
        
    
if __name__ == '__main__':

    unittest.main()    

        
                
