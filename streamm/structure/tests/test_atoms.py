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
    import streamm.structure.atoms as atoms
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    sys.path.append(os.path.join(os.path.dirname(__file__),'..',''))
    import atoms

class TestParticleAtom(unittest.TestCase):
    
    
    def setUp(self):
        self.atom = atoms.Atom("C")
        
        self.atom.mol = 112
        self.atom.ring  = 3
        self.atom.residue  = 387
        self.atom.resname = "THIO"
        self.atom.qgroup = 23469

    def test_checktype(self):
        self.assertEqual(self.atom.type,"atom")

    def test_properties(self):
        self.assertEqual(self.atom.tag,"C")
        
        self.assertEqual(self.atom.mol,112)
        self.assertEqual(self.atom.ring ,3)
        self.assertEqual(self.atom.residue ,387)
        self.assertEqual(self.atom.resname,"THIO")
        self.assertEqual(self.atom.qgroup,23469)
        
        # Test element setting from pymatgen
        self.assertEqual(self.atom.element.atomic_weight,12.0107)
        self.assertEqual(self.atom.element.covalent_radius,0.67)
        self.assertEqual(self.atom.element.vdw_radius,1.7)
    
    def test_str(self):
        self.assertEqual(str(self.atom)," C")

    def tearDown(self):
        del self.atom 
        
        

if __name__ == '__main__':
    unittest.main()
        