# coding: utf-8
# Copyright (c) STREAMM Development Team.
# Distributed under the terms of the MIT License.


'''
Unit tests for the particles module
'''

from __future__ import division, unicode_literals

import unittest2 as unittest
import pickle
from copy import deepcopy




import unittest, os , random 
import numpy as np

from streamm import structure
from streamm import periodictable


class TestParticle(unittest.TestCase):
    def setUp(self):
	self.part = structure.Particle("genpart")
        
    def test_checktype(self):
        self.assertEqual(self.part.type,"genpart")
        
    def test_str(self):
        self.assertEqual(str(self.part)," genpart ")

    def tearDown(self):
        del self.part 
        self.part = None

class TestParticleAtom(unittest.TestCase):
    def setUp(self):
	self.part = structure.Atom(symbol="Ir")
        self.part.properties["mol"] = 1
        self.part.properties["fftype"] = "CIr"
        
    def test_checktype(self):
        self.assertEqual(self.part.type,"atom")


    def test_properties(self):
        self.assertEqual(self.part.properties["symbol"],"Ir")
        self.assertEqual(self.part.properties["number"],77)
        self.assertEqual(self.part.properties["mass"],192.217)
        self.assertEqual(self.part.properties["vdw_radii"],2.0)
        self.assertEqual(self.part.properties["cov_radii"],1.41)
        self.assertEqual(self.part.properties["mol"],1)
        self.assertEqual(self.part.properties["fftype"],"CIr")

    def test_str(self):
        self.assertEqual(str(self.part)," atom ")

    def tearDown(self):
        del self.part 
        self.part = None


if __name__ == '__main__':
    unittest.main()
