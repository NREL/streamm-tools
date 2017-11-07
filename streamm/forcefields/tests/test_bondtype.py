# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"


'''
Unit tests for the particles module
'''

import logging
logger = logging.getLogger(__name__)

import unittest
import os

import streamm.forcefields.bondtype as bondtype
            

from streamm_testutil import *


    
class Testbondtype(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.bondtype_i = bondtype.Bondtype("C","H")
        self.bondtype_i.r0 = 0.56
        self.bondtype_i.kb = 24.023

    def test_bondstr(self):
        bond_str = ' bond  C - H type harmonic \n  harmonic r_0 = 0.560000 K = 24.023000 lammps index 0  gromacs index 0  '
        self.assertEqual(str(self.bondtype_i),bond_str)
        
    def test_setharmonic(self):
        self.bondtype_i.setharmonic(0.76,33.33)
        bond_str = ' bond  C - H type harmonic \n  harmonic r_0 = 0.760000 K = 33.330000 lammps index 0  gromacs index 0  '
        self.assertEqual(str(self.bondtype_i),bond_str)

    def test_save(self):
        json_data = self.bondtype_i.export_json()
        del self.bondtype_i
        self.bondtype_i = bondtype.Bondtype("X","X")
        self.bondtype_i.import_json(json_data)
        
        self.assertEqual(self.bondtype_i.fftype1,'C')
        self.assertEqual(self.bondtype_i.fftype2,'H')
        self.assertEqual(self.bondtype_i.type,'harmonic')
        self.assertEqual(self.bondtype_i.r0,0.56)
        self.assertEqual(self.bondtype_i.kb,24.023)
        
        
    @tearDown_streamm 
    def tearDown(self):
        del self.bondtype_i 

if __name__ == '__main__':

    unittest.main()    

        
