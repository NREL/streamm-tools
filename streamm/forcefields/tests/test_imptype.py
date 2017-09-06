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

import streamm.forcefields.imptype as imptype


HOME_DIR = os.getcwd()
RELATIVE_TEST_DIR = os.path.join(os.path.dirname(__file__))
TEST_DIR = os.path.join(HOME_DIR,RELATIVE_TEST_DIR)

class Testimproper(unittest.TestCase):
    
    def setUp(self):
        self.imptype_i = imptype.Imptype("C1","C2","C3","C4",type="improper")
        self.imptype_i.e0 = 180.0
        self.imptype_i.ke = 67.3
        self.imptype_i.pn = 4.0
        

    def test_impstr(self):
        imp_str = ' improper  C1 - C2 - C3 - C4 type improper \n  imp e0 = 180.000000 ke = 67.300000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.imptype_i),imp_str)
        
    def tearDown(self):
        del self.imptype_i 
        self.imptype_i = None


if __name__ == '__main__':
    os.chdir(TEST_DIR)
    unittest.main()    
    os.chdir(HOME_DIR)
        
                
