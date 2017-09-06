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
Unit tests for the bonds module
'''


import logging
logger = logging.getLogger(__name__)

import unittest
import os

import streamm.structures.dihedral as dihedral


from streamm_testutil import *



class TestDihedral(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.dih_i = dihedral.Dihedral(2,0,1,5)
        
    def test_str(self):
        dih_str = ' 2 - 0 - 1 - 5'
        self.assertEqual(str(self.dih_i),dih_str)

    @tearDown_streamm 
    def tearDown(self):
        del self.dih_i 


    
if __name__ == '__main__':

    unittest.main()    

        
                
