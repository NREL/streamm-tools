# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Ph.D."
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3.4"
__email__ = "organicelectronics@nrel.gov"
__status__ = "Beta"


'''
Unit tests for the particles module
'''

import logging
logger = logging.getLogger(__name__)

import unittest
import os

import streamm.forcefields.particletype as particletype

from streamm_testutil import *



    
class Test_FFparticle(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.ff = particletype.Particletype("CC")
        self.ff.epsilon = 1.05
        self.ff.sigma = 3.25
        
        

    def test_str(self):        
        self.assertEqual(str(self.ff) ,' CC epsilon:1.05 sigma:3.25')
        
    def test_properties(self):
        self.assertEqual(self.ff.epsilon,1.05)
        self.assertEqual(self.ff.sigma,3.25)
        

    def test_save(self):
        json_data = self.ff.export_json()
        del self.ff
        self.ff = particletype.Particletype("X")
        self.ff.import_json(json_data)
        self.assertEqual(self.ff.epsilon,1.05)
        self.assertEqual(self.ff.sigma,3.25)
                
    @tearDown_streamm 
    def tearDown(self):
        del self.ff 


if __name__ == '__main__':

    unittest.main()    

        
                


