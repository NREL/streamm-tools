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
Unit tests for the bonds module
'''


import logging
logger = logging.getLogger(__name__)

import unittest
import os

import streamm.structures.angle as angle


from streamm_testutil import *



class TestAngle(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.angle_i = angle.Angle(2,0,1)
        self.angle_i.cosine  = 0.888754
        
    def test_str(self):
        angle_str = ' 2 - 0 - 1'
        self.assertEqual(str(self.angle_i),angle_str)
        self.assertEqual(self.angle_i.pkey1,2)
        self.assertEqual(self.angle_i.pkey2,0)
        self.assertEqual(self.angle_i.pkey3,1)
        self.assertEqual(self.angle_i.cosine,0.888754)


    def test_save(self):
        json_data = self.angle_i.export_json()
        del self.angle_i
        self.angle_i = angle.Angle(0,0,0)
        self.angle_i.import_json(json_data)
        self.assertEqual(self.angle_i.pkey1,2)
        self.assertEqual(self.angle_i.pkey2,0)
        self.assertEqual(self.angle_i.pkey3,1)
        self.assertEqual(self.angle_i.cosine,0.888754)
    
    @tearDown_streamm 
    def tearDown(self):
        del self.angle_i 
        self.angle_i = None

    
if __name__ == '__main__':

    unittest.main()    

        
                
