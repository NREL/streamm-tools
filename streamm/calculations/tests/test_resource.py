
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
import copy 
import math
import sys
import os
import shutil

from streamm_testutil import *



import streamm.calculations.resource as resource 
    
PRECISION = 4

class TestResource(unittest.TestCase):
    
    @setUp_streamm 
    def setUp(self):
        
        self.res_tag = 'my_test_res'  # Change this to remote to run the calculations remotely 
        self.res_i = resource.Resource(self.res_tag )
        
    def test_make_dir(self):
        self.res_i.make_dir()

    def test_json(self):
        json_data = self.res_i.export_json()
        home_dir = self.res_i.dir['home']
        del self.res_i
        self.res_i = resource.Resource(self.res_tag)
        self.res_i.import_json(json_data)
        
        self.assertEqual(self.res_i.meta['type'],'local')
        self.assertEqual(home_dir,self.res_i.dir['home'])
        
        
    @tearDown_streamm 
    def tearDown(self):
        del self.res_i
        
        

if __name__ == '__main__':

    unittest.main()    

        
