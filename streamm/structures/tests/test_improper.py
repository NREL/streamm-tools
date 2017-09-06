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

import streamm.structures.improper as improper

    

from streamm_testutil import *




class TestImproper(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.imp_i = improper.Improper(0,1,2,3)
        
    def test_str(self):
        imp_str = ' 0 - 1 - 2 - 3'
        self.assertEqual(str(self.imp_i),imp_str)

    @tearDown_streamm 
    def tearDown(self):
        del self.imp_i 
        self.imp_i = None


if __name__ == '__main__':

    unittest.main()    

        
                
