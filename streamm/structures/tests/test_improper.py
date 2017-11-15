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

import streamm.structures.improper as improper

    

from streamm_testutil import *

class TestImproper(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.imp_i = improper.Improper(0,1,2,3)
        self.imp_i.cosine  = -0.4875
        
    def test_str(self):
        imp_str = ' 0 - 1 - 2 - 3'
        self.assertEqual(str(self.imp_i),imp_str)

    def test_save(self):
        json_data = self.imp_i.export_json()
        del self.imp_i
        self.imp_i = improper.Improper(0,0,0,0)
        self.imp_i.import_json(json_data)
        self.assertEqual(self.imp_i.pkey1,0)
        self.assertEqual(self.imp_i.pkey2,1)
        self.assertEqual(self.imp_i.pkey3,2)
        self.assertEqual(self.imp_i.pkey4,3)
        self.assertEqual(self.imp_i.cosine,-0.4875)
        
    @tearDown_streamm 
    def tearDown(self):
        del self.imp_i 
        self.imp_i = None


if __name__ == '__main__':

    unittest.main()    

        
                
