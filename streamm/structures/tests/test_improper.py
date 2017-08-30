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

try:
    import streamm.structure.impropers as impropers
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    sys.path.append(os.path.join(os.path.dirname(__file__),'..',''))
    import impropers


class TestImproper(unittest.TestCase):
    def setUp(self):
        self.imp_i = impropers.Improper(0,1,2,3)
        
    def test_str(self):
        imp_str = ' 0 - 1 - 2 - 3'
        self.assertEqual(str(self.imp_i),imp_str)

    def tearDown(self):
        del self.imp_i 
        self.imp_i = None


if __name__ == '__main__':
    unittest.main()
        
        