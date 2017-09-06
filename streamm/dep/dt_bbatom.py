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


try:
    import streamm.buildingblock.bbatom as bbatom

except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    import bbatom

class TestBBatom(unittest.TestCase):
    def setUp(self):
        self.particle_i = bbatom.BBatom("H")
        self.particle_i.bbid = 'R'

    def test_name(self):
        self.assertEqual(str(self.particle_i)," H (R)")

    def tearDown(self):
        del self.particle_i 
        self.particle_i = None


if __name__ == '__main__':
    unittest.main()
        