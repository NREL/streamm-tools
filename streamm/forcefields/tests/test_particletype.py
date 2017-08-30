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
    import streamm.forcefield.particletype as particletype
            
    
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    import particletype
    
    
class Test_FFparticle(unittest.TestCase):
    def setUp(self):
        self.ff = particletype.Particletype("CC")
        self.ff.epsilon = 1.05
        self.ff.sigma = 3.25
        
        

    def test_str(self):        
        self.assertEqual(str(self.ff) ,' CC epsilon:1.05 sigma:3.25')
        
    def test_properties(self):
        self.assertEqual(self.ff.epsilon,1.05)
        self.assertEqual(self.ff.sigma,3.25)
        
    def tearDown(self):
        del self.ff 


if __name__ == '__main__':
    unittest.main()
        
                

