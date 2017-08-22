# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0


'''
Unit tests for the particles module
'''

from __future__ import division, unicode_literals

import logging
logger = logging.getLogger(__name__)

import unittest

try:
    import streamm.structure.particles as particles
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    import particles
    
    
class TestParticle(unittest.TestCase):
    def setUp(self):
        self.part = particles.Particle("genpart")
        
    def test_checktype(self):
        self.assertEqual(self.part.type,"genpart")
        
    def test_str(self):
        self.assertEqual(str(self.part)," genpart")

    def tearDown(self):
        del self.part 
        self.part = None



if __name__ == '__main__':
    unittest.main()
