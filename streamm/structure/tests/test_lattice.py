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
    import streamm.structure.lattices as lattices
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    sys.path.append(os.path.join(os.path.dirname(__file__),'..',''))
    import lattices

class TestLattice(unittest.TestCase):
    def setUp(self):
        # fixture 
        matrix = [ 100,0,0,0,100,0,0,0,100 ]
        self.lat = lattices.Lattice()
        self.lat.set_matrix(matrix)
        
        self.assertEqual(self.lat._lengths[0],100.0)
        self.assertEqual(self.lat._lengths[1],100.0)
        self.assertEqual(self.lat._lengths[2],100.0)

        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        
    def test_changematrix(self):
        matrix = [ 132,0,0,0,127,0,0,0,150 ]
        self.lat.set_matrix(matrix)
        self.assertEqual(self.lat._lengths[0],132.0)
        self.assertEqual(self.lat._lengths[1],127.0)
        self.assertEqual(self.lat._lengths[2],150.0)

        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
             
    def test_LCtoLV(self):
        # Cubic 
        box = [ 12,12,12,90.0,90.0,90.0 ]
        self.lat.LCtoLV(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],12.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self._matrix[0][0],12.0)
        self.assertEqual(self._matrix[0][1],12.0)
        self.assertEqual(self._matrix[0][2],12.0)
        # 
        self.assertEqual(self._matrix[1][0],12.0)
        self.assertEqual(self._matrix[1][1],12.0)
        self.assertEqual(self._matrix[1][2],12.0)
        # 
        self.assertEqual(self._matrix[2][0],12.0)
        self.assertEqual(self._matrix[2][1],12.0)
        self.assertEqual(self._matrix[2][2],12.0)
        #
        # Tetragonal
        #
        box = [ 12,12,15,90.0,90.0,90.0 ]
        self.lat.LCtoLV(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self._matrix[0][0],12.0)
        self.assertEqual(self._matrix[0][1],12.0)
        self.assertEqual(self._matrix[0][2],12.0)
        # 
        self.assertEqual(self._matrix[1][0],12.0)
        self.assertEqual(self._matrix[1][1],12.0)
        self.assertEqual(self._matrix[1][2],12.0)
        # 
        self.assertEqual(self._matrix[2][0],12.0)
        self.assertEqual(self._matrix[2][1],12.0)
        self.assertEqual(self._matrix[2][2],12.0)
        #
        # Rhombodedral
        #
        box = [ 12,12,12,60.0,120.0,40.0 ]
        self.lat.LCtoLV(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self._matrix[0][0],12.0)
        self.assertEqual(self._matrix[0][1],12.0)
        self.assertEqual(self._matrix[0][2],12.0)
        # 
        self.assertEqual(self._matrix[1][0],12.0)
        self.assertEqual(self._matrix[1][1],12.0)
        self.assertEqual(self._matrix[1][2],12.0)
        # 
        self.assertEqual(self._matrix[2][0],12.0)
        self.assertEqual(self._matrix[2][1],12.0)
        self.assertEqual(self._matrix[2][2],12.0)        
        #
        # Orthorhombic  
        #
        box = [ 12,8,15,90.0,90.0,90.0 ]
        self.lat.LCtoLV(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self._matrix[0][0],12.0)
        self.assertEqual(self._matrix[0][1],12.0)
        self.assertEqual(self._matrix[0][2],12.0)
        # 
        self.assertEqual(self._matrix[1][0],12.0)
        self.assertEqual(self._matrix[1][1],12.0)
        self.assertEqual(self._matrix[1][2],12.0)
        # 
        self.assertEqual(self._matrix[2][0],12.0)
        self.assertEqual(self._matrix[2][1],12.0)
        self.assertEqual(self._matrix[2][2],12.0)           
        #
        # Monoclinic  
        #
        box = [ 12,8,15,90.0,60.0,90.0 ]
        self.lat.LCtoLV(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self._matrix[0][0],12.0)
        self.assertEqual(self._matrix[0][1],12.0)
        self.assertEqual(self._matrix[0][2],12.0)
        # 
        self.assertEqual(self._matrix[1][0],12.0)
        self.assertEqual(self._matrix[1][1],12.0)
        self.assertEqual(self._matrix[1][2],12.0)
        # 
        self.assertEqual(self._matrix[2][0],12.0)
        self.assertEqual(self._matrix[2][1],12.0)
        self.assertEqual(self._matrix[2][2],12.0)        
        #
        # Triclinic   
        #
        box = [ 12,8,15,60.0,120.0,80.0 ]
        self.lat.LCtoLV(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self._matrix[0][0],12.0)
        self.assertEqual(self._matrix[0][1],12.0)
        self.assertEqual(self._matrix[0][2],12.0)
        # 
        self.assertEqual(self._matrix[1][0],12.0)
        self.assertEqual(self._matrix[1][1],12.0)
        self.assertEqual(self._matrix[1][2],12.0)
        # 
        self.assertEqual(self._matrix[2][0],12.0)
        self.assertEqual(self._matrix[2][1],12.0)
        self.assertEqual(self._matrix[2][2],12.0)        
                                                                
    def tearDown(self):
        del self.lat 

if __name__ == '__main__':
    unittest.main()
        