# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0


'''
Unit tests for the Lattice module
'''

from __future__ import division, unicode_literals

import logging
logger = logging.getLogger(__name__)

import unittest
import numpy as np
import numpy.testing.utils as nptu

try:
    import streamm.structure.lattices as lattices
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
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
                

    def test_set_box(self):
        # Cubic 
        box = [ 12,12,12,90.0,90.0,90.0 ]
        self.lat.set_box(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],12.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self.lat._matrix[0][0],12.0)
        self.assertEqual(self.lat._matrix[0][1],0.0)
        self.assertEqual(self.lat._matrix[0][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[1][0],0.0)
        self.assertEqual(self.lat._matrix[1][1],12.0)
        self.assertEqual(self.lat._matrix[1][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[2][0],0.0)
        self.assertEqual(self.lat._matrix[2][1],0.0)
        self.assertEqual(self.lat._matrix[2][2],12.0)
        #
        # Tetragonal
        #
        box = [ 12,12,15,90.0,90.0,90.0 ]
        self.lat.set_box(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self.lat._matrix[0][0],12.0)
        self.assertEqual(self.lat._matrix[0][1],0.0)
        self.assertEqual(self.lat._matrix[0][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[1][0],0.0)
        self.assertEqual(self.lat._matrix[1][1],12.0)
        self.assertEqual(self.lat._matrix[1][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[2][0],0.0)
        self.assertEqual(self.lat._matrix[2][1],0.0)
        self.assertEqual(self.lat._matrix[2][2],15.0)
        #
        # Rhombodedral
        #
        box = [ 12,12,12,60.0,120.0,40.0 ]
        self.lat.set_box(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],12.0)
        self.assertEqual(self.lat._lengths[2],12.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],60.0)
        self.assertEqual(self.lat._angles[1],120.0)
        self.assertEqual(self.lat._angles[2],40.0)
        # Check matrix
        self.assertEqual(self.lat._matrix[0][0],12.0)
        self.assertEqual(self.lat._matrix[0][1],0.0)
        self.assertEqual(self.lat._matrix[0][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[1][0],9.19253332)
        self.assertEqual(self.lat._matrix[1][1],7.71345132)
        self.assertEqual(self.lat._matrix[1][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[2][0],-6.0)
        self.assertEqual(self.lat._matrix[2][1],-2.3092942)
        self.assertEqual(self.lat._matrix[2][2],10.13248046)        
        #
        # Orthorhombic  
        #
        box = [ 12,8,15,90.0,90.0,90.0 ]
        self.lat.set_box(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],8.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self.lat._matrix[0][0],12.0)
        self.assertEqual(self.lat._matrix[0][1],0.0)
        self.assertEqual(self.lat._matrix[0][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[1][0],0.0)
        self.assertEqual(self.lat._matrix[1][1],8.0)
        self.assertEqual(self.lat._matrix[1][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[2][0],0.0)
        self.assertEqual(self.lat._matrix[2][1],0.0)
        self.assertEqual(self.lat._matrix[2][2],15.0)           
        #
        # Monoclinic  
        #
        box = [ 12,8,15,90.0,60.0,90.0 ]
        self.lat.set_box(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],8.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],60.0)
        self.assertEqual(self.lat._angles[2],90.0)
        # Check matrix
        self.assertEqual(self.lat._matrix[0][0],12.0)
        self.assertEqual(self.lat._matrix[0][1],0.0)
        self.assertEqual(self.lat._matrix[0][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[1][0],0.0)
        self.assertEqual(self.lat._matrix[1][1],8.0)
        self.assertEqual(self.lat._matrix[1][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[2][0],7.5)
        self.assertEqual(self.lat._matrix[2][1],0.0)
        self.assertEqual(self.lat._matrix[2][2],12.99038106)        
        #
        # Triclinic   
        #
        box = [ 12,8,15,60.0,120.0,80.0 ]
        self.lat.set_box(box)
        # Check length
        self.assertEqual(self.lat._lengths[0],12.0)
        self.assertEqual(self.lat._lengths[1],8.0)
        self.assertEqual(self.lat._lengths[2],15.0)
        # Check angles 
        self.assertEqual(self.lat._angles[0],60.0)
        self.assertEqual(self.lat._angles[1],120.0)
        self.assertEqual(self.lat._angles[2],80.0)
        # Check matrix
        self.assertEqual(self.lat._matrix[0][0],12.0)
        self.assertEqual(self.lat._matrix[0][1],0.0)
        self.assertEqual(self.lat._matrix[0][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[1][0], 1.38918542)
        self.assertEqual(self.lat._matrix[1][1],7.87846202)
        self.assertEqual(self.lat._matrix[1][2],0.0)
        # 
        self.assertEqual(self.lat._matrix[2][0],-7.5)
        self.assertEqual(self.lat._matrix[2][1],5.4622685)
        self.assertEqual(self.lat._matrix[2][2],11.78616234)        
    
    def test_deltasq_pos(self):
        pos_i  = [123.9,1298.0,93.762]
        pos_j = [832.123,112.127398,9374.9123]
        
        dr_ij_correct = np.array([708.223   , -1185.872602,  9281.1503 ])
        dr_ij = self.lat.deltasq_pos(pos_i,pos_j)
        
        nptu.assert_almost_equal(dr_ij,dr_ij_correct)
        
    def tearDown(self):
        del self.lat 
        self.lat = None

if __name__ == '__main__':
    unittest.main()
        
        