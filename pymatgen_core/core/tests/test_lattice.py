# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import itertools
from pymatgen_core.core.lattice import Lattice
import numpy as np


import unittest #2 as unittest
import numpy.testing.utils as nptu

# from pymatgen.core.operations import SymmOp



class TestLattice(unittest.TestCase):
    def setUp(self):
        # fixture 
        matrix = [ 100,0,0,0,100,0,0,0,100 ]
        self.lat = Lattice(matrix)
        
        self.assertEqual(self.lat.lengths[0],100.0)
        self.assertEqual(self.lat.lengths[1],100.0)
        self.assertEqual(self.lat.lengths[2],100.0)

        self.assertEqual(self.lat.angles[0],90.0)
        self.assertEqual(self.lat.angles[1],90.0)
        self.assertEqual(self.lat.angles[2],90.0)
        
    def test_changematrix(self):
        self.lat.matrix = [ 132,0,0,0,127,0,0,0,150 ]

        self.assertEqual(self.lat.lengths[0],132.0)
        self.assertEqual(self.lat.lengths[1],127.0)
        self.assertEqual(self.lat.lengths[2],150.0)

        self.assertEqual(self.lat.angles[0],90.0)
        self.assertEqual(self.lat.angles[1],90.0)
        self.assertEqual(self.lat.angles[2],90.0)
        

    def tearDown(self):
        del self.lat 
        self.lat = None
        

class TestLattice(unittest.TestCase):
    def setUp(self):
        # fixture 
        matrix = [ 100,0,0,0,100,0,0,0,100 ]
        self.lat = Lattice()
        self.lat.matrix = matrix
        
        self.assertAlmostEqual(self.lat.lengths[0],100.0)
        self.assertAlmostEqual(self.lat.lengths[1],100.0)
        self.assertAlmostEqual(self.lat.lengths[2],100.0)

        self.assertAlmostEqual(self.lat.angles[0],90.0)
        self.assertAlmostEqual(self.lat.angles[1],90.0)
        self.assertAlmostEqual(self.lat.angles[2],90.0)
        
    def test_changematrix(self):
        self.lat.matrix =  [ 132,0,0,0,127,0,0,0,150 ]
        self.assertAlmostEqual(self.lat.lengths[0],132.0)
        self.assertAlmostEqual(self.lat.lengths[1],127.0)
        self.assertAlmostEqual(self.lat.lengths[2],150.0)

        self.assertAlmostEqual(self.lat.angles[0],90.0)
        self.assertAlmostEqual(self.lat.angles[1],90.0)
        self.assertAlmostEqual(self.lat.angles[2],90.0)
                

    def test_set_box(self):
        # Cubic 
        box = [ 12,12,12,90.0,90.0,90.0 ]
        self.lat.constants = box
        # Check length
        self.assertAlmostEqual(self.lat.lengths[0],12.0)
        self.assertAlmostEqual(self.lat.lengths[1],12.0)
        self.assertAlmostEqual(self.lat.lengths[2],12.0)
        # Check angles 
        self.assertAlmostEqual(self.lat.angles[0],90.0)
        self.assertAlmostEqual(self.lat.angles[1],90.0)
        self.assertAlmostEqual(self.lat.angles[2],90.0)
        # Check matrix
        self.assertAlmostEqual(self.lat.matrix[0][0],12.0)
        self.assertAlmostEqual(self.lat.matrix[0][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[0][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[1][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[1][1],12.0)
        self.assertAlmostEqual(self.lat.matrix[1][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[2][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][2],12.0)
        #
        # Tetragonal
        #
        self.lat.constants = [ 12,12,15,90.0,90.0,90.0 ]
        # Check length
        self.assertAlmostEqual(self.lat.lengths[0],12.0)
        self.assertAlmostEqual(self.lat.lengths[1],12.0)
        self.assertAlmostEqual(self.lat.lengths[2],15.0)
        # Check angles 
        self.assertAlmostEqual(self.lat.angles[0],90.0)
        self.assertAlmostEqual(self.lat.angles[1],90.0)
        self.assertAlmostEqual(self.lat.angles[2],90.0)
        # Check matrix
        self.assertAlmostEqual(self.lat.matrix[0][0],12.0)
        self.assertAlmostEqual(self.lat.matrix[0][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[0][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[1][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[1][1],12.0)
        self.assertAlmostEqual(self.lat.matrix[1][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[2][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][2],15.0)
        #
        # Rhombodedral
        #
        self.lat.constants = [ 12,12,12,60.0,120.0,40.0 ]
        # Check length
        self.assertAlmostEqual(self.lat.lengths[0],12.0)
        self.assertAlmostEqual(self.lat.lengths[1],12.0)
        self.assertAlmostEqual(self.lat.lengths[2],12.0)
        # Check angles 
        self.assertAlmostEqual(self.lat.angles[0],60.0)
        self.assertAlmostEqual(self.lat.angles[1],120.0)
        self.assertAlmostEqual(self.lat.angles[2],40.0)
        # Check matrix
        self.assertAlmostEqual(self.lat.matrix[0][0],12.0)
        self.assertAlmostEqual(self.lat.matrix[0][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[0][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[1][0],9.19253332)
        self.assertAlmostEqual(self.lat.matrix[1][1],7.71345132)
        self.assertAlmostEqual(self.lat.matrix[1][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[2][0],-6.0)
        self.assertAlmostEqual(self.lat.matrix[2][1],-2.3092942)
        self.assertAlmostEqual(self.lat.matrix[2][2],10.13248046)        
        #
        # Orthorhombic  
        #
        self.lat.constants = [ 12,8,15,90.0,90.0,90.0 ]
        # Check length
        self.assertAlmostEqual(self.lat.lengths[0],12.0)
        self.assertAlmostEqual(self.lat.lengths[1],8.0)
        self.assertAlmostEqual(self.lat.lengths[2],15.0)
        # Check angles 
        self.assertAlmostEqual(self.lat.angles[0],90.0)
        self.assertAlmostEqual(self.lat.angles[1],90.0)
        self.assertAlmostEqual(self.lat.angles[2],90.0)
        # Check matrix
        self.assertAlmostEqual(self.lat.matrix[0][0],12.0)
        self.assertAlmostEqual(self.lat.matrix[0][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[0][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[1][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[1][1],8.0)
        self.assertAlmostEqual(self.lat.matrix[1][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[2][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][2],15.0)           
        #
        # Monoclinic  
        #
        self.lat.constants = [ 12,8,15,90.0,60.0,90.0 ]
        # Check length
        self.assertAlmostEqual(self.lat.lengths[0],12.0)
        self.assertAlmostEqual(self.lat.lengths[1],8.0)
        self.assertAlmostEqual(self.lat.lengths[2],15.0)
        # Check angles 
        self.assertAlmostEqual(self.lat.angles[0],90.0)
        self.assertAlmostEqual(self.lat.angles[1],60.0)
        self.assertAlmostEqual(self.lat.angles[2],90.0)
        # Check matrix
        self.assertAlmostEqual(self.lat.matrix[0][0],12.0)
        self.assertAlmostEqual(self.lat.matrix[0][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[0][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[1][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[1][1],8.0)
        self.assertAlmostEqual(self.lat.matrix[1][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[2][0],7.5)
        self.assertAlmostEqual(self.lat.matrix[2][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][2],12.99038106)        
        #
        # Triclinic   
        #
        self.lat.constants = [ 12,8,15,60.0,120.0,80.0 ]
        # Check length
        self.assertAlmostEqual(self.lat.lengths[0],12.0)
        self.assertAlmostEqual(self.lat.lengths[1],8.0)
        self.assertAlmostEqual(self.lat.lengths[2],15.0)
        # Check angles 
        self.assertAlmostEqual(self.lat.angles[0],60.0)
        self.assertAlmostEqual(self.lat.angles[1],120.0)
        self.assertAlmostEqual(self.lat.angles[2],80.0)
        # Check matrix
        self.assertAlmostEqual(self.lat.matrix[0][0],12.0)
        self.assertAlmostEqual(self.lat.matrix[0][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[0][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[1][0], 1.38918542)
        self.assertAlmostEqual(self.lat.matrix[1][1],7.87846202)
        self.assertAlmostEqual(self.lat.matrix[1][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[2][0],-7.5)
        self.assertAlmostEqual(self.lat.matrix[2][1],5.4622685)
        self.assertAlmostEqual(self.lat.matrix[2][2],11.78616234)        
    
    def test_deltasq_pos(self):
        pos_i  = [123.9,1298.0,93.762]
        pos_j = [832.123,112.127398,9374.9123]
        
        dr_ij_correct = np.array([708.223   , -1185.872602,  9281.1503 ])
        dr_ij = self.lat.deltasq_pos(pos_i,pos_j)
        
        nptu.assert_almost_equal(dr_ij,dr_ij_correct)
        
    def test_deltasq_pos_c(self):
        pos_i  = [123.9,1298.0,93.762]
        pos_j = [832.123,112.127398,9374.9123]
        
        dr_ij_correct = np.array([708.223   , -1185.872602,  9281.1503 ])
        dr_ij = self.lat.deltasq_pos(pos_i,pos_j)
        
        nptu.assert_almost_equal(dr_ij,dr_ij_correct)
        
    def test_norm_delta_pos_c(self):
        pos_i  = [123.9,1298.0,93.762]
        pos_j = [832.123,112.127398,9374.9123]
        
        dr_ij_correct = np.array([0.3295766,  0.5662239, -0.7554931])
        dr_ij = self.lat.norm_delta_pos_c(pos_i,pos_j)
        
        nptu.assert_almost_equal(dr_ij,dr_ij_correct)
        
        
        
    def test_delta_pos_c(self):
        pos_i  = [123.9,1298.0,93.762]
        pos_j = [832.123,112.127398,9374.9123]
        
        dr_ij_correct = np.array([708.223   , -1185.872602,  9281.1503 ])
        mag_dr_ij_correct  = 9383.3695726584974
        dr_ij,mag_dr_ij = self.lat.delta_pos(pos_i,pos_j)
        
        nptu.assert_almost_equal(dr_ij,dr_ij_correct)
        self.assertAlmostEqual(mag_dr_ij,mag_dr_ij_correct)
        
    def test_delta_npos(self):
        
        npos_i = []
        npos_i.append([-123.9,1298.0,93.762])
        npos_i.append([23487.885,-364,702673.0121])
        npos_i.append([  32482.299, 77.917240,-12378.88234])   
        
        npos_j = []
        npos_j.append([13.9,-23487.12,-3289.12834])
        npos_j.append([918034.1234,-12648.0,-5.5])
        npos_j.append([  487.1234, 7959236.314,-12378.883])
        
        
        nd_ij_correct = []
        nd_ij_correct.append([ 25015.2975059,   918263.9362965,  7957948.1117718])
        nd_ij_correct.append([706730.6774047,  1137594.3743965,  7991687.2690533])
        nd_ij_correct.append([  41135.5265133,   885729.6897238,  7959222.7055145])
        
        # nd_ij_correct = np.array([ 708.223   , -1185.872602,  9281.1503 ])
        npos_ij,nd_ij = self.lat.delta_npos(npos_i,npos_j)
        
        # nptu.assert_almost_equal(npos_ij,npos_ij_correct)
        nptu.assert_almost_equal(nd_ij,nd_ij_correct)
        
    def test_proximitycheck(self):
        
        npos_i = []
        npos_i.append([-123.9,1298.0,93.762])
        npos_i.append([13.9,-23487.12,-3289.12834])
        npos_i.append([23487.885,-364,702673.0121])
        
        npos_j = []
        npos_j.append([  32482.299, 77.917240,-12378.88234])   
        npos_j.append([6.234,-23467.12,-3299.12834])
        npos_j.append([  487.1234, 7959236.314,-12378.883])
            
        pos_cut = 10.0 
        overlap = self.lat.proximitycheck(npos_i,npos_j,pos_cut)
        self.assertTrue(overlap)
        
        pos_cut = 100.0 
        overlap = self.lat.proximitycheck(npos_i,npos_j,pos_cut)
        self.assertFalse(overlap)
        
    def test_expandmatrix(self):
        
        matrix = [ 100,0,0,0,100,0,0,0,100 ]
        self.lat.setmatrix(matrix)
        
        exlat_frac = 0.0233
        self.lat.expandmatrix(exlat_frac)
        #
        # Check matrix
        # 
        self.assertAlmostEqual(self.lat.matrix[0][0],102.33)
        self.assertAlmostEqual(self.lat.matrix[0][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[0][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[1][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[1][1],102.33)
        self.assertAlmostEqual(self.lat.matrix[1][2],0.0)
        # 
        self.assertAlmostEqual(self.lat.matrix[2][0],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][1],0.0)
        self.assertAlmostEqual(self.lat.matrix[2][2],102.33)
        
    def test_random_pos(self):
        
        pos_o = self.lat.random_pos()
        

    def test_fractional(self):
        
        matrix_i = self.lat.matrix
        matrix_i[0][0] = 50.0 
        matrix_i[0][1] = 50.0 
        matrix_i[1][0] = 50.0 
        matrix_i[1][1] = -50.0 
        matrix_i[2][0] = 10.0 
        matrix_i[2][1] = 0.0 
        matrix_i[2][2] = 60.0 
        self.lat.matrix = matrix_i
        
        frac_o = np.array([0.5,0.5,0.5])
        pos_o = self.lat.fractoreal(frac_o)
        nptu.assert_almost_equal(pos_o,[ 55. ,  0. , 30.])
                                 
        frac_o = np.array([0.5,0.0,0.5])
        pos_o = self.lat.fractoreal(frac_o)
        nptu.assert_almost_equal(pos_o,[ 30. , 25. , 30.])
                                 
        frac_o = np.array([0.0,-0.50,0.5])
        pos_o = self.lat.fractoreal(frac_o)
        nptu.assert_almost_equal(pos_o,[-20. , 25. , 30.])
        
    def test_lat_cubic(self):
        len_o = 256.16
        matrix_correct = []
        matrix_correct.append([len_o,0.0,0.0])
        matrix_correct.append([0.0,len_o,0.0])
        matrix_correct.append([0.0,0.0,len_o])
        self.lat.set_cubic(len_o)
        nptu.assert_almost_equal(self.lat.matrix,matrix_correct)
        
    def tearDown(self):
        del self.lat 
        self.lat = None
        
if __name__ == '__main__':
    import unittest
    unittest.main()
