# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals


'''
Unit tests for the particles module
'''

import logging
logger = logging.getLogger(__name__)

import unittest
import os

import streamm.structures.nblist as nblist

from streamm_testutil import *



    
class Test_NBlist(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.nblist = nblist.NBlist()
        
        # Structure of vanillin 
        self.particles = {}
        self.particles[0] = {'symbol':'C'}
        self.particles[1] = {'symbol':'C'}
        self.particles[2] = {'symbol':'C'}
        self.particles[3] = {'symbol':'C'}
        self.particles[4] = {'symbol':'C'}
        self.particles[5] = {'symbol':'C'}
        self.particles[6] = {'symbol':'C'}
        self.particles[7] = {'symbol':'H'}
        self.particles[8] = {'symbol':'H'}
        self.particles[9] = {'symbol':'H'}
        self.particles[10] = {'symbol':'H'}
        self.particles[11] = {'symbol':'O'}
        self.particles[12] = {'symbol':'O'}
        self.particles[13] = {'symbol':'C'}
        self.particles[14] = {'symbol':'H'}
        self.particles[15] = {'symbol':'H'}
        self.particles[16] = {'symbol':'H'}
        self.particles[17] = {'symbol':'O'}
        self.particles[18] = {'symbol':'H'}
        
        self.nblist.list = [1, 5, 6, 0, 2, 9, 1, 3, 12, 4, 2, 17, 3, 5, 8, 4, 0, 7, 10, 11, 0, 5, 4, 1, 6, 6, 2, 13, 12, 14, 15, 16, 13, 13, 13, 3, 18, 17]
        self.nblist.index = [0, 3, 6, 9, 12, 15, 18, 21, 22, 23, 24, 25, 26, 28, 32, 33, 34, 35, 37, 38]
        
    def test_str(self):
        self.assertEqual(str(self.nblist)," NBlist of 20 particle with 38 connections")
                
    def test_getnbs(self):
        self.sub_check = []
        self.sub_check.append( [1, 5, 6] )
        self.sub_check.append( [0, 2, 9] )
        self.sub_check.append( [1, 3, 12] )
        self.sub_check.append( [4, 2, 17] )
        self.sub_check.append( [3, 5, 8] )
        self.sub_check.append( [4, 0, 7] )
        self.sub_check.append( [10, 11, 0] )
        self.sub_check.append( [5] )
        self.sub_check.append( [4] )
        self.sub_check.append( [1] )
        self.sub_check.append( [6] )
        self.sub_check.append( [6] )
        self.sub_check.append( [2, 13] )
        self.sub_check.append( [12, 14, 15, 16] )
        self.sub_check.append( [13] )
        self.sub_check.append( [13] )
        self.sub_check.append( [13] )
        self.sub_check.append( [3, 18] )
        self.sub_check.append( [17] )

        for pkey,particle_i  in self.particles.iteritems():
            nblist_i = self.nblist.getnbs(pkey)
            self.assertListEqual(nblist_i,self.sub_check[pkey])

    def test_calc_nnab(self):
        self.nnab_check = []
        self.nnab_check.append( 3 )
        self.nnab_check.append( 3 )
        self.nnab_check.append( 3 )
        self.nnab_check.append( 3 )
        self.nnab_check.append( 3 )
        self.nnab_check.append( 3 )
        self.nnab_check.append( 3 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 2 )
        self.nnab_check.append( 4 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 1 )
        self.nnab_check.append( 2 )
        self.nnab_check.append( 1 )

        for pkey,particle_i  in self.particles.iteritems():
            nblist_i = self.nblist.getnbs(pkey)
            nnab = self.nblist.calc_nnab(pkey)
            self.assertEqual(nnab,self.nnab_check[pkey])
            self.assertEqual(len(nblist_i),self.nnab_check[pkey])
            
            
    @tearDown_streamm 
    def tearDown(self):
        del self.nblist 

    
if __name__ == '__main__':

    unittest.main()    

        
                
