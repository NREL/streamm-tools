# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3.4"
__email__ = "organicelectronics@nrel.gov"
__status__ = "Beta"


'''
Unit tests for the particles module
'''

import logging
logger = logging.getLogger(__name__)

import unittest
import os
import numpy as np
import random
import numpy.testing.utils as nptu


import streamm.mpi.mpiBase as mpiBase 

from streamm_testutil import *



class TestMPIBase(unittest.TestCase):
    # 
    @setUp_streamm 
    def setUp(self):

        self.p = mpiBase.getMPIObject()
        # p = mpiBase.getMPIObject(True, localVerbose=True)
    
    def test_bcast(self):
        # MPI setup
        rank = self.p.getRank()
        size = self.p.getCommSize()
    
        x = 0              # Sets x on all processors
        self.p.barrier()
        if rank == 0:      # Resets x on proc 0
            x = 12.3456
        self.p.barrier()
    
        # Check x setting on each processor
        for proc in range(size):
            if proc == rank:
                self.assertEqual(x,12.3456) #, " on processor ", proc
            self.p.barrier()
        self.p.barrier()
    
        # Broadcast x to all processors (default root proc is '0')
        xAll = self.p.bcast(x)
        self.p.barrier()
    
        # Check x setting on each processor
        for proc in range(size):
            if proc == rank:
                self.assertEqual(xAll,12.3456) #, " on processor ", proc
            self.p.barrier()
        self.p.barrier()
        
    @tearDown_streamm 
    def tearDown(self):
        del self.p         

if __name__ == '__main__':

    unittest.main()    

        
                
