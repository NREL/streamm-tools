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
import os
import numpy as np
import random
import numpy.testing.utils as nptu
import copy
import shutil


from streamm.calculations.project import Project
from streamm.calculations.gaussian import Gaussian
from streamm.calculations.nwchem import NWChem
from streamm.calculations.lammps import LAMMPS


import streamm.forcefields.particletype as particletype
import streamm.forcefields.bondtype as bondtype
import streamm.forcefields.angletype as angletype
import streamm.forcefields.dihtype as dihtype
import streamm.forcefields.imptype as imptype


from streamm.structures.buildingblock import Buildingblock 
import streamm.structures.particle as particle

from streamm.calculations.resource import Resource

import streamm.util.replicate as replicate

from streamm_testutil import *

class TestReplicate(unittest.TestCase):
    # 
    @setUp_streamm
    def setUp(self):

        self.Th = Buildingblock('thiophene')
        symbols = ['C','C','C','C','S','H','H','H','H']
        positions = [ ]
        positions.append([-1.55498576,-1.91131218,-0.00081000])
        positions.append([-0.17775976,-1.91131218,-0.00081000])
        positions.append([0.34761524,-0.57904218,-0.00081000])
        positions.append([-0.65884476,0.36101082,0.00000000])
        positions.append([-2.16948076,-0.35614618,-0.00000800])
        positions.append([-2.18966076,-2.79526518,-0.00132100])
        positions.append([0.45389024,-2.80145418,-0.00106400])
        positions.append([1.41682424,-0.35961818,-0.00138200])
        positions.append([-0.51943676,1.44024682,0.00064700])
        for i in range(len(symbols)):
            pt_i = particle.Particle(symbol=symbols[i])
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.Th.add_partpos(pt_i,pos_i)

        for pkey,p in self.Th.particles.iteritems():
            print  pkey,p
            if( p.symbol == 'C' ):
                p.ffkey = 'CA'
            elif( p.symbol == 'H' ):
                p.ffkey = 'HA'
            elif( p.symbol == 'S' ):
                p.ffkey = 'S'
                
        self.Th.particles[5].rsite = 'TH'
        self.Th.particles[6].rsite = 'FH'
        self.Th.particles[8].rsite = 'TH'
        self.Th.bonded_nblist = self.Th.guess_nblist(0,radii_buffer=1.25)
        self.Th.find_rsites()
        
        
    def test_add_struct(self):
        

        seed = 82343
        blank_strucC =  Buildingblock()
        self.strucC = replicate.add_struc(blank_strucC,self.Th,10,seed)
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            mol_i = int(pkey_i/self.Th.n_particles) 
            self.assertEqual(str(particle_i.mol),str(mol_i))

        self.assertEqual(self.strucC.n_molecules(),9)
            
        file_i = os.path.join(TEST_DIR, "th_x10.xyz")
        self.strucC.write_xyz(file_i)
            
        
        
        
    def test_add_struct_on_grid(self):


        blank_strucC =  Buildingblock()
        self.strucC = replicate.add_struc_grid(blank_strucC,self.Th,10)
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            mol_i = int(pkey_i/self.Th.n_particles) 
            self.assertEqual(str(particle_i.mol),str(mol_i))

        self.assertEqual(self.strucC.n_molecules(),9)
            
        file_i = os.path.join(TEST_DIR, "th_x10.xyz")
        self.strucC.write_xyz(file_i)
            
        
    @tearDown_streamm 
    def tearDown(self):
        del self.Th
        
if __name__ == '__main__':
    unittest.main()    
        
                   
        