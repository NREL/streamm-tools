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


from streamm.calculations.nwchem import NWChem
from streamm.buildingblocks.container import Container as BBCont
import streamm.structures.particle as particle

    

class Test_Gaussian(unittest.TestCase):
    # 
    def setUp(self):
        self.calc_i = NWChem('nwchem_thiophene_SP')
        

        self.Th = BBCont('thiophene')
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

        self.Th.particles[5].rsite = 'TH'
        self.Th.particles[6].rsite = 'FH'
        self.Th.particles[8].rsite = 'TH'
        self.Th.bonded_nblist = self.Th.guess_nblist(0,radii_buffer=1.25)
        self.Th.find_rsites()
    
    def test_addstrucC(self):
        self.calc_i.add_strucC(self.Th)


    def test_writejson(self):
        self.calc_i.dump_json()
        tag_i = self.calc_i.tag
        del self.calc_i
        self.calc_i = NWChem(tag_i)
        self.calc_i.load_json()
        # Clean up files 
        os.remove(self.calc_i.files['data']['json'])


    def tearDown(self):
        del self.calc_i         
        del self.Th         



if __name__ == '__main__':
    unittest.main()        
                