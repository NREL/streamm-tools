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
Unit tests for the gaussian module
'''

import logging
logger = logging.getLogger(__name__)

import unittest
import numpy as np
import random
import numpy.testing.utils as nptu

from streamm.calculations.gaussian import Gaussian
from streamm.structures.buildingblock import Buildingblock 
import streamm.structures.particle as particle

from streamm.calculations.resource import Resource

from streamm_testutil import * 

class Test_Gaussian(unittest.TestCase):
    #
    @setUp_streamm
    def setUp(self):
        
        self.calc_i = Gaussian('gaussian_thiophene_SP')
        

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

        self.Th.particles[5].rsite = 'TH'
        self.Th.particles[6].rsite = 'FH'
        self.Th.particles[8].rsite = 'TH'
        self.Th.bonded_nblist = self.Th.guess_nblist(0,radii_buffer=1.25)
        self.Th.find_rsites()
        self.calc_i.add_strucC(self.Th)
    

        self.res_tag = 'local'  # Change this to remote to run the calculations remotely 
        self.res_i = Resource(self.res_tag )
        self.res_i.dir['templates'] = TEMPLATE_PATH
        self.res_i.make_dir()
        self.res_i.export_json()

        self.calc_i.set_resource(self.res_i)
        
        self.calc_i.make_dir()
        

    def test_writeinput(self):
        file_type = 'templates'
        file_key = 'run'
        file_name = "gaussian.sh"
        from_dirkey = 'templates'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
        file_type = 'templates'
        file_key = 'com'
        file_name = "gaussian.com"
        from_dirkey = 'templates'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
        os.chdir(self.calc_i.dir['scratch'])
        
        self.calc_i.load_str('templates','com')
        
        self.calc_i.load_str('templates','run')
                
        self.calc_i.properties['commands'] = 'HF/3-21G SP'
        self.calc_i.properties['charge'] = 0
        self.calc_i.properties['spin_mult'] = 1
        self.calc_i.properties['coord'] = self.calc_i.strucC.write_coord()
        
        self.calc_i.replacewrite_prop('com','input','com','%s.com'%(self.calc_i.tag))
        
        self.calc_i.properties['input_com'] = self.calc_i.files['input']['com']
        self.calc_i.replacewrite_prop('run','scripts','run','%s.sh'%(self.calc_i.tag))

        os.chdir(self.calc_i.dir['home'])
        self.calc_i.export_json()
        
        os.chdir(self.calc_i.dir['scratch'])
        self.calc_i.check()
                
        self.calc_i.strucC.calc_charge()
        
        self.assertEqual(self.calc_i.strucC.charge,0.0)
        
        os.chdir(self.calc_i.dir['home'])


    def test_writejson(self):
        self.calc_i.export_json()
        tag_i = self.calc_i.tag
        del self.calc_i
        self.calc_i = Gaussian(tag_i)
        self.calc_i.import_json()
        
    @tearDown_streamm
    def tearDown(self):
        del self.calc_i         
        del self.Th
        

if __name__ == '__main__':
    unittest.main()
    
        
                
                
