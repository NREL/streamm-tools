
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
import copy 
import math
import sys
import os
import shutil

import streamm.calculations.calculation as calculation
import streamm.calculations.resource as resource 
from streamm.structures.buildingblock import Buildingblock 
import streamm.structures.particle as particle
    


from streamm_testutil import *



PRECISION = 4

class TestResource(unittest.TestCase):
    
    @setUp_streamm 
    def setUp(self):
        
        self.res_tag = 'local'  # Change this to remote to run the calculations remotely 
        self.res_i = resource.Resource(self.res_tag )
        self.res_i.dump_json()
        
    def test_make_dir(self):
        self.res_i.make_dir()

    def test_json(self):
        self.res_i.dump_json()
        self.res_i.load_json()
        

    @tearDown_streamm 
    def tearDown(self):
        del self.res_i
        
        
class TestAddStruc(unittest.TestCase):
    
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

        self.Th.particles[5].rsite = 'TH'
        self.Th.particles[6].rsite = 'FH'
        self.Th.particles[8].rsite = 'TH'
        self.Th.bonded_nblist = self.Th.guess_nblist(0,radii_buffer=1.25)
        self.Th.find_rsites()
    
        self.Th.write_xyz()
        self.calc = calculation.Calculation("calc001")
        

    def test_Thdist(self):

        self.calc_i = resource.CalculationRes("ThDist")
        file_type = 'input'
        file_key = 'xyz'
        file_name = "thiophene.xyz"
        self.calc_i.add_file(file_type,file_key,file_name)
        self.res_i = resource.Resource('res1',res_type = 'local')
        self.res_i.make_dir()
        self.struc_i = Buildingblock()
        self.calc_i.set_resource(self.res_i)
        self.calc_i.make_dir()
        
        file_dir = os.path.dirname(__file__)
        file_type = 'input'
        file_key = 'xyz'
        file_name = "thiophene.xyz"
        from_dirkey = 'home'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        self.calc_i.compress_dirfiles('input','scratch')
        
        os.chdir(self.calc_i.dir['scratch'])
        input_file = self.calc_i.files['input']['xyz']
        self.struc_i.read_xyz(input_file)
        self.struc_i.lat_cubic(10.0)
        self.struc_i.bonded_nblist = self.struc_i.guess_nblist(0,radii_buffer=1.25)
        
        # Set up log file 
        log_file = "%s.log"%self.calc_i.tag
        self.calc_i.add_file('output','log',log_file)
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        hdlr = logging.FileHandler(log_file,mode='w')
        hdlr.setLevel(logging.DEBUG)
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        prop = 'cosine'
        units = None
        self.calc_i.units[prop] = units
        logger.info(" Setting prop %s to units %s "%(prop,str(units)))
        dihcalc_tag = 'all_cos'
        self.struc_i.bonded_dih()
        self.struc_i.calc_dihedrals()
        outf = '%s.csv'%(dihcalc_tag)
        self.calc_i.add_file('data',dihcalc_tag,outf)
        self.struc_i.write_dihedrals(outf)
        # Get output
        file_key = self.calc_i.properties['comp_key']
        from_dirkey = 'scratch'
        to_dirkey = 'storage'
        file_type = 'output'
        self.calc_i.compress_dirfiles(file_type,from_dirkey)
        file_name = self.calc_i.files[file_type][self.calc_i.properties['comp_key']]
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        file_type = 'data'
        self.calc_i.compress_dirfiles(file_type,from_dirkey)
        file_name = self.calc_i.files[file_type][self.calc_i.properties['comp_key']]
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
        os.chdir(self.calc_i.dir['home'])

        

    @tearDown_streamm 
    def tearDown(self):
        del self.calc
        del self.Th


if __name__ == '__main__':

    unittest.main()    

        
