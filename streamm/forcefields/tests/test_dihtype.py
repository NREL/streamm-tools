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

import streamm.forcefields.dihtype as dihtype
            
from streamm_testutil import *




class Testdihtypeharmonic(unittest.TestCase):
    
    @setUp_streamm 
    def setUp(self):
        self.dihtype_i = dihtype.Dihtype("HC","CH","CH","HC",type="harmonic")
        self.dihtype_i.d = 4.0
        self.dihtype_i.mult = 3.0
        self.dihtype_i.theat_s = 45.0
        self.dihtype_i.kb = 80.6
        

    def test_dihstr(self):
        dih_str = ' dihedral  HC - CH - CH - HC type harmonic \n  harmonic d = 4.000000 mult = 3.000000 K = 80.600000 theat_s = 45.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    @tearDown_streamm 
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None


class Testdihtypemultiharmonic(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.dihtype_i = dihtype.Dihtype("HC","CH","CH","HC",type="multiharmonic")
        self.dihtype_i.d = 4.0
        self.dihtype_i.mult = 3.0
        self.dihtype_i.theat_s = 45.0
        self.dihtype_i.kb = 80.6
        

    def test_dihstr(self):
        dih_str = ' dihedral  HC - CH - CH - HC type multiharmonic \n  harmonic d = 4.000000 mult = 3.000000 K = 80.600000 theat_s = 45.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    @tearDown_streamm 
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None


class Testdihtypeopls(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.dihtype_i = dihtype.Dihtype("HC","CH","CH","HC",type="opls")
        self.dihtype_i.setopls(14.0,1.0,45.0,100.0)

    def test_dihstropls(self):
        dih_str = ' dihedral  HC - CH - CH - HC type opls \n  k1 = 14.000000 k2 = 1.000000 k3 = 45.000000 k4 = 100.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)

    def test_dihstrrb(self):
        self.dihtype_i.type = "rb"
        dih_str = ' dihedral  HC - CH - CH - HC type rb \n  C0 = 30.500000  C1 = 60.500000 C2 = 179.000000 C3 = -90.000000 C4 = -400.000000  C5 = 0.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
                
    @tearDown_streamm 
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None

class Testdihtyperb(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.dihtype_i = dihtype.Dihtype("HC","CH","CH","HC",type="rb")        
        self.dihtype_i.setrb(0.1,23.4,73.1,32.5,66.7,55.0)        

        
    def test_dihstrrb(self):
        dih_str = ' dihedral  HC - CH - CH - HC type rb \n  C0 = 0.100000  C1 = 23.400000 C2 = 73.100000 C3 = 32.500000 C4 = 66.700000  C5 = 55.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    def test_dihstropls(self):
        self.dihtype_i.type = "opls"
        dih_str = ' dihedral  HC - CH - CH - HC type opls \n  k1 = -95.550000 k2 = -139.800000 k3 = -16.250000 k4 = -16.675000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    @tearDown_streamm 
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None


if __name__ == '__main__':

    unittest.main()    

        
                

