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

import streamm.forcefields.parameters as parameters 
import streamm.forcefields.particletype as particletype
import streamm.forcefields.bondtype as bondtype
import streamm.forcefields.angletype as angletype
import streamm.forcefields.dihtype as dihtype
import streamm.forcefields.imptype as imptype
                        
                        
from streamm_testutil import *



class TestParameter(unittest.TestCase):
    @setUp_streamm 
    def setUp(self):
        self.paramC = parameters.Parameters()

    def test_str(self):
        empty_paramC_str = '\n    Parameters \n      LJ parameters 0 \n      Bond parameters 0 \n      Angle parameters 0 \n      Dihedral parameters 0 \n      Imporper Dihedral parameters 0 \n'
        self.assertEqual(str(self.paramC) ,empty_paramC_str)


    def test_LJ(self):
        particle_str = []
        self.particletype_i = particletype.Particletype("Ir")
        self.particletype_i.epsilon = 2.35
        self.particletype_i.sigma = 4.15
        self.paramC.add_particletype(self.particletype_i)
        particle_str.append(' Ir epsilon:2.35 sigma:4.15')
        self.particletype_i = particletype.Particletype('C*')
        self.particletype_i.epsilon = 1.05
        self.particletype_i.sigma = 3.25
        self.paramC.add_particletype(self.particletype_i)
        particle_str.append(' C* epsilon:1.05 sigma:3.25')
    	self.particletype_i = particletype.Particletype("HH")
        self.particletype_i.epsilon = 0.75
        self.particletype_i.sigma = 3.15
        self.paramC.add_particletype(self.particletype_i)
        particle_str.append(' HH epsilon:0.75 sigma:3.15')
        for particletkey_i, particletype_i  in self.paramC.particletypes.iteritems():
            self.assertEqual(str(particletype_i),particle_str[particletkey_i])
        


    def test_bond(self):
        bond_str = []
        self.bondtype_i = bondtype.Bondtype("Ir","C")
        self.bondtype_i.r0 = 1.02
        self.bondtype_i.kb = 13.563
        self.paramC.add_bondtype(self.bondtype_i)
        bond_str.append(' bond  Ir - C type harmonic \n  harmonic r_0 = 1.020000 K = 13.563000 lammps index 0  gromacs index 0  ')

        self.bondtype_i = bondtype.Bondtype("C","C")
        self.bondtype_i.r0 = 0.56
        self.bondtype_i.kb = 24.023
        self.paramC.add_bondtype(self.bondtype_i)
        bond_str.append(' bond  C - C type harmonic \n  harmonic r_0 = 0.560000 K = 24.023000 lammps index 0  gromacs index 0  ')

        self.bondtype_i = bondtype.Bondtype("C","H")
        self.bondtype_i.r0 = 0.43
        self.bondtype_i.kb = 65.123
        self.paramC.add_bondtype(self.bondtype_i)
        bond_str.append(' bond  C - H type harmonic \n  harmonic r_0 = 0.430000 K = 65.123000 lammps index 0  gromacs index 0  ')
        for btkey_i,bondtype_i  in self.paramC.bondtypes.iteritems():
            self.assertEqual(str(bondtype_i),bond_str[btkey_i])


    def test_angle(self):
        angle_str = []
        self.angletype_i = angletype.Angletype("H","C","H")
        self.angletype_i.theta0 = 120.0
        self.angletype_i.kb = 4.56
        angle_str.append(' angle  H - C - H type harmonic \n  harmonic theta_0 = 120.000000 K = 4.560000 lammps index 0  gromcas index 0  ')
        
        self.angletype_i = angletype.Angletype("C","Ir","C")
        self.angletype_i.theta0 = 90.0
        self.angletype_i.kb = 2.86
        angle_str.append(' angle  C - Ir - C type harmonic \n  harmonic theta_0 = 90.000000 K = 2.860000 lammps index 0  gromcas index 0  ')

        self.angletype_i = angletype.Angletype("Ir","C","H")
        self.angletype_i.theta0 = 120.0
        self.angletype_i.kb = 1.73
        angle_str.append(' angle  Ir - C - H type harmonic \n  harmonic theta_0 = 120.000000 K = 1.730000 lammps index 0  gromcas index 0  ')

        for atkey_i,angletype_i  in self.paramC.angletypes.iteritems():
            self.assertEqual(str(angletype_i),angle_str[atkey_i])


    def test_dih(self):
        dih_str = []
        self.dihtype_i = dihtype.Dihtype("H","C","Ir","H",type="harmonic")
        self.dihtype_i.d = 4.0
        self.dihtype_i.mult = 3.0
        self.dihtype_i.theat_s = 45.0
        self.dihtype_i.kb = 80.6
        dih_str.append(' dihedral  HC - CH - CH - HC type harmonic \n  harmonic d = 4.000000 mult = 3.000000 K = 80.600000 theat_s = 45.000000 lammps index 0  gromcas index 0 ')
        for dtkey_i, dihtype_i  in self.paramC.dihtypes.iteritems():
            self.assertEqual(str(dihtype_i),dih_str[dtkey_i])
            
        
    def test_imp(self):
        imp_str = []
        self.imptype_i = imptype.Imptype("Ir","C","C","C",type="harmonic")
        self.imptype_i.e0 = 180.0
        self.imptype_i.ke = 67.3
        self.imptype_i.pn = 4.0
        imp_str.append(' improper  Ir - C - C - C type harmonic \n  imp e0 = 180.000000 ke = 67.300000 lammps index 0  gromcas index 0 ')
        for itkey_i, imptype_i  in self.paramC.imptypes.iteritems():    
            self.assertEqual(str(imptype_i),imp_str[itkey_i])
        
    @tearDown_streamm 
    def tearDown(self):
        del self.paramC 
        self.paramC = None


if __name__ == '__main__':

    unittest.main()    

        
                
