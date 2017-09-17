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

import streamm.structures.particle as particle
import streamm.structures.dihedral as dihedral

import streamm.structures.buildingblock as buildingblock

from streamm_testutil import *

class Test_attach(unittest.TestCase):

    @setUp_streamm 
    def setUp(self):
        
                                            
        self.Th = buildingblock.Buildingblock('thiophene')
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
            

        self.Hx = buildingblock.Buildingblock('hexane')
        
        symbols = []
        symbols.append('C')
        symbols.append('H')
        symbols.append('H')
        symbols.append('H')
        symbols.append('C')
        symbols.append('H')
        symbols.append('H')
        symbols.append('C')
        symbols.append('H')
        symbols.append('H')
        symbols.append('C')
        symbols.append('H')
        symbols.append('H')
        symbols.append('C')
        symbols.append('H')
        symbols.append('H')
        symbols.append('C')
        symbols.append('H')
        symbols.append('H')
        symbols.append('H')
        positions = []
        positions.append([-6.410969,-0.381641,-0.000031])
        positions.append([-7.310084,0.245311,-0.000038])
        positions.append([-6.456117,-1.028799,0.884636])
        positions.append([-6.456111,-1.028812,-0.884689])
        positions.append([-5.135268,0.467175,-0.000033])
        positions.append([-5.135484,1.128782,0.877977])
        positions.append([-5.135479,1.128771,-0.87805])
        positions.append([-3.850566,-0.371258,-0.000024])
        positions.append([-3.85112,-1.033978,0.87841])
        positions.append([-3.851114,-1.033987,-0.878451])
        positions.append([-2.567451,0.469603,-0.000024])
        positions.append([-2.567784,1.132155,0.8784])
        positions.append([-2.567776,1.132146,-0.878455])
        positions.append([-1.283527,-0.370234,-0.000013])
        positions.append([-1.28337,-1.032804,0.87836])
        positions.append([-1.28336,-1.032812,-0.87838])
        positions.append([0.00482234,0.47342231,-0.00000898])
        positions.append([0.02595107,1.09220686,0.87266464])
        positions.append([0.85585781,-0.17514133,0.00194589])
        positions.append([0.02780957,1.08937798,-0.87463473])

        for i in range(len(symbols)):
            pt_i = particle.Particle(symbol=symbols[i])
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.Hx.add_partpos(pt_i,pos_i)

        self.Hx.particles[1].rsite = 'R'
        self.Hx.bonded_nblist = self.Hx.guess_nblist(0,radii_buffer=1.25)
        self.Hx.find_rsites()
        
    def test_tag(self):
        self.assertEqual(self.Hx.tag,'hexane')
        self.assertEqual(self.Th.tag,'thiophene')

    def test_parsecplytag(self):
        self.assertEqual(self.Th.particles[0].rsite,'')
        self.assertEqual(self.Th.particles[1].rsite,'')
        self.assertEqual(self.Th.particles[2].rsite,'')
        self.assertEqual(self.Th.particles[3].rsite,'')
        self.assertEqual(self.Th.particles[4].rsite,'')
        self.assertEqual(self.Th.particles[5].rsite,'TH')
        self.assertEqual(self.Th.particles[6].rsite,'FH')
        self.assertEqual(self.Th.particles[7].rsite,'')
        self.assertEqual(self.Th.particles[8].rsite,'TH')
        

        self.assertEqual(self.Hx.particles[1].rsite,'R')

    def test_cnt(self):
        
        self.assertEqual(self.Hx.n_func,1)
        self.assertEqual(self.Th.n_func,3)
    
    def test_show_rsites(self):
        self.assertEqual(self.Hx.show_rsites(),'rsite:R[ paticle:atom[None] H (H) index:1 n_bonds:1] \n')
        self.assertEqual(self.Th.show_rsites(),'rsite:FH[ paticle:atom[None] H (H) index:6 n_bonds:1] \nrsite:TH[ paticle:atom[None] H (H) index:5 n_bonds:1] \nrsite:TH[ paticle:atom[None] H (H) index:8 n_bonds:1] \n')
        
    def test_get_rsite(self):
        Rkey_i,Xkey_i = self.Hx.get_rsite('R')
        self.assertEqual(Rkey_i,1)
        self.assertEqual(Xkey_i,0)
        
        Rkey_i,Xkey_i = self.Th.get_rsite('FH')
        self.assertEqual(Rkey_i,6)
        self.assertEqual(Xkey_i,1)
        
        Rkey_i,Xkey_i = self.Th.get_rsite('TH',n_i=0)
        self.assertEqual(Rkey_i,5)
        self.assertEqual(Xkey_i,0)
        
        Rkey_i,Xkey_i = self.Th.get_rsite('TH',n_i=1)
        self.assertEqual(Rkey_i,8)
        self.assertEqual(Xkey_i,3)
        
    def test_attachprpe(self):
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.Th.tag)
        self.Th.write_xyz(file_i)
        os.remove(file_i)
        
        self.Th2 = copy.deepcopy(self.Th)
        
        self.bb_i = self.Th2.prepattach("TH",0,0,dir=-1)
        self.bb_i.tag += "_"

        angle_rad = 90.0*math.pi/180.0
        self.bb_j = self.Th.prepattach("TH",1,0,dir=1,yangle=angle_rad)
        self.bb_j.tag += "_"
        

        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.bb_i.tag)
        self.bb_i.write_xyz(file_i)
        os.remove(file_i)
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.bb_j.tag)
        self.bb_j.write_xyz(file_i)
        os.remove(file_i)

        self.bbC_i,self.bbC_j =  buildingblock.shiftprep(self.bb_i,self.bb_j )
        self.overlap_found =  buildingblock.checkprep(self.bb_i,self.bb_j )

        self.assertTrue(self.overlap_found)


        self.Th2_Th =  buildingblock.attachprep(self.bbC_i,self.bbC_j )
        self.Th2_Th.tag = self.bb_i.tag + self.bb_j.tag + "v2"
        self.Th2_Th.lat_cubic(100.0)
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.Th2_Th.tag)
        self.Th2_Th.write_xyz(file_i)
        os.remove(file_i)
                
                

        n_bonds = self.Th2_Th.n_bonds
        bond_i = self.Th2_Th.bonds[n_bonds-1]
        self.key_i = bond_i.pkey1
        self.key_j = bond_i.pkey2

        self.assertEqual(self.key_i,0)
        self.assertEqual(self.key_j,11)


        self.dr_ij,self.mag_dr_ij = self.Th2_Th.lat.delta_pos(self.Th2_Th.positions[self.key_i],self.Th2_Th.positions[self.key_j])
        self.assertEqual(round(self.mag_dr_ij,2),1.34)
        
        nn_i =  self.Th2_Th.bonded_nblist.calc_nnab(self.key_i)
        nn_j =  self.Th2_Th.bonded_nblist.calc_nnab(self.key_j)

        self.assertEqual(nn_i,3)
        self.assertEqual(nn_j,3)
        
        # Find attached carbons 
        for self.key_k in self.Th2_Th.bonded_nblist.getnbs(self.key_i):
            if( self.key_k != self.key_j ):
                particle_k = self.Th2_Th.particles[self.key_k]
                if( particle_k.symbol == 'C' ): break
        for self.key_l in self.Th2_Th.bonded_nblist.getnbs(self.key_j):
            if( self.key_l != self.key_i ):
                particle_l = self.Th2_Th.particles[self.key_l]
                if( particle_l.symbol == 'C' ): break

            
        self.assertEqual(self.key_k,1)
        self.assertEqual(self.key_l,10)
        # Find angle
        dih_i = dihedral.Dihedral(self.key_k,self.key_i,self.key_j,self.key_l)
        
        self.cos_kijl = self.Th2_Th.calc_dihedral(dih_i)

        self.assertEqual(round(self.cos_kijl,6),0.0)
                
        

    def test_cat1(self):
        self.bblockC_p3htn1 = buildingblock.attach(self.Th,self.Hx,"FH",0,"R",0)
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.bblockC_p3htn1.tag)
        self.bblockC_p3htn1.write_xyz(file_i)
        os.remove(file_i)

    def test_cat2(self):

        self.bblockC_ptn2 = buildingblock.attach(self.Th,self.Th,"TH",0,"TH",1,tag = "P3HT_n1")
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.bblockC_ptn2.tag)
        self.bblockC_ptn2.write_xyz(file_i)
        os.remove(file_i)

        
    def test_cat3(self):
        self.bblockC_p3htn1 = buildingblock.attach(self.Th,self.Hx,"FH",0,"R",0,tag = "P3HT_nX")
        self.bblockC_p3htnX = buildingblock.Buildingblock()
        self.bblockC_p3htnX += self.bblockC_p3htn1 # Same as self.bblockC_p3htnX = copy.deepcopy( self.bblockC_p3htn1 )
        # for n in range(1,xN):
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"TH",0,"TH",1)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"TH",0,"TH",1)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"TH",1,"TH",0)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"TH",0,"TH",1)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"TH",0,"TH",1)

        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.Hx,"TH",0,"R",0)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.Hx,"TH",0,"R",0)
        
        self.bblockC_p3htnX.tag = "P3HT_nX"
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.bblockC_p3htnX.tag)
        self.bblockC_p3htnX.write_xyz(file_i)
        os.remove(file_i)
        
    def test_p3ht_nX(self):

        self.p3htn1 = buildingblock.attach(self.Th,self.Hx,"FH",0,"R",0)
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.p3htn1.tag)
        self.p3htn1.write_xyz(file_i)
        os.remove(file_i)
        
        # Create a xN membered chain
        xN = 5
        #self.bb_p3ht_nX = buildingblock.Buildingblock()
        #self.bb_p3ht_nX.tag = "p3ht_n%d"%(xN)
        self.bb_p3ht_nX = copy.deepcopy(self.p3htn1)

        for n in range(1,xN):
            # print " Adding %d "%(n)
            self.bb_p3ht_nX = buildingblock.attach(self.bb_p3ht_nX,self.p3htn1,"TH",0,"TH",1,tag="p3ht_n%d"%(n))

        # Print xyz file and cply file 
        file_i = os.path.join(TEST_DIR, "%s.xyz"%self.bb_p3ht_nX.tag)
        self.bb_p3ht_nX.write_xyz(file_i)
        os.remove(file_i)

    @tearDown_streamm 
    def tearDown(self):
        del self.Th
        del self.Hx



if __name__ == '__main__':
    unittest.main()    
        
