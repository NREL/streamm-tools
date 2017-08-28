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

try:
    import streamm.structure.containers as containers
    import streamm.structure.groups as groups
            
    import streamm.structure.lattices as lattices
    import streamm.structure.nblists as nblists
    
    import streamm.structure.atoms as atoms
    import streamm.structure.bonds as bonds 
    import streamm.structure.angles as angles
    import streamm.structure.dihedrals as dihedrals
    import streamm.structure.impropers as impropers 
    
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    import containers
    import groups
    import lattices
    import nblists
    import atoms
    import bonds 
    import angles 
    import dihedrals 
    import impropers     
    

class TestGroupsProps(unittest.TestCase):
    # 
    def setUp(self):
        
        self.th = containers.Container("thiophene")
            
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
            pt_i = atoms.Atom(symbols[i])
            if( symbols[i] == 'C' ):
                pt_i.charge = -0.75
            elif(  symbols[i] == 'H' ):
                pt_i.charge = 0.80
            
            pt_i.mol = 0
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.th.add_partpos(pt_i,pos_i)
        
        
        self.th.lat_cubic(100.0)
        #
        for pkey_i, particle_i  in self.th.particles.iteritems():
            if( particle_i.tag == 'C' ):
                particle_i.resname = "SCP2"
                particle_i.residue = 1
            if( particle_i.tag == 'S' ):
                particle_i.resname = "ThS"
                particle_i.residue = 2
            if( particle_i.tag == 'H' ):
                particle_i.resname = "HA"
                particle_i.residue = 3
        self.strucC = containers.Container()
        self.strucC.lat_cubic(100.0)
        seed = 82343
        self.strucC = self.strucC.add_struc(self.th,10,seed)

    def test_molnumbers(self):
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            mol_i = int(pkey_i/self.th.n_particles) 
            self.assertEqual(str(particle_i.mol),str(mol_i))

    def test_groupmol(self):
        
        group_tag = 'mol'
        groupContainer_i = groups.Container(group_tag,self.strucC )
        groupContainer_i.group_prop('mol',group_tag)
        
        self.assertEqual(str(len(groupContainer_i.groups)),str(10))
        groupContainer_i.calc_cent_mass()
        groupContainer_i.calc_cent_mass()
        
        cm = []
        cm.append('[ 61.022463  12.212374  55.579404]')
        cm.append('[ 94.589545   0.548985  40.058567]')
        cm.append('[ 13.025619  22.458819  96.090279]')
        cm.append('[ 45.93974   30.752004  73.031331]')
        cm.append('[ 28.945124  70.792119  10.476723]')
        cm.append('[ 26.732501  56.981684  23.793239]')
        cm.append('[ 48.205917  63.191955  94.038944]')
        cm.append('[ 28.343741  95.032088  28.668735]')
        cm.append('[ 83.906182   8.100332  26.885987]')
        cm.append('[ 97.987557  38.078986  85.843074]')

        groupContainer_i.calc_radius()
        groupContainer_i.calc_radius()
        groupContainer_i.calc_radius()
        for gkey,group_i in groupContainer_i.groups.iteritems():
            self.assertEqual(str(group_i.cent_mass),str(cm[gkey]))
            self.assertEqual(str(group_i.radius),'2.57775210944')
            self.assertEqual(str(group_i.r_gy_sq),'3.6041389371')
            # print "r_gy.append(\'%s\')"%str(group_i.properties)
        for gkey in groupContainer_i.groups.keys:
            self.assertEqual(str(groupContainer_i.cent_mass[gkey]),str(cm[gkey]))
            self.assertEqual(str(groupContainer_i.radius[gkey]),'2.57775210944')
            self.assertEqual(str(groupContainer_i.r_gy_sq[gkey]),'3.6041389371')
            
        groupContainer_i.group_pbcs()

        os.chdir(os.path.dirname(__file__))
        groupContainer_i.write_cm_xyz()
        groupContainer_i.write_xyzs()
        groupContainer_i.dump_json()

        
    def test_groupres(self):
        group_tag = 'residue'
        groupContainer_i = groups.Container(group_tag,self.strucC )
        groupContainer_i.group_prop('residue',group_tag)
        
        self.assertEqual(str(len(groupContainer_i.groups)),str(30))
        
        groupContainer_i.calc_cent_mass()
        groupContainer_i.calc_radius()

        #for gkey,group_i in groupContainer_i.groups.iteritems():
        self.assertEqual(round(groupContainer_i.radius[2],6),2.587885)
        self.assertEqual(round(groupContainer_i.r_gy_sq[2],6),4.967159)
        self.assertEqual(round(groupContainer_i.Q_mn[2][0][0],6),0.005067)
        self.assertEqual(round(groupContainer_i.Rgy_eignval[0][0],6),1.002185)
        self.assertEqual(round(groupContainer_i.Rgy_eignval[0][1],6),0.410354)
        self.assertEqual(round(groupContainer_i.Rgy_eignval[0][2],6),0.0)
        self.assertEqual(round(groupContainer_i.A_sphere[0],6),0.381661)
        self.assertEqual(round(groupContainer_i.A_sphere_num[0],6),1.52303)
        self.assertEqual(round(groupContainer_i.A_sphere_dem[0],6),1.995267)
        
        groupContainer_i.calc_dl()

        self.assertEqual(round(groupContainer_i.dl_sq[0],6),5.966521)
        self.assertEqual(round(groupContainer_i.dl_sq[2],6),20.729214)
                
        os.chdir(os.path.dirname(__file__))
        groupContainer_i.write_cm_xyz()
        groupContainer_i.write_xyzs()
        groupContainer_i.dump_json()
        
    def tearDown(self):
        del self.th         
        del self.strucC
        
class TestGroupsHtermSp2(unittest.TestCase):
    # 
    def setUp(self):
        
        
        self.th = containers.Container("thiophene")
        

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
            pt_i = atoms.Atom(symbols[i])
            pos_i = positions[i]
            self.th.add_partpos(pt_i,pos_i)
        
        
        self.th.lat_cubic(100.0)
        # If no bonds guess based on radii 
        self.strucC.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
        # Build bonds from nblist for reference 
        self.th.bonded_bonds()
        #
        for pkey_i, particle_i  in self.th.particles.iteritems():
            if( particle_i.tag == 'C' ):
                particle_i.resname = "ThSC"
                particle_i.residue = 1
            if( particle_i.tag == 'S' ):
                particle_i.resname = "ThSC"
                particle_i.residue = 1
            if( particle_i.tag == 'H' ):
                particle_i.resname = "HA"
                particle_i.residue = 3
        group_tag = 'residue'
        groupContainer_i = groups.Container(group_tag,self.th )
        groupContainer_i.group_prop('residue',group_tag)

        self.assertEqual(len(self.groupContainer_i.groups),2)  
        
    def test_hterm(self):

        os.chdir(os.path.dirname(__file__))
        group_i =  self.groupContainer_i.groups[0]
        group_i.write_xyz('Th_SC.xyz')
        hterm_i = group_i.hterm_group()
        hterm_i.write_xyz('Th_SC_hterm.xyz')
        
    def tearDown(self):
        del self.th         

class TestGroupsHtermSp3(unittest.TestCase):
    # 
    def setUp(self):
        self.struc_i = containers.Container('ethane')
        symbols = ['C','C','H','H','H','H','H','H']
        positions = [ ]
        positions.append([-3.29091,-1.65766,-0.00000])
        positions.append([-2.35783,-0.47894,-0.00000])
        positions.append([-4.16830,-1.39014,0.58763])
        positions.append([-2.76492,-2.50106,0.44575])
        positions.append([-3.56295,-1.86934,-1.03338])
        positions.append([-2.08579,-0.26727,1.03338])
        positions.append([-2.88382,0.36446,-0.44575])
        positions.append([-1.48044,-0.74646,-0.58763])
        for i in range(len(symbols)):
            pt_i = atoms.Atom(symbols[i])
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.struc_i.add_partpos(pt_i,pos_i)
        self.struc_i.lat_cubic(100.0)
        # If no bonds guess based on radii 
        self.struc_i.bonded_nblist = self.struc_i.guess_nblist(0,radii_buffer=1.25)
        # Build bonds from nblist for reference 
        self.struc_i.bonded_bonds()
        #
    def test_hterm1(self):

        rmHcnt = 0 
        for pkey_i, particle_i  in self.struc_i.particles.iteritems():
            particle_i.resname = 'CRES'
            particle_i.residue = 1
            if( particle_i.tag == 'H' and rmHcnt < 1 ):
                particle_i.resname = 'HRES'
                particle_i.residue = 3
                rmHcnt += 1
                
        group_tag = 'residue'
        groupContainer_i = groups.Container(group_tag,self.struc_i )
        groupContainer_i.group_prop('residue',group_tag)
        self.assertEqual(len(self.groupContainer_i.groups),2)  
        
        os.chdir(os.path.dirname(__file__))
        group_i =  self.groupContainer_i.groups[0]
        group_i.write_xyz('Eth_C.xyz')
        hterm_i = group_i.hterm_group()

        for pkey_i, particle_i  in hterm_i.particles.iteritems():
            if( particle_i.tag == 'C' ):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),4)  
            if( particle_i.tag == 'H'):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),1)  
                        
        hterm_i.write_xyz('Eth_C_hterm1.xyz')
        
    def test_hterm2(self):

        rmHcnt = 0 
        for pkey_i, particle_i  in self.struc_i.particles.iteritems():
            particle_i.resname = 'CRES'
            particle_i.residue = 1
            if( particle_i.tag == 'H' and rmHcnt < 2 ):
                particle_i.resname = 'HRES'
                particle_i.residue = 3
                rmHcnt += 1
            
        group_tag = 'residue'
        groupContainer_i = groups.Container(group_tag,self.struc_i )
        groupContainer_i.group_prop('residue',group_tag)
        self.assertEqual(len(self.groupContainer_i.groups),2)  
        
        os.chdir(os.path.dirname(__file__))
        group_i =  self.groupContainer_i.groups[0]
        group_i.write_xyz('Eth_C.xyz')
        hterm_i = group_i.hterm_group()

        for pkey_i, particle_i  in hterm_i.particles.iteritems():
            if( particle_i.tag == 'C' ):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),4)  
            if( particle_i.tag == 'H'):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),1)  
                        
        hterm_i.write_xyz('Eth_C_hterm2.xyz')
        
    def tearDown(self):
        del self.struc_i         


class TestGroup_dr(unittest.TestCase):
    # 
    def setUp(self):
        self.th = containers.Container("thiophene")

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
            pt_i = atoms.Atom(symbols[i])
            pos_i = positions[i]
            self.th.add_partpos(pt_i,pos_i)
        

        self.th.lat_cubic(100.0)
        #
        for pkey_i, particle_i  in self.th.particles.iteritems():
            if( particle_i.tag == 'C' ):
                particle_i.resname = "SCP2"
                particle_i.residue = 1
            if( particle_i.tag == 'S' ):
                particle_i.resname = "ThS"
                particle_i.residue = 1
            if( particle_i.tag == 'H' ):
                particle_i.resname = "HA"
                particle_i.residue = 1
                
        self.strucC = containers.Container('th_x2')
        self.strucC.lat_cubic(30.0)
        seed = 82343
        self.strucC = self.strucC.add_struc(self.th,3,seed)
        self.strucC.tag = 'th_x3'
        #
    def test_finddr(self):
        self.strucC.lat_cubic(300.0)
        self.strucC.write_xyz()
        #self.strucC.write_cply()
        self.list_i = []
        for pkey,par_i in self.strucC.particles.iteritems():
            # print  pkey,par_i.mol,par_i.symbol 
            if( par_i.symbol == 'C' or par_i.symbol == 'S' ):
                self.list_i.append(pkey)
                print pkey ,par_i.mol , par_i.symbol 
                
        # self.strucC.bonded_nblist.build_nblist(self.strucC.particles,self.strucC.bonds)
        self.strucC.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)

        group_tag = 'mol'
        groupContainer_i = groups.Container(group_tag,self.struc_i )
        groupContainer_i.group_prop('mol',group_tag)
        
        groupContainer_i.calc_cent_mass()
        groupContainer_i.write_cm_xyz()
        groupContainer_i.calc_radius()     

        list_i = groupContainer_i.groups.keys()
        list_j = groupContainer_i.groups.keys()        
        
        group_i = groupContainer_i.groups[0]
        group_i.write_xyz()
        
        #pairbuffer = 2.5        
        pairs_ij = groupContainer_i.find_pairs(list_i,list_j,mol_inter=True,mol_intra=False)
        
        
        r_cut = 25.0
        bin_size = 0.10            
        close_contacts = True
        
        n_bins = int(r_cut/bin_size) + 1 
        bin_r = np.zeros(n_bins)    
        bin_r_nn = np.zeros(n_bins)    # Nearest neighbor count 
        bin_r_pp = np.zeros(n_bins)    
        probabilityperpair = 1
        volumes = []
        #
        N_i = len(list_i)
        N_j = len(list_j)
        
        
        self.strucC.calc_volume()
        volumes.append(self.strucC.volume)
        
        npos_i = groupContainer_i.cent_mass
        npos_j = groupContainer_i.cent_mass        
        npos_ij,nd_ij = self.strucC.lat.delta_npos(npos_i,npos_j)
        
        
        for ref_i in range(N_i):
            a_i_hasnieghbor = False
            r_ij_nn = r_cut   # Nearest Neighbor distance  
            g_i = list_i[ref_i]
            for ref_j in range(N_j):
                if(  pairs_ij[ref_i][ref_j] > 0.0 ):
                    dr_ij =  nd_ij[ref_i,ref_j]
                    if(  dr_ij <= r_cut ):
                            # bin distance =
                            bin_index = int( round( dr_ij / bin_size) )
                            #
                            # print " dist / bin / bin_sit", dist[ref_i,ref_j],bin_index,bin_size*float(bin_index)
                            #
                            bin_r[bin_index] += probabilityperpair
                            # Find nearest neighbor distance 
                            a_i_hasnieghbor = True
                            if( dr_ij < r_ij_nn ):
                                r_ij_nn = dr_ij
                                p_ij_nn = pairs_ij[ref_i][ref_j]
                            # 
                            if( close_contacts ):
                                g_j = list_i[ref_j]
                                dr_pi_pj = groupContainer_i.dr_particles(g_i,g_j,r_cut)
                                bin_pp_index = int( round( dr_pi_pj / bin_size) )
                                bin_r_pp[bin_pp_index] += probabilityperpair

            # Record nearest neighbor distance 
            if( a_i_hasnieghbor ):
                bin_nn_index = int( round( r_ij_nn /bin_size) )
                bin_r_nn[bin_nn_index] += p_ij_nn  
     
    def tearDown(self):
        del self.th         
        del self.strucC         


if __name__ == '__main__':
    unittest.main()
        
                