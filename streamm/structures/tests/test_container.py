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
    import streamm.structures.container as container
            
    import streamm.structures.lattice as lattice
    import streamm.structures.nblist as nblist
    
    import streamm.structures.particle as particle
    import streamm.structures.bond as bond
    import streamm.structures.angle as angle
    import streamm.structures.dihedral as dihedral
    import streamm.structures.improper as improper
    
except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    import container
    import lattice
    import nblist
    import particle
    import bond 
    import angle
    import dihedral 
    import improper
    

class TestContainer(unittest.TestCase):
    def setUp(self):
        self.strucC = container.Container("Test_structure1")
        
    def test_str(self):
        self.assertEqual(str(self.strucC)," Test_structure1")

    def test_Lattice(self):
        matrix = [ 132,0,0,0,127,0,0,0,150 ]
        self.strucC.lat.set_matrix(matrix)
        
        self.assertEqual(self.strucC.lat._lengths[0],132.0)
        self.assertEqual(self.strucC.lat._lengths[1],127.0)
        self.assertEqual(self.strucC.lat._lengths[2],150.0)

        self.assertEqual(self.strucC.lat._angles[0],90.0)
        self.assertEqual(self.strucC.lat._angles[1],90.0)
        self.assertEqual(self.strucC.lat._angles[2],90.0)

    def test_particles(self):
        self.part = particle.Particle(symbol="Ir")
        self.strucC.add_particle(self.part)
        pos_add = []
        pos_i = [ 0.0,0.0,0.0]
        pos_add.append(pos_i)
        self.strucC.add_position(pos_i)        
        self.assertEqual(self.strucC.n_particles,1)

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,0.0,1.5]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,2)

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,1.50,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,3)


        self.part = particle.Particle(symbol="C")
        pos_i = [ 1.50,0.0,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,0.0,-1.5]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,-1.50,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = particle.Particle(symbol="C")
        pos_i = [ -1.50,0.0,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))


        self.part = particle.Particle(symbol="H")
        pos_i = [ 0.0,0.750,1.75]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))
        
        self.part = particle.Particle(symbol="H")
        pos_i = [ 0.750,0.0,1.75]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = particle.Particle(symbol="H")
        pos_i = [- 0.53,-0.53,1.75]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        p_cnt = 0
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            pos_i = self.strucC.positions[pkey_i]
            
            self.assertEqual(pkey_i,p_cnt,"Particle keys are non sequential" )
            self.assertEqual(pos_i[0],pos_add[pkey_i][0])
            self.assertEqual(pos_i[1],pos_add[pkey_i][1])
            self.assertEqual(pos_i[2],pos_add[pkey_i][2])

            p_cnt +=1

    def test_bonds(self):
        bond_str = []        
        self.bond_i = bond.Bond(0,1)
        bond_str.append(' 0 - 1')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,2)
        bond_str.append(' 0 - 2')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,3)
        bond_str.append(' 0 - 3')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,4)
        bond_str.append(' 0 - 4')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,5)
        bond_str.append(' 0 - 5')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,6)
        bond_str.append(' 0 - 6')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(1,7)
        bond_str.append(' 1 - 7')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(1,8)
        bond_str.append(' 1 - 8')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(1,9)
        bond_str.append(' 1 - 9')
        self.strucC.add_bond(self.bond_i)
        for bkey_i, bond_i  in self.strucC.bonds.iteritems():
            self.assertEqual(str(bond_i),bond_str[bkey_i])

    def test_nblist(self):

        self.part = particle.Particle(symbol="Ir")
        pos_i = [ 0.0,0.0,0.0]
        self.strucC.add_partpos(self.part,pos_i)

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,0.0,1.5]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,1.50,0.0]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = particle.Particle(symbol="C")
        pos_i = [ 1.50,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,0.0,-1.5]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,-1.50,0.0]   
        self.strucC.add_partpos(self.part,pos_i)


        self.part = particle.Particle(symbol="C")
        pos_i = [ -1.50,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)


        self.part = particle.Particle(symbol="H")
        pos_i = [ 0.0,0.750,1.75]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.part = particle.Particle(symbol="H")
        pos_i = [ 0.750,0.0,1.75]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = particle.Particle(symbol="H")
        pos_i = [- 0.53,-0.53,1.75]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.bond_i = bond.Bond(0,1)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,2)
        self.strucC.add_bond(self.bond_i)
        # 
        self.bond_i = bond.Bond(0,3)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,4)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,5)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(0,6)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(1,7)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(1,8)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = bond.Bond(1,9)
        self.strucC.add_bond(self.bond_i)

        self.strucC.bonded_nblist = self.strucC.build_nblist()
        
        nb_str = []
        nb_str.append('0  -  1')
        nb_str.append('0  -  2')
        nb_str.append('0  -  3')
        nb_str.append('0  -  4')
        nb_str.append('0  -  5')
        nb_str.append('0  -  6')
        nb_str.append('1  -  0')
        nb_str.append('1  -  7')
        nb_str.append('1  -  8')
        nb_str.append('1  -  9')
        nb_str.append('2  -  0')
        nb_str.append('3  -  0')
        nb_str.append('4  -  0')
        nb_str.append('5  -  0')
        nb_str.append('6  -  0')
        nb_str.append('7  -  1')
        nb_str.append('8  -  1')
        nb_str.append('9  -  1')

        nb_cnt = 0 
        for pkey_i  in self.strucC.particles.keys():
            for pkey_j in self.strucC.bonded_nblist.getnbs(pkey_i):
                self.assertEqual(str(pkey_i)+"  -  "+str(pkey_j),nb_str[nb_cnt])
                nb_cnt += 1

    def test_angles(self):
        angle_str = []        
        self.angle_i = angle.Angle(2,0,1)
        angle_str.append(' 2 - 0 - 1')
        self.strucC.add_angle(self.angle_i)
        self.angle_i = angle.Angle(3,0,4)
        angle_str.append(' 3 - 0 - 4')
        self.strucC.add_angle(self.angle_i)
        self.angle_i = angle.Angle(5,0,6)
        angle_str.append(' 5 - 0 - 6')
        self.strucC.add_angle(self.angle_i)
        for akey_i, angle_i  in self.strucC.angles.iteritems():
            self.assertEqual(str(angle_i),angle_str[akey_i])


    def test_dihedrals(self):
        dihedral_str = []        
        self.dihedral_i = dihedral.Dihedral(2,0,1,7)
        dihedral_str.append(' 2 - 0 - 1 - 7')
        self.strucC.add_dihedral(self.dihedral_i)
        self.dihedral_i = dihedral.Dihedral(2,0,1,8)
        dihedral_str.append(' 2 - 0 - 1 - 8')
        self.strucC.add_dihedral(self.dihedral_i)
        self.dihedral_i = dihedral.Dihedral(2,0,1,9)
        dihedral_str.append(' 2 - 0 - 1 - 9')
        self.strucC.add_dihedral(self.dihedral_i)
        for dkey_i, dihedral_i  in self.strucC.dihedrals.iteritems():
            self.assertEqual(str(dihedral_i),dihedral_str[dkey_i])


    def test_impropers(self):
        improper_str = []        
        self.improper_i = improper.Improper(0,1,2,3)
        improper_str.append(' 0 - 1 - 2 - 3')
        self.strucC.add_improper(self.improper_i)
        self.improper_i = improper.Improper(0,4,5,6)
        improper_str.append(' 0 - 4 - 5 - 6')
        self.strucC.add_improper(self.improper_i)
        for ikey_i, improper_i  in self.strucC.impropers.iteritems():
            self.assertEqual(str(improper_i),improper_str[ikey_i])
        
        

    def test_shift(self):
        self.part = particle.Particle(symbol="C")
        pos_i = [ 1.0,2.0,1.5]   
        self.strucC.add_partpos(self.part,pos_i)
        self.strucC.shift( 0, np.array([-0.5,-1.5,-1.0]) )

        pos_j = self.strucC.positions[0]
        self.assertEqual(pos_j[0],0.50)
        self.assertEqual(pos_j[1],0.50)
        self.assertEqual(pos_j[2],0.50)


    def test_shift_pos(self):
        
        self.part = particle.Particle(symbol="C")
        pos_i = [ 1.0,2.0,1.5]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.part = particle.Particle(symbol="C")
        pos_i = [ 6.0,7.0,6.5]   
        self.strucC.add_partpos(self.part,pos_i)


        self.strucC.shift_pos( np.array([-0.5,-1.5,-1.0]) )

        pos_j = self.strucC.positions[0]
        self.assertEqual(pos_j[0],0.50)
        self.assertEqual(pos_j[1],0.50)
        self.assertEqual(pos_j[2],0.50)

        pos_j = self.strucC.positions[1]
        self.assertEqual(pos_j[0],5.50)
        self.assertEqual(pos_j[1],5.50)
        self.assertEqual(pos_j[2],5.50)
                
    def tearDown(self):
        del self.strucC 
        self.strucC = None



class Test_guessnbs(unittest.TestCase):
    
    def setUp(self):
        
        self.strucC = container.Container()

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)        
        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.750,0.750,0.750]   
        self.strucC.add_partpos(self.part,pos_i)
        self.part = particle.Particle(symbol="C")
        pos_i = [ -0.750,-0.750,0.750]   
        self.strucC.add_partpos(self.part,pos_i)
        self.part = particle.Particle(symbol="C")
        pos_i = [ -1.250,-1.250,0.0]   
        self.strucC.add_partpos(self.part,pos_i)

    def test_guessnb(self):
        # self.strucC.bonded_nblist.build_nblist(self.strucC.particles,self.strucC.bonds )
        str_nbs_list = []

        str_nbs_list.append(' 0 - 1')
        str_nbs_list.append(' 0 - 2')
        str_nbs_list.append(' 1 - 0')
        str_nbs_list.append(' 2 - 0')
        str_nbs_list.append(' 2 - 3')
        str_nbs_list.append(' 3 - 2')
        cnt = 0
        self.strucC.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
        
        for pkey_i, particle_i in self.strucC.particles.iteritems():
             for pkey_j in   self.strucC.bonded_nblist.getnbs(pkey_i):
                 # print " str_nbs_list.append(\' %d - %d \')"%(pkey_i,pkey_j)
                 self.assertEqual(str_nbs_list[cnt],' %d - %d'%(pkey_i,pkey_j) )
                 cnt += 1
        
    def tearDown(self):
        del self.strucC 
        self.strucC = None

class TestBuildThiophene(unittest.TestCase):

    def setUp(self):
        
        self.strucC = container.Container('thiophene')
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
            if( symbols[i] == 'C' ):
                pt_i.charge = -0.75
            elif(  symbols[i] == 'H' ):
                pt_i.charge = 0.80
            
            pt_i.mol = 0
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.strucC.add_partpos(pt_i,pos_i)
        

        matrix_i = self.strucC.lat._matrix
        matrix_i[0][0] = 100.0 
        matrix_i[1][1] = 100.0 
        matrix_i[2][2] = 100.0 
        self.strucC.lat.set_matrix(matrix_i)
        self.strucC.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
        #
        #
    def test_random(self):
        seed = 440514
        random.seed(seed)
        self.strucC.lat.random_pos()
        
    def test_rotate(self):
        self.strucC.calc_mass()
        self.strucC.calc_center_mass()
        self.strucC.shift_pos(-1.0*self.strucC.center_mass)  # Place center of mass at origin
        file_i = os.path.join(os.path.dirname(__file__), "rot_0.xyz")
        self.strucC.write_xyz(file_i)
        
        rot_angle = 1.5707963267949  # Angle in radians
        self.strucC.rotate_xy(rot_angle)
        file_i = os.path.join(os.path.dirname(__file__), "rot_1.xyz")
        self.strucC.write_xyz(file_i)
        
        rot_angle = -1.5707963267949  # Angle in radians
        self.strucC.rotate_xz(rot_angle)
        file_i = os.path.join(os.path.dirname(__file__), "rot_2.xyz")
        self.strucC.write_xyz(file_i)
        
        rot_angle = 1.5707963267949  # Angle in radians
        self.strucC.rotate_yz(rot_angle)
        file_i = os.path.join(os.path.dirname(__file__), "rot_3.xyz")
        self.strucC.write_xyz(file_i)

    def test_calcel(self):

        n_el = self.strucC.calc_elcnt(0,self.strucC.bonded_nblist)
        self.assertEqual(n_el[1],1)
        self.assertEqual(n_el[6],1)
        self.assertEqual(n_el[16],1)
        n_el = self.strucC.calc_elcnt(1,self.strucC.bonded_nblist)
        self.assertEqual(n_el[1],1)
        self.assertEqual(n_el[6],2)
        n_el = self.strucC.calc_elcnt(4,self.strucC.bonded_nblist)
        self.assertEqual(n_el[6],2)
        

    def test_bonds(self):
        self.strucC.bonded_bonds()
        self.strucC.calc_bonds()
        self.strucC.write_bonds('bonds_all.csv')
        self.assertEqual(round(self.strucC.bonds[0].length,6),1.377226)
        self.assertEqual(round(self.strucC.bonds[1].length,6),1.672168)
        self.assertEqual(round(self.strucC.bonds[2].length,6),1.088203)
        self.assertEqual(round(self.strucC.bonds[3].length,6),1.432118)
        self.assertEqual(round(self.strucC.bonds[4].length,6),1.091483)
        self.assertEqual(round(self.strucC.bonds[5].length,6),1.377194)


    def test_listbondlengths(self):
        PRECISION =4

        self.strucC.bonded_bonds()
        # Get C-C bond lengths
        list_i = []
        list_j = []
        for k,p in self.strucC.particles.iteritems():
            if( p.symbol == 'C' ):
                list_i.append(k)
                list_j.append(k)
        CC_bondkeys = self.strucC.find_bonds(list_i,list_j)
        self.strucC.calc_bonds(CC_bondkeys)
        CC_bondlengths = self.strucC.get_list_bonds_lengths(CC_bondkeys)
               
        self.assertEqual(round(CC_bondlengths[0],PRECISION),1.3772)
        self.assertEqual(round(CC_bondlengths[1],PRECISION),1.4321)
        self.assertEqual(round(CC_bondlengths[2],PRECISION),1.3772)
        outf = 'CCbonds.csv'
        self.strucC.write_bonds(outf,CC_bondkeys)
        # Get C-S bond lengths
        list_i = []
        list_j = []
        #
        for k,p in self.strucC.particles.iteritems():
                if( p.symbol == 'C' ):
                    list_i.append(k)
                if( p.symbol == 'S' ):
                    list_j.append(k)
        # 
        CS_bondkeys = self.strucC.find_bonds(list_i,list_j)
        self.strucC.calc_bonds(CS_bondkeys)
        CS_bondlengths = self.strucC.get_list_bonds_lengths(CS_bondkeys)
        
        self.assertEqual(round(CS_bondlengths[0],PRECISION),1.6722)
        self.assertEqual(round(CS_bondlengths[0],PRECISION),1.6722)
        outf = 'CSbonds.csv'
        self.strucC.write_bonds(outf,CS_bondkeys)
        
    def test_angles(self):
        self.strucC.bonded_angles()
        self.strucC.calc_angles()
        self.strucC.write_angles('angles_all.csv')
        self.assertEqual(round(self.strucC.angles[0].cosine,6),0.367484)
        self.assertEqual(round(self.strucC.angles[12].cosine,6),0.066883)
        
    def test_listbondangles(self):
        PRECISION =4

        self.strucC.bonded_angles()
        # Get C-C-C angles
        list_k = []
        list_i = []
        list_j = []
        for k,p in self.strucC.particles.iteritems():
                if( p.symbol == 'C' ):
                        list_k.append(k)
                        list_i.append(k)
                        list_j.append(k)
        CCC_keys = self.strucC.find_angles(list_k,list_i,list_j)
        self.strucC.calc_angles(CCC_keys)
        CCC_cos = self.strucC.get_list_angles_cosine(CCC_keys) # structure.prop_list('cosine',CCC_keys,self.strucC.angles)
        self.assertEqual(round(CCC_cos[0],PRECISION),0.3669)
        self.assertEqual(round(CCC_cos[1],PRECISION),0.3669)
        outf = 'CCC_cos.csv'
        self.strucC.write_angles(outf,CCC_keys)
        # Get C-S bond lengths
        list_k = []
        list_i = []
        list_j = []
        for k,p in self.strucC.particles.iteritems():
                if( p.symbol == 'C' ):
                        list_k.append(k)
                        list_j.append(k)
                if( p.symbol == 'S' ):
                        list_i.append(k)
        CSC_keys = self.strucC.find_angles(list_k,list_i,list_j)
        self.strucC.calc_angles(CSC_keys)
        
        CSC_cos = self.strucC.get_list_angles_cosine(CSC_keys) # structure.prop_list('cosine',CSC_keys,self.strucC.angles)
        self.assertEqual(round(CSC_cos[0],PRECISION),0.0669)
        outf = 'CSC_cos.csv'
        self.strucC.write_angles(outf,CSC_keys)
        
        
    def test_dih(self):
        self.strucC.bonded_dih()
        self.strucC.calc_dihedrals()
        self.strucC.write_dihedrals('dih_all.csv')
        self.assertEqual(round(self.strucC.dihedrals[8].cosine,6),-1.0)
        
    def test_calcdihc(self):
        dih_i = dihedral.Dihedral(3,2,1,0)
        self.strucC.calc_dihedral(dih_i)
        self.assertEqual(dih_i.cosine,0.99999980013324519)

    def test_listdihangle(self):
        PRECISION =4

        self.strucC.bonded_dih()
                     
        # Select some particles 
        list_k = []
        list_i = []
        list_j = []
        list_l = []

        for k,p in self.strucC.particles.iteritems():
                if( p.symbol == 'C' ):
                        list_k.append(k)
                        list_i.append(k)
                        list_j.append(k)
                        list_l.append(k)

        CCCC_keys = self.strucC.find_dihedrals(list_k,list_i,list_j,list_l)
        self.strucC.calc_dihedrals(CCCC_keys)
        CCCC_cos = self.strucC.get_list_dihedral_cosine(CCCC_keys) # structure.prop_list('cosine',CCCC_keys,self.strucC.dihedrals)
        self.assertEqual(round(CCCC_cos[0],PRECISION),1.0)
        outf = 'CCCC_cos.csv'
        self.strucC.write_dihedrals(outf,CCCC_keys)
        

    def test_expand_matrix(self):
        exlat_frac = 0.050  # 5%
        self.strucC.lat.expand_matrix(exlat_frac)
        self.assertEqual(self.strucC.lat._lengths[0],105.0)
        self.assertEqual(self.strucC.lat._lengths[1],105.0)
        self.assertEqual(self.strucC.lat._lengths[2],105.0)
            
    def test_tag(self):
        self.assertEqual(self.strucC.tag,'thiophene')
        #el_cnt = calc_elcnt
        
    def test_n_molecules(self):
        self.assertEqual(self.strucC.n_molecules(),0)

    def test_mol_mult(self): 
        self.strucC.mol_mult()
        self.assertEqual(self.strucC.mol_multiplier,100.0)
        

    def test_maxtags(self): 
        self.strucC.maxtags()
        self.assertEqual(self.strucC.mol_max,0)

    def test_dump_pickle(self):
        self.strucC.dump_pickle()
                                      

    def test_write_xyz_str(self):
        xyz_str = self.strucC.write_xyz_str()
        self.assertEqual(xyz_str,' 9 \n thiophene \n     C      -1.55498576      -1.91131218      -0.00081000 \n     C      -0.17775976      -1.91131218      -0.00081000 \n     C       0.34761524      -0.57904218      -0.00081000 \n     C      -0.65884476       0.36101082       0.00000000 \n     S      -2.16948076      -0.35614618      -0.00000800 \n     H      -2.18966076      -2.79526518      -0.00132100 \n     H       0.45389024      -2.80145418      -0.00106400 \n     H       1.41682424      -0.35961818      -0.00138200 \n     H      -0.51943676       1.44024682       0.00064700 \n')
        

    def test_write_xyz_list(self):
        list_i = [2,3,4,8]
        self.strucC.write_xyz_list(list_i,xyz_file='test_write_xyz_list.xyz')
                

    def test_write_xyz(self):
        # NoteTK os.chdir(os.path.dirname(__file__))
        self.strucC.write_xyz()

    def test_read_write_xyz(self): 
        self.strucC.write_xyz()
        self.strucC.read_xyz()

    def test_shift_pos(self): 
        vec = [-3.478324,0.234798234,-234.234]
        positions_correct = []
        self.strucC.shift_pos(vec)
        
        nptu.assert_almost_equal(self.strucC.positions[0] ,[  -5.03330976  , -1.67651395 ,-234.23481   ])
        nptu.assert_almost_equal(self.strucC.positions[1] ,[  -3.65608376  , -1.67651395 , -234.23481   ])
        nptu.assert_almost_equal(self.strucC.positions[2] ,[  -3.13070876  , -0.34424395 ,-234.23481   ])
        nptu.assert_almost_equal(self.strucC.positions[3] ,[  -4.13716876   , 0.59580905 ,-234.234     ])
        nptu.assert_almost_equal(self.strucC.positions[4] ,[ -5.64780476e+00 , -1.21347946e-01 , -2.34234008e+02])
        nptu.assert_almost_equal(self.strucC.positions[5] ,[  -5.66798476   ,-2.56046695 ,-234.235321  ])

                   
    def test_shift(self):
        pkey = 2
        vec = [34.32,13.1234,-218.3]
        self.strucC.shift(pkey,vec)
        nptu.assert_almost_equal(self.strucC.positions[pkey] ,[  34.6676152,   12.5443578, -218.30081 ])
                                        

    def test_calc_properties(self): 
        self.strucC.calc_mass()
        self.assertEqual(self.strucC.mass,84.13956)
 
        self.strucC.calc_charge()
        self.assertEqual(self.strucC.charge,0.2)
        self.strucC.lat_cubic(12.0)
        self.strucC.calc_volume()
        self.assertEqual(self.strucC.volume,1728.0)
        
        self.strucC.calc_density()
        self.assertEqual(self.strucC.density,0.048691875)
        self.strucC.calc_center_mass()
        
        nptu.assert_almost_equal(self.strucC.center_mass,[ -1.12858935e+00  , -7.66617741e-01 , -3.87300502e-04])
        
    def test_calc_composition(self): 
        self.strucC.calc_composition()
        self.assertEqual(self.strucC.composition[1],4)
        self.assertEqual(self.strucC.composition[6],4)
        self.assertEqual(self.strucC.composition[16],1)
        
        
    def test_calc_formula(self): 
        self.strucC.calc_formula()
        self.assertEqual(self.strucC.chemicalformula,'C4H4S1')
        
    def test_sum_prop(self):
        pkey_i = 5
        pkeyp_j = 0
        self.strucC.sum_prop(pkey_i,pkeyp_j)
        self.assertEqual(self.strucC.particles[pkey_i].charge,0.05)
                     

    def test_change_mass(self): 
        self.strucC.change_mass('C',13.02)
        self.assertEqual(self.strucC.particles[0].mass,13.02)
        self.assertEqual(self.strucC.particles[1].mass,13.02)
        self.assertEqual(self.strucC.particles[2].mass,13.02)
        self.assertEqual(self.strucC.particles[3].mass,13.02)

        
    def test_get_pos(self): 
        list_i = [0,2]
        npos_i_correct = np.array([[ -1.5549858e+00,  -1.9113122e+00,  -8.1000000e-04], 
       [  3.4761524e-01,  -5.7904218e-01,  -8.1000000e-04]])
        npos_i = self.strucC.get_pos(list_i)
        nptu.assert_almost_equal(npos_i,npos_i_correct)
        
    
    def test_writestruc(self):
        # print "self.strucC",self.strucC.n_particles
        part_xyz = []
        # 
        part_xyz.append('  C -1.554986 -1.911312 -0.000810  ')
        part_xyz.append('  C -0.177760 -1.911312 -0.000810  ')
        part_xyz.append('  C 0.347615 -0.579042 -0.000810  ')
        part_xyz.append('  C -0.658845 0.361011 0.000000  ')
        part_xyz.append('  S -2.169481 -0.356146 -0.000008  ')
        part_xyz.append('  H -2.189661 -2.795265 -0.001321  ')
        part_xyz.append('  H 0.453890 -2.801454 -0.001064  ')
        part_xyz.append('  H 1.416824 -0.359618 -0.001382  ')
        part_xyz.append('  H -0.519437 1.440247 0.000647  ')
        # 
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            pos_i = self.strucC.positions[pkey_i]
            particle_line = "  %s %f %f %f  "%(particle_i.symbol,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]))
            self.assertEqual(part_xyz[pkey_i],particle_line)
            # print " part_xyz.append(\'",particle_line,"\')"
        #                 
    def tearDown(self):
        del self.strucC 

class TestBuildEthane(unittest.TestCase):

    def setUp(self):
        
        self.Eth = container.Container('ethane')
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
            pt_i = particle.Particle(symbol=symbols[i])
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.Eth.add_partpos(pt_i,pos_i)

    def test_write_xyz(self):
        # NoteTK os.chdir(os.path.dirname(__file__))
        self.Eth.write_xyz()
                    
    def tearDown(self):
        del self.Eth 
            
class Testiadd(unittest.TestCase):
    
    def setUp(self):
        self.strucC = container.Container()
        self.strucC_j = container.Container()

        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.0,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.part = particle.Particle(symbol="C")
        pos_i = [ 0.750,0.750,0.750]   
        self.strucC.add_partpos(self.part,pos_i)
        # 
        self.strucC.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
        self.strucC.bonded_bonds()
        
        self.part = particle.Particle(symbol="C")
        pos_i = [ -0.750,-0.750,0.750]   
        self.strucC_j.add_partpos(self.part,pos_i)

        self.part = particle.Particle(symbol="C")
        pos_i = [ -1.250,-1.250,0.0]   
        self.strucC_j.add_partpos(self.part,pos_i)

        self.strucC_j.bonded_nblist = self.strucC_j.guess_nblist(0,radii_buffer=1.25)
        self.strucC_j.bonded_bonds()
        
        self.strucC += self.strucC_j

    def test_iaddx3(self):
        self.strucC += self.strucC_j
        self.strucC += self.strucC_j
        self.strucC += self.strucC_j
        self.strucC += self.strucC_j

        
    def test_join(self):

        pos_str_list = []
        pos_str_list.append(' 0 - C 0.000000 0.000000 0.000000 ')
        pos_str_list.append(' 1 - C 0.750000 0.750000 0.750000 ')
        pos_str_list.append(' 2 - C -0.750000 -0.750000 0.750000 ')
        pos_str_list.append(' 3 - C -1.250000 -1.250000 0.000000 ')
        
        str_nbs_list = []

        str_nbs_list.append(' 0 - 1 ')
        str_nbs_list.append(' 1 - 0 ')
        str_nbs_list.append(' 2 - 3 ')
        str_nbs_list.append(' 3 - 2 ')

        cnt = 0
        for pkey_i, particle_i in self.strucC.particles.iteritems():
            self.assertEqual(pos_str_list[pkey_i],' %d - %s %f %f %f '%(pkey_i,particle_i.symbol,self.strucC.positions[pkey_i][0],self.strucC.positions[pkey_i][1],self.strucC.positions[pkey_i][2]))
            for pkey_j in   self.strucC.bonded_nblist.getnbs(pkey_i):
                #print " str_nbs_list.append(\' %d - %d \')"%(pkey_i,pkey_j)
                self.assertEqual(str_nbs_list[cnt],' %d - %d '%(pkey_i,pkey_j) )
                cnt += 1

    def tearDown(self):
        del self.strucC 
        del self.strucC_j
        self.strucC = None
        self.strucC_j = None

class TestProximityCheck(unittest.TestCase):
    # 
    def setUp(self):
        self.strucC1 = container.Container()
        
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
            if( symbols[i] == 'C' ):
                pt_i.charge = -0.75
            elif(  symbols[i] == 'H' ):
                pt_i.charge = 0.80
            
            pt_i.mol = 0
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.strucC1.add_partpos(pt_i,pos_i)
        
        
        matrix_i = self.strucC1.lat._matrix
        matrix_i[0][0] = 1000.0 
        matrix_i[1][1] = 1000.0 
        matrix_i[2][2] = 1000.0 
        self.strucC1.lat.set_matrix(matrix_i)
        
        self.strucC2 = container.Container()

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
            if( symbols[i] == 'C' ):
                pt_i.charge = -0.75
            elif(  symbols[i] == 'H' ):
                pt_i.charge = 0.80
            
            pt_i.mol = 0
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.strucC2.add_partpos(pt_i,pos_i)
        
        matrix_i = self.strucC2.lat._matrix
        matrix_i[0][0] = 100.0 
        matrix_i[1][1] = 100.0 
        matrix_i[2][2] = 100.0 
        self.strucC2.lat.set_matrix(matrix_i)
        
        

                    
        #
    def test_poxcheck1(self):
        npos_i = self.strucC1.positions
        npos_j = self.strucC2.positions
        pos_cut = 2.0 # minimum distance between particles of added structure and current s
        poxflag =  self.strucC1.lat.proximitycheck(npos_i,npos_j,pos_cut)
        self.assertFalse(poxflag)
                                    
        #
    def test_poxcheck2(self):
        self.strucC1.shift_pos(np.array([25.0,0.0,0.0]))
        npos_i = self.strucC1.positions
        npos_j = self.strucC2.positions
        pos_cut = 2.0 # minimum distance between particles of added structure and current s
        poxflag =  self.strucC1.lat.proximitycheck(npos_i,npos_j,pos_cut)
        self.assertTrue(poxflag)
                

    def test_del_particle(self): 
        self.strucC1.del_particle(4)
        
    def tearDown(self):
        del self.strucC1 
        del self.strucC2
        
        

if __name__ == '__main__':
    unittest.main()
        
                