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

try:
    import streamm.structure.containers as containers
            
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
    
    import lattices
    import nblists
    import atoms
    import bonds 
    import angles 
    import dihedrals 
    import impropers     
    

class TestContainer(unittest.TestCase):
    def setUp(self):
        self.strucC = containers.Container("Test_structure1")
        
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
        self.part = atoms.Atom(symbol="Ir")
        self.strucC.add_particle(self.part)
        pos_add = []
        pos_i = [ 0.0,0.0,0.0]
        pos_add.append(pos_i)
        self.strucC.add_position(pos_i)        
        self.assertEqual(self.strucC.n_particles,1)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,0.0,1.5]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,2)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,1.50,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,3)


        self.part = atoms.Atom(symbol="C")
        pos_i = [ 1.50,0.0,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,0.0,-1.5]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,-1.50,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = atoms.Atom(symbol="C")
        pos_i = [ -1.50,0.0,0.0]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))


        self.part = atoms.Atom(symbol="H")
        pos_i = [ 0.0,0.750,1.75]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))
        
        self.part = atoms.Atom(symbol="H")
        pos_i = [ 0.750,0.0,1.75]   
        pos_add.append(pos_i)
        self.strucC.add_partpos(self.part,pos_i)
        self.assertEqual(self.strucC.n_particles,len(pos_add))

        self.part = atoms.Atom(symbol="H")
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
        self.bond_i = structure.Bond(0,1)
        bond_str.append(' 0 - 1  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,2)
        bond_str.append(' 0 - 2  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,3)
        bond_str.append(' 0 - 3  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,4)
        bond_str.append(' 0 - 4  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,5)
        bond_str.append(' 0 - 5  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,6)
        bond_str.append(' 0 - 6  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(1,7)
        bond_str.append(' 1 - 7  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(1,8)
        bond_str.append(' 1 - 8  ')
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(1,9)
        bond_str.append(' 1 - 9  ')
        self.strucC.add_bond(self.bond_i)
        for bkey_i, bond_i  in self.strucC.bonds.iteritems():
            self.assertEqual(str(bond_i),bond_str[bkey_i])

    def test_nblist(self):

        self.part = atoms.Atom(symbol="Ir")
        pos_i = [ 0.0,0.0,0.0]
        self.strucC.add_partpos(self.part,pos_i)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,0.0,1.5]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,1.50,0.0]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 1.50,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,0.0,-1.5]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,-1.50,0.0]   
        self.strucC.add_partpos(self.part,pos_i)


        self.part = atoms.Atom(symbol="C")
        pos_i = [ -1.50,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)


        self.part = atoms.Atom(symbol="H")
        pos_i = [ 0.0,0.750,1.75]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.part = atoms.Atom(symbol="H")
        pos_i = [ 0.750,0.0,1.75]   
        self.strucC.add_partpos(self.part,pos_i)

        self.part = atoms.Atom(symbol="H")
        pos_i = [- 0.53,-0.53,1.75]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.bond_i = structure.Bond(0,1)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,2)
        self.strucC.add_bond(self.bond_i)
        # 
        self.bond_i = structure.Bond(0,3)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,4)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,5)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(0,6)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(1,7)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(1,8)
        self.strucC.add_bond(self.bond_i)
        self.bond_i = structure.Bond(1,9)
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
        self.angle_i = structure.Angle(2,0,1)
        angle_str.append(' 2 - 0 - 1  ')
        self.strucC.add_angle(self.angle_i)
        self.angle_i = structure.Angle(3,0,4)
        angle_str.append(' 3 - 0 - 4  ')
        self.strucC.add_angle(self.angle_i)
        self.angle_i = structure.Angle(5,0,6)
        angle_str.append(' 5 - 0 - 6  ')
        self.strucC.add_angle(self.angle_i)
        for akey_i, angle_i  in self.strucC.angles.iteritems():
            self.assertEqual(str(angle_i),angle_str[akey_i])


    def test_dihedrals(self):
        dihedral_str = []        
        self.dihedral_i = structure.Dihedral(2,0,1,7)
        dihedral_str.append(' 2 - 0 - 1 - 7 ')
        self.strucC.add_dihedral(self.dihedral_i)
        self.dihedral_i = structure.Dihedral(2,0,1,8)
        dihedral_str.append(' 2 - 0 - 1 - 8 ')
        self.strucC.add_dihedral(self.dihedral_i)
        self.dihedral_i = structure.Dihedral(2,0,1,9)
        dihedral_str.append(' 2 - 0 - 1 - 9 ')
        self.strucC.add_dihedral(self.dihedral_i)
        for dkey_i, dihedral_i  in self.strucC.dihedrals.iteritems():
            self.assertEqual(str(dihedral_i),dihedral_str[dkey_i])


    def test_impropers(self):
        improper_str = []        
        self.improper_i = structure.Improper(0,1,2,3)
        improper_str.append(' 0 - 1 - 2 - 3 ')
        self.strucC.add_improper(self.improper_i)
        self.improper_i = structure.Improper(0,4,5,6)
        improper_str.append(' 0 - 4 - 5 - 6 ')
        self.strucC.add_improper(self.improper_i)
        for ikey_i, improper_i  in self.strucC.impropers.iteritems():
            self.assertEqual(str(improper_i),improper_str[ikey_i])


    def test_dump_pickle(self):
        self.strucC.dump_pickle()
        
    def test_write_coord(self):
        coord = self.strucC.write_coord()
        print coord 
        self.assertEqual(coord,'')
        
    def test_write_xyz_str(self):
        xyz_str = self.strucC.write_xyz_str()
        self.assertEqual(xyz_str,'')
                
    def test_write_xyz_list(self):
        list_i = [2,3,4,8]
        self.strucC.write_xyz_list(list_i)
        
    def test_read_write_xyz(self): 
        self.strucC.write_xyz()
        self.strucC.read_xyz()
        
    def test_write_list(self): 
        list_i = [2,3,4,8]
        self.strucC.write_list()
        
    def test_shift(self):
        pkey = 2
        vec = [34.32,13.1234,-218.3]
        self.strucC.shift(pkey,vec)
        nptu.assert_almost_equal(self.positions[pkey] ,[ 55. ,  0. , 30.])
        
    def test_shift_pos(self): 
        vec = [-3478.324,234798.234,-234.234]
        positions_correct = []
        for i in self.strucC.n_particles:
            positions_correct.append([0,0,0])
        self.strucC.shift_pos(vec)
        nptu.assert_almost_equal(self.positions ,positions_correct)
        
    def test_pbc_pos(self):
        positions_correct = []
        for i in self.strucC.n_particles:
            positions_correct.append([0,0,0])
        self.strucC.pbc_pos()
        nptu.assert_almost_equal(self.positions ,positions_correct)
        
    def test_lat_cubic(self):
        len_o = 256.16
        matrix_correct = []
        matrix_correct.append([len_o,0.0,0.0])
        matrix_correct.append([0.0,len_o,0.0])
        matrix_correct.append([0.0,0.0,len_o])
        self.strucC.lat_cubic(len_o)
        nptu.assert_almost_equal(self.strucC.lat._matrix,matrix_correct)
        
        
    def test_calc_mass(self): 
        self.strucC.calc_mass()
        self.assertEqual(self.strucC.mass,0.0)
 
    def test_calc_charge(self): 
        self.strucC.calc_charge()
        self.assertEqual(self.strucC.charge,0.0)
        
    def test_calc_volume(self): 
        self.strucC.calc_volume()
        self.assertEqual(self.strucC.volume,0.0)
        
    def test_calc_density(self): 
        self.strucC.calc_density()
        self.assertEqual(self.strucC.density,0.0)
        
    def test_calc_composition(self): 
        self.strucC.calc_composition()
        self.assertEqual(self.strucC.composition,0.0)
        
        
    def test_calc_formula(self): 
        self.strucC.calc_formula()
        self.assertEqual(self.strucC.chemicalformula,0.0)
        
    def test_calc_center_mass(self): 
        self.strucC.calc_center_mass()
        self.assertEqual(self.strucC.center_mass,0.0)
        
    def test_sum_prop(self):
        propkey = 'charge'
        pkey_i = 3
        pkeyp_j = 8
        self.strucC.sum_prop(propkey,pkey_i,pkeyp_j)
        self.assertEqual(self.strucC.particles[pkeyp_j].properties[propkey],0.0)
        
        
    def test_maxtags(self): 
        self.strucC.maxtags()
        self.assertEqual(self.strucC.mol_max,0.0)
        
    def test_mol_mult(self): 
        self.strucC.mol_mult()
        self.assertEqual(self.strucC.mol_multiplier,0.0)
        
        
    def test_n_molecules(self): 
        self.assertEqual(self.strucC.n_molecules(),1)
        
    def test_get_max(self): 
        self.strucC.get_max()
        self.assertEqual(self.strucC.max_mol,0.0)
        self.assertEqual(self.strucC.max_residue,0.0)
        self.assertEqual(self.strucC.max_qgroup,0.0)
        self.assertEqual(self.strucC.max_ring,0.0)
        
    def test_add_mol(self): 
        self.strucC.add_mol(max_ref_mol)
        self.assertEqual(self.strucC.max_ring,0.0)
        
        
    def test_shift_tag(self): 
        self.strucC.shift_tag('charge',0.2)
        self.assertEqual(self.strucC.particles[0].properties["charge"],0.2)
    
    def test_change_mass(self): 
        self.strucC.change_mass('C',13.02)
        self.assertEqual(self.strucC.particles[0].properties["mass"],13.02)
        
    def test_findbond_key(self): 
        b_indx = self.strucC.findbond_key(5,6)
        self.assertEqual(b_indx,3)
        
    def test_pairs(self):
        list_i = [0,2]
        list_j = [1,3]
        
        pairvalue_ij_correct =  np.zeros((2,2), dtype=np.float64)  
        
        pairvalue_ij = self.strucC.find_pairs(list_i,list_j)
        
        nptu.assert_almost_equal(pairvalue_ij ,pairvalue_ij_correct)
        
        bin_r,bin_r_nn,volumes = self.strucC.distbin_pairs(list_i,list_j,pairvalue_ij,0.1,5.0,0)

        nptu.assert_almost_equal(bin_r ,pairvalue_ij_correct)
        nptu.assert_almost_equal(bin_r_nn ,pairvalue_ij_correct)
        nptu.assert_almost_equal(volumes ,pairvalue_ij_correct)
        
    def test_get_pos(self): 
        list_i = [0,2]
        npos_i_correct = np.zeros((2,3), dtype=np.float64)  
        npos_i = self.strucC.get_pos()
        
        nptu.assert_almost_equal(npos_i,npos_i_correct)
        
        
    def test_propcompile_particles(self): 
        self.strucC.propcompile_particles()
        prop_particles_correct = {}
        
        self.assertDictEqual(self.strucC.prop_particles,prop_particles_correct)
        
    def test_write_particles(self): 
        self.strucC.write_particles()
        
    def test_proc_dihedrals(self): 
        keys = [0,1]
        self.strucC.proc_dihedrals(keys)
        dih_dic_correct = {}
        #
        self.assertDictEqual(self.dih_dic,dih_dic_correct)
        
        
        
    def test_shift(self):
        self.part = atoms.Atom(symbol="C")
        pos_i = [ 1.0,2.0,1.5]   
        self.strucC.add_partpos(self.part,pos_i)
        self.strucC.shift( 0, np.array([-0.5,-1.5,-1.0]) )

        pos_j = self.strucC.positions[0]
        self.assertEqual(pos_j[0],0.50)
        self.assertEqual(pos_j[1],0.50)
        self.assertEqual(pos_j[2],0.50)


    def test_shift_pos(self):
        
        self.part = atoms.Atom(symbol="C")
        pos_i = [ 1.0,2.0,1.5]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.part = atoms.Atom(symbol="C")
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
        
        self.strucC = containers.Container()

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)        
        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.750,0.750,0.750]   
        self.strucC.add_partpos(self.part,pos_i)
        self.part = atoms.Atom(symbol="C")
        pos_i = [ -0.750,-0.750,0.750]   
        self.strucC.add_partpos(self.part,pos_i)
        self.part = atoms.Atom(symbol="C")
        pos_i = [ -1.250,-1.250,0.0]   
        self.strucC.add_partpos(self.part,pos_i)

    def test_guessnb(self):
        # self.strucC.bonded_nblist.build_nblist(self.strucC.particles,self.strucC.bonds )
        str_nbs_list = []

        str_nbs_list.append(' 0 - 1 ')
        str_nbs_list.append(' 0 - 2 ')
        str_nbs_list.append(' 1 - 0 ')
        str_nbs_list.append(' 2 - 0 ')
        str_nbs_list.append(' 2 - 3 ')
        str_nbs_list.append(' 3 - 2 ')
        cnt = 0
        self.strucC_j.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
        
        for pkey_i, particle_i in self.strucC.particles.iteritems():
             for pkey_j in   self.strucC.bonded_nblist.getnbs(pkey_i):
                 # print " str_nbs_list.append(\' %d - %d \')"%(pkey_i,pkey_j)
                 self.assertEqual(str_nbs_list[cnt],' %d - %d '%(pkey_i,pkey_j) )
                 cnt += 1
        

        
    def tearDown(self):
        del self.strucC 
        self.strucC = None

class TestBuildThiophene(unittest.TestCase):

    def setUp(self):
        
        self.Th = structure.Container('thiophene')
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
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.Th.add_partpos(pt_i,pos_i)
            
    def test_tag(self):
        self.assertEqual(self.Th.tag,'thiophene')
        #el_cnt = calc_elcnt
        
    def test_write_xyz(self):
        os.chdir(os.path.dirname(__file__))
        self.Th.write_xyz()

    def tearDown(self):
        del self.Th 

class TestBuildEthane(unittest.TestCase):

    def setUp(self):
        
        self.Eth = structure.Container('ethane')
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
            self.Eth.add_partpos(pt_i,pos_i)

    def test_write_xyz(self):
        os.chdir(os.path.dirname(__file__))
        self.Eth.write_xyz()
                    
    def tearDown(self):
        del self.Eth 
            
class TestReadXYZ(unittest.TestCase):
    # 
    def setUp(self):
        self.strucC = structure.Container("thiophene")
        file_i = os.path.join(os.path.dirname(__file__), "thiophene.xyz")
        # file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.strucC.tag)
        self.strucC.read_xyz(file_i)
        # 
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
            particle_line = "  %s %f %f %f  "%(particle_i.properties["symbol"] ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]))
            self.assertEqual(part_xyz[pkey_i],particle_line)
            # print " part_xyz.append(\'",particle_line,"\')"
        #         
    def tearDown(self):
        del self.strucC 
    
class Testiadd(unittest.TestCase):
    
    def setUp(self):
        self.strucC = containers.Container()
        self.strucC_j = containers.Container()

        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.0,0.0,0.0]   
        self.strucC.add_partpos(self.part,pos_i)
        
        self.part = atoms.Atom(symbol="C")
        pos_i = [ 0.750,0.750,0.750]   
        self.strucC.add_partpos(self.part,pos_i)
        # 
        self.strucC_j.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
        self.strucC.bonded_bonds()
        
        self.part = atoms.Atom(symbol="C")
        pos_i = [ -0.750,-0.750,0.750]   
        self.strucC_j.add_partpos(self.part,pos_i)

        self.part = atoms.Atom(symbol="C")
        pos_i = [ -1.250,-1.250,0.0]   
        self.strucC_j.add_partpos(self.part,pos_i)

        self.strucC_j.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
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
            self.assertEqual(pos_str_list[pkey_i],' %d - %s %f %f %f '%(pkey_i,particle_i.properties["symbol"],self.strucC.positions[pkey_i][0],self.strucC.positions[pkey_i][1],self.strucC.positions[pkey_i][2]))
            for pkey_j in   self.strucC.bonded_nblist.getnbs(pkey_i):
                #print " str_nbs_list.append(\' %d - %d \')"%(pkey_i,pkey_j)
                self.assertEqual(str_nbs_list[cnt],' %d - %d '%(pkey_i,pkey_j) )
                cnt += 1

    def tearDown(self):
        del self.strucC 
        del self.strucC_j
        self.strucC = None
        self.strucC_j = None

class TestProperties(unittest.TestCase):
    # 
    def setUp(self):
        self.strucC = containers.Container()
        self.strucC.tag = "thiophene"
        file_i = os.path.join(os.path.dirname(__file__), "thiophene.xyz")
        # file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.strucC.tag)
        self.strucC.read_xyz(file_i)
        matrix_i = self.strucC.lat._matrix
        matrix_i[0][0] = 100.0 
        matrix_i[1][1] = 100.0 
        matrix_i[2][2] = 100.0 
        self.strucC.lat.set_matrix(matrix_i)
        self.strucC_j.bonded_nblist = self.strucC.guess_nblist(0,radii_buffer=1.25)
        #
    def test_calc(self):
        #
        self.strucC.calc_mass()
        self.assertEqual(self.strucC.mass,84.14199999999998)
        self.strucC.calc_volume()
        self.assertEqual(self.strucC.volume,1000000.0)
        #
        self.strucC.calc_center_mass()
        self.assertEqual(self.strucC.center_mass[0],-1.1285902921242665)
        self.assertEqual(self.strucC.center_mass[1],-0.76661736952485104)
        self.assertEqual(self.strucC.center_mass[2],-0.00038730025433196259)
        #
        self.strucC.calc_density()
        self.assertEqual(str(self.strucC.density),'8.4142e-05')
        #
        self.strucC.calc_composition()
        el_n_list = [ n for n in self.strucC.composition  if  n > 0  ]
        self.assertEqual(str(el_n_list),'[4, 4, 1]')
        #
        self.strucC.calc_formula()
        self.assertEqual(str(self.strucC.chemicalformula),'C4H4S1')
        #
        self.strucC.calc_center_mass()
        self.assertEqual(str(self.strucC.center_mass),'[ -1.12859029e+00  -7.66617370e-01  -3.87300254e-04]')
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
        
    def test_guess_oplsa(self):
        self.strucC.guess_oplsa()
        self.assertEqual(self.strucC.particles[0].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[1].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[2].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[3].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[4].properties['fftype'],'S')
        self.assertEqual(self.strucC.particles[5].properties['fftype'],'HA')
        self.assertEqual(self.strucC.particles[6].properties['fftype'],'HA')
        self.assertEqual(self.strucC.particles[7].properties['fftype'],'HA')
        self.assertEqual(self.strucC.particles[8].properties['fftype'],'HA')
        

    def test_bonds(self):
        self.strucC.bonded_bonds()
        self.strucC.calc_bonds()
        self.strucC.write_bonds('bonds_all.csv')
        self.assertEqual(round(self.strucC.bonds[0].properties['length'],6),1.377226)
        self.assertEqual(round(self.strucC.bonds[1].properties['length'],6),1.672168)
        self.assertEqual(round(self.strucC.bonds[2].properties['length'],6),1.088203)
        self.assertEqual(round(self.strucC.bonds[3].properties['length'],6),1.432118)
        self.assertEqual(round(self.strucC.bonds[4].properties['length'],6),1.091483)
        self.assertEqual(round(self.strucC.bonds[5].properties['length'],6),1.377194)

    def test_angles(self):
        self.strucC.bonded_angles()
        self.strucC.calc_angles()
        self.strucC.write_angles('angles_all.csv')
        self.assertEqual(round(self.strucC.angles[0].properties['cosine'],6),0.367484)
        self.assertEqual(round(self.strucC.angles[12].properties['cosine'],6),0.066883)

    def test_dih(self):
        self.strucC.bonded_dih()
        self.strucC.calc_dihedrals()
        self.strucC.write_dihedrals('dih_all.csv')
        self.assertEqual(round(self.strucC.dihedrals[8].properties['cosine'],6),-1.0)

    def test_expand_matrix(self):
        exlat_frac = 0.05  # 5%
        self.strucC.lat.expand_matrix(exlat_frac)
        self.assertEqual(self.strucC.lat._lengths[0],105.0)
        self.assertEqual(self.strucC.lat._lengths[1],105.0)
        self.assertEqual(self.strucC.lat._lengths[2],105.0)

        
    def tearDown(self):
        del self.strucC 



class TestDihcalc(unittest.TestCase):
    # 
    def setUp(self):
        self.strucC = structure.Container("thiophene")
        file_i = os.path.join(os.path.dirname(__file__), "thiophene.xyz")
        # file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.strucC.tag)
        self.strucC.read_xyz(file_i)
        matrix_i = self.strucC.lat._matrix
        matrix_i[0][0] = 100.0 
        matrix_i[1][1] = 100.0 
        matrix_i[2][2] = 100.0 
        self.strucC.lat.set_matrix(matrix_i)
        #
    def test_calcdihc(self):
        dih_i = structure.Dihedral(3,2,1,0)
        self.strucC.calc_dihedral(dih_i)
        self.assertEqual(dih_i.properties['cosine'],0.99999980013324519)

    def tearDown(self):
        del self.strucC 


class TestProximityCheck(unittest.TestCase):
    # 
    def setUp(self):
        self.strucC1 = containers.Container()
        self.strucC1.tag = "thiophene"
        file_i = os.path.join(os.path.dirname(__file__), "thiophene.xyz")
        # file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.strucC1.tag)
        self.strucC1.read_xyz(file_i)
        matrix_i = self.strucC1.lat._matrix
        matrix_i[0][0] = 1000.0 
        matrix_i[1][1] = 1000.0 
        matrix_i[2][2] = 1000.0 
        self.strucC1.lat.set_matrix(matrix_i)
        
        self.strucC2 = containers.Container()
        self.strucC2.tag = "thiophene"
        file_i = os.path.join(os.path.dirname(__file__), "thiophene.xyz")
        # file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.strucC2.tag)
        self.strucC2.read_xyz(file_i)
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
        self.strucC.del_particle(4)
        
    def tearDown(self):
        del self.strucC1 
        
        
class TestGroupsProps(unittest.TestCase):
    # 
    def setUp(self):
        self.th = structure.Container("thiophene")
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.th.tag)
        self.th.read_xyz(file_i)
        self.th.lat_cubic(100.0)
        #
        for pkey_i, particle_i  in self.th.particles.iteritems():
            if( particle_i.properties['symbol'] == 'C' ):
                particle_i.properties['resname'] = "SCP2"
                particle_i.properties['residue'] = 1
            if( particle_i.properties['symbol'] == 'S' ):
                particle_i.properties['resname'] = "ThS"
                particle_i.properties['residue'] = 2
            if( particle_i.properties['symbol'] == 'H' ):
                particle_i.properties['resname'] = "HA"
                particle_i.properties['residue'] = 3
        self.strucC = containers.Container()
        self.strucC.lat_cubic(100.0)
        seed = 82343
        self.strucC = self.strucC.add_struc(self.th,10,seed,verbose=False)

    def test_molnumbers(self):
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            mol_i = int(pkey_i/self.th.n_particles) 
            self.assertEqual(str(particle_i.properties['mol']),str(mol_i))

    def test_groupmol(self):
        group_tag = 'mol'
        self.strucC.group_prop('mol',group_tag)
        groupset_i = self.strucC.groupsets[group_tag]
        self.assertEqual(str(len(groupset_i.groups)),str(10))
        groupset_i.calc_cent_mass()
        groupset_i.calc_cent_mass()
        
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

        groupset_i.calc_radius()
        groupset_i.calc_radius()
        groupset_i.calc_radius()
        for gkey,group_i in groupset_i.groups.iteritems():
            self.assertEqual(str(group_i.properties['cent_mass']),str(cm[gkey]))
            self.assertEqual(str(group_i.properties['radius']),'2.57775210944')
            self.assertEqual(str(group_i.properties['r_gy_sq']),'3.6041389371')
            # print "r_gy.append(\'%s\')"%str(group_i.properties)
        for gkey in groupset_i.properties['keys']:
            self.assertEqual(str(groupset_i.properties['cent_mass'][gkey]),str(cm[gkey]))
            self.assertEqual(str(groupset_i.properties['radius'][gkey]),'2.57775210944')
            self.assertEqual(str(groupset_i.properties['r_gy_sq'][gkey]),'3.6041389371')
            
        groupset_i.group_pbcs()

        os.chdir(os.path.dirname(__file__))
        groupset_i.write_cm_xyz()
        groupset_i.write_xyzs()
        groupset_i.dump_json()

        
    def test_groupres(self):
        group_tag = 'residue'
        self.strucC.group_prop('residue',group_tag)
        groupset_i = self.strucC.groupsets[group_tag]
        self.assertEqual(str(len(groupset_i.groups)),str(30))
        
        groupset_i.calc_cent_mass()
        groupset_i.calc_radius()

        #for gkey,group_i in groupset_i.groups.iteritems():
        self.assertEqual(round(groupset_i.properties['radius'][2],6),2.587885)
        self.assertEqual(round(groupset_i.properties['r_gy_sq'][2],6),4.967159)
        self.assertEqual(round(groupset_i.properties['Q_mn'][2][0][0],6),0.005067)
        self.assertEqual(round(groupset_i.properties['Rgy_eignval'][0][0],6),1.002185)
        self.assertEqual(round(groupset_i.properties['Rgy_eignval'][0][1],6),0.410354)
        self.assertEqual(round(groupset_i.properties['Rgy_eignval'][0][2],6),0.0)
        self.assertEqual(round(groupset_i.properties['A_sphere'][0],6),0.381661)
        self.assertEqual(round(groupset_i.properties['A_sphere_num'][0],6),1.52303)
        self.assertEqual(round(groupset_i.properties['A_sphere_dem'][0],6),1.995267)
        
        groupset_i.calc_dl()

        self.assertEqual(round(groupset_i.properties['dl_sq'][0],6),5.966521)
        self.assertEqual(round(groupset_i.properties['dl_sq'][2],6),20.729214)
                
        os.chdir(os.path.dirname(__file__))
        groupset_i.write_cm_xyz()
        groupset_i.write_xyzs()
        groupset_i.dump_json()
        

                
    def tearDown(self):
        del self.th         
        del self.strucC         

if __name__ == '__main__':
    unittest.main()
        
                