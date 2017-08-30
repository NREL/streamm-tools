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
    import streamm.buildingblock.container as container
    import streamm.buildingblock.bbatom as bbatom

except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    import container
    import  bbatom


class TestBuildThiophene(unittest.TestCase):

    def setUp(self):
        
        self.Th = container.Container('thiophene')
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
            pt_i = bbatom.BBatom(symbols[i])
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.Th.add_partpos(pt_i,pos_i)

        self.Th.particles[0].cplytag = 'term_'
        self.Th.particles[1].cplytag = 'func_'
        self.Th.particles[3].cplytag = 'term_'
        self.Th.particles[5].cplytag = 'termcap_'
        self.Th.particles[6].cplytag = 'funccap_'
        self.Th.particles[8].cplytag = 'termcap_'
        self.Th.parse_cplytag()
            
    def test_tag(self):
        self.assertEqual(self.Th.tag,'thiophene')
        #el_cnt = calc_elcnt

    def test_parsecplytag(self):
        self.assertEqual(self.Th.particles[0].bbid,'X')
        self.assertEqual(self.Th.particles[1].bbid,'X')
        self.assertEqual(self.Th.particles[2].bbid,'')
        self.assertEqual(self.Th.particles[3].bbid,'X')
        self.assertEqual(self.Th.particles[4].bbid,'')
        self.assertEqual(self.Th.particles[5].bbid,'T')
        self.assertEqual(self.Th.particles[6].bbid,'R')
        self.assertEqual(self.Th.particles[7].bbid,'')
        self.assertEqual(self.Th.particles[8].bbid,'T')
        
    def test_write_cply(self):
        self.Th.write_cply(write_ff=False)


    def tearDown(self):
        del self.Th

class TestBuildHexane(unittest.TestCase):

    def setUp(self):
        
        self.Hx = container.Container('hexane')
        
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
            pt_i = bbatom.BBatom(symbols[i])
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.Hx.add_partpos(pt_i,pos_i)

        self.Hx.particles[0].cplytag = 'rg_'
        self.Hx.particles[1].cplytag = 'rgcap_'
        self.Hx.parse_cplytag()


    def test_tag(self):
        self.assertEqual(self.Hx.tag,'hexane')
        #el_cnt = calc_elcnt

    def test_parsecplytag(self):
        self.assertEqual(self.Hx.particles[0].bbid,'rg_')
        self.assertEqual(self.Hx.particles[1].bbid,'R')

    def test_write_cply(self):
        self.Hx.write_cply(write_ff=False)
                
    def tearDown(self):
        del self.Hx


class TestContainer(unittest.TestCase):
    def setUp(self):
        self.bblockC = container.Container()
        self.bblockC2 = container.Container()
        
    def test_read_cply(self):
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene.cply')
        self.bblockC.read_cply(cply_file)
        self.assertEqual(str(self.bblockC)," thiophene")
        self.assertEqual(self.bblockC.n_term,2)
        self.assertEqual(self.bblockC.n_func,1)

        func_list = []
        func_list.append(' term( H 5)-( C 0) ')
        func_list.append(' term( H 8)-( C 3) ')
        func_list.append(' func( H 6)-( C 1) ')

        funcp_cnt = 0 
        for pkey_i in self.bblockC.terms:
            term_str = " term(%s%d)"%(str(self.bblockC.particles[pkey_i]),pkey_i)
            for pkey_j in self.bblockC.bonded_nblist.getnbs(pkey_i):
                term_str += "-(%s%d) "%(str(self.bblockC.particles[pkey_j]),pkey_j)
            self.assertEqual(term_str,func_list[funcp_cnt])
            funcp_cnt += 1 
        
        for pkey_i in self.bblockC.funcs:
            term_str = " func(%s%d)"%(str(self.bblockC.particles[pkey_i]),pkey_i)
            for pkey_j in self.bblockC.bonded_nblist.getnbs(pkey_i):
                term_str += "-(%s%d) "%(str(self.bblockC.particles[pkey_j]),pkey_j)
            self.assertEqual(term_str,func_list[funcp_cnt])
            funcp_cnt += 1 
            
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bblockC.tag)
        self.bblockC.write_xyz(file_i)


    def test_write_cply(self):
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene.cply')
        self.bblockC.read_cply(cply_file)
        self.bblockC.tag = 'pt_n1'
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bblockC.tag)
        self.bblockC.write_cply(file_i)


    def test_writeread_cply(self):
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene.cply')
        self.bblockC.read_cply(cply_file)
        self.bblockC.tag = 'pt_n1'
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bblockC.tag)
        self.bblockC.write_cply(file_i)
        self.bblockC2.read_cply(file_i)
        self.bblockC2.tag = 'pt_n1v2'
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bblockC2.tag)
        self.bblockC2.write_cply(file_i)

        
    def tearDown(self):
        del self.bblockC 
        self.bblockC = None


class TestJoin(unittest.TestCase):
    def setUp(self):
        self.bblockC_i = container.Container()
        self.bblockC_j = container.Container()
        self.bblockC_p3htnX = container.Container()
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene.cply')
        self.bblockC_i.read_cply(cply_file)        
        self.bblockC_i.write_cply('thiophene_v2.cply')        

        cply_file = os.path.join(os.path.dirname(__file__), 'hexane.cply')
        self.bblockC_j.read_cply(cply_file)
        
        self.th_v2 = container.Container()
        self.th_v2.tag = "testthiov2"
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene_v2.cply')
        self.th_v2.read_cply(cply_file)        
        
    def test_cat1(self):
        self.bblockC_p3htn1 = buildingblock.attach(self.bblockC_i,self.bblockC_j,"R",0,"R",0)
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bblockC_p3htn1.tag)
        self.bblockC_p3htn1.write_xyz(file_i)
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bblockC_p3htn1.tag)
        self.bblockC_p3htn1.write_cply(file_i)

    def test_cat2(self):

        self.bblockC_ptn2 = buildingblock.attach(self.bblockC_i,self.bblockC_i,"T",0,"T",1,tag = "P3HT_n1")
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bblockC_ptn2.tag)
        self.bblockC_ptn2.write_xyz(file_i)
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bblockC_ptn2.tag)
        self.bblockC_ptn2.write_cply(file_i)

        
    def test_cat3(self):
        self.bblockC_p3htn1 = buildingblock.attach(self.bblockC_i,self.bblockC_j,"R",0,"R",0,tag = "P3HT_nX")
        self.bblockC_p3htnX += self.bblockC_p3htn1 # Same as self.bblockC_p3htnX = copy.deepcopy( self.bblockC_p3htn1 )
        # for n in range(1,xN):
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"T",0,"T",1)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"T",0,"T",1)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"T",1,"T",0)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"T",0,"T",1)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_p3htn1,"T",0,"T",1)

        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_j,"T",0,"R",0)
        self.bblockC_p3htnX = buildingblock.attach(self.bblockC_p3htnX,self.bblockC_j,"T",0,"R",0)
        
        self.bblockC_p3htnX.tag = "P3HT_nX"
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bblockC_p3htnX.tag)
        self.bblockC_p3htnX.write_xyz(file_i)
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bblockC_p3htnX.tag)
        self.bblockC_p3htnX.write_cply(file_i)

    def test_cat4(self):
        self.th_v2.tag = "testthiov2"        
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.th_v2.tag)
        self.th_v2.write_xyz(file_i)
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.th_v2.tag)
        self.th_v2.write_cply(file_i)
        

    def tearDown(self):
        del self.bblockC_i 
        del self.bblockC_j
        del self.bblockC_p3htnX
        self.bblockC_i = None
        self.bblockC_j = None


class Testattachprep(unittest.TestCase):
    def setUp(self):

        # Read in thiophene building block 
        self.Th = container.Container('thiophene')
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene.cply')
        self.Th.read_cply(cply_file)
        
        self.Th2 = container.Container('Th2')
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene.cply')
        self.Th2.read_cply(cply_file)

        self.bb_i = self.Th2.prepattach("T",0,0,dir=-1)
        self.bb_i.tag += "_"

        angle_rad = 90.0*math.pi/180.0
        self.bb_j = self.Th.prepattach("T",1,0,dir=1,yangle=angle_rad)
        self.bb_j.tag += "_"

        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bb_i.tag)
        self.bb_i.write_xyz(file_i)
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bb_j.tag)
        self.bb_j.write_xyz(file_i)

        self.bbC_i,self.bbC_j =  buildingblock.shiftprep(self.bb_i,self.bb_j )
        self.overlap_found =  buildingblock.checkprep(self.bb_i,self.bb_j )

        self.Th2_Th =  buildingblock.attachprep(self.bbC_i,self.bbC_j )
        self.Th2_Th.tag = self.bb_i.tag + self.bb_j.tag + "v2"
        self.Th2_Th.properties['deptag'] =  self.bb_i.properties['deptag'] + self.bb_j.properties['deptag'] + "v2"
        self.Th2_Th.lat_cubic(100.0)
        
    def test_overlap_found(self):
        self.assertEqual(str(self.overlap_found),'False')
    
    def test_tag(self):
        self.assertEqual(str(self.Th2_Th.tag),'Th2_thiophene_v2')
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.Th2_Th.tag)
        self.Th2_Th.write_xyz(file_i)


        
    def test_dih1(self):
        n_bonds = self.Th2_Th.n_bonds
        bond_i = self.Th2_Th.bonds[n_bonds-1]
        self.key_i = bond_i.pkey1
        self.key_j = bond_i.pkey2

        self.assertEqual(self.key_i,0)
        self.assertEqual(self.key_j,11)


        self.dr_ij,self.mag_dr_ij = self.Th2_Th.lat.delta_pos(self.Th2_Th.positions[self.key_i],self.Th2_Th.positions[self.key_j])
        self.assertEqual(self.mag_dr_ij,1.36)
        
        nn_i =  self.Th2_Th.bonded_nblist.calc_nnab(self.key_i)
        nn_j =  self.Th2_Th.bonded_nblist.calc_nnab(self.key_j)

        self.assertEqual(nn_i,3)
        self.assertEqual(nn_j,3)
        
        # Find attached carbons 
        for self.key_k in self.Th2_Th.bonded_nblist.getnbs(self.key_i):
            if( self.key_k != self.key_j ):
                particle_k = self.Th2_Th.particles[self.key_k]
                if( particle_k.properties['number'] == 6 ): break
        for self.key_l in self.Th2_Th.bonded_nblist.getnbs(self.key_j):
            if( self.key_l != self.key_i ):
                particle_l = self.Th2_Th.particles[self.key_l]
                if( particle_l.properties['number'] == 6 ): break

            
        self.assertEqual(self.key_k,1)
        self.assertEqual(self.key_l,10)
        # Find angle
        dih_i = structure.Dihedral(self.key_k,self.key_i,self.key_j,self.key_l)
        
        self.cos_kijl = self.Th2_Th.calc_dihedral(dih_i)

        self.assertEqual(round(self.cos_kijl,6),0.0)
        
        
    def tearDown(self):
        del self.Th2
        del self.Th

        
class TestP3HT(unittest.TestCase):
    def setUp(self):

        # Read in thiophene building block 
        self.bb_thiophene = container.Container()
        cply_file = os.path.join(os.path.dirname(__file__), 'thiophene.cply')
        self.bb_thiophene.read_cply(cply_file)        

        # Read in hexane building block 
        self.bb_R_hexane = container.Container()
        cply_file = os.path.join(os.path.dirname(__file__), 'hexane.cply')
        self.bb_R_hexane.read_cply(cply_file)

        # Attach hexane to functional group point 0 on thiophene
        self.bb_p3ht = buildingblock.attach(self.bb_thiophene,self.bb_R_hexane,"R",0,"R",0,tag = "p3ht")

        # Create a xN membered chain
        xN = 5
        #self.bb_p3ht_nX = container.Container()
        #self.bb_p3ht_nX.tag = "p3ht_n%d"%(xN)
        self.bb_p3ht_nX = copy.deepcopy(self.bb_p3ht)

        for n in range(1,xN):
            # print " Adding %d "%(n)
            self.bb_p3ht_nX = buildingblock.attach(self.bb_p3ht_nX,self.bb_p3ht,"T",1,"T",0,tag="p3ht_n%d"%(n))

        # Print xyz file and cply file 
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bb_p3ht_nX.tag)
        self.bb_p3ht_nX.write_xyz(file_i)
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bb_p3ht_nX.tag)
        self.bb_p3ht_nX.write_cply(file_i)

    def test_add_struc(self):
        repX = 3 
        # Replicate repX times
        self.bb_p3ht_n10_x10 = container.Container()
        # self.bb_p3ht_n10_x10.tag= "%s_x%d"%(self.bb_p3ht_nX.tag,repX)
        seed = 256250

        matrix_i = self.bb_p3ht_n10_x10.lat._matrix
        matrix_i[0][0] = 300.0 
        matrix_i[1][1] = 300.0 
        matrix_i[2][2] = 300.0 
        self.bb_p3ht_n10_x10.lat.set_matrix(matrix_i)
        #
        name_i = "%s_x%d"%(self.bb_p3ht_nX.tag,repX)
        self.bb_p3ht_n10_x10 = self.bb_p3ht_n10_x10.add_struc(self.bb_p3ht_nX,repX,seed,tag=name_i,verbose=False)

        repX = 10
        name_i= "%s_%s_x%d"%(self.bb_p3ht_n10_x10.tag,self.bb_R_hexane.tag,repX)
        self.bb_p3ht_n10_x10 = self.bb_p3ht_n10_x10.add_struc_grid(self.bb_R_hexane,repX,tag=name_i,verbose=False)

        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.bb_p3ht_n10_x10.tag)
        self.bb_p3ht_n10_x10.write_xyz(file_i)
        file_i = os.path.join(os.path.dirname(__file__), "%s.cply"%self.bb_p3ht_n10_x10.tag)
        self.bb_p3ht_n10_x10.write_cply(file_i)


    def tearDown(self):
        del self.bb_thiophene 
        del self.bb_R_hexane

if __name__ == '__main__':
    unittest.main()
        