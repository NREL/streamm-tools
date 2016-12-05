import unittest

from streamm import parameters

class TestLJtype(unittest.TestCase):
    def setUp(self):
	self.ljtype_i = parameters.LJtype("C")
        self.ljtype_i.epsilon = 1.05
        self.ljtype_i.sigma = 3.25

    def test_LJtype(self):        
        self.assertEqual(str(self.ljtype_i) ,' LJ  C  epsilon 1.050000 sigma 3.250000 ')
        
    def tearDown(self):
        del self.ljtype_i 
        self.ljtype_i = None

class Testbondtype(unittest.TestCase):
    def setUp(self):
	self.bondtype_i = parameters.Bondtype("C","H")
        self.bondtype_i.r0 = 0.56
        self.bondtype_i.kb = 24.023

    def test_bondstr(self):
        bond_str = ' bond  C - H type harmonic \n  harmonic r_0 = 0.560000 K = 24.023000 lammps index 0  gromcas index 0  '
        self.assertEqual(str(self.bondtype_i),bond_str)
        
    def tearDown(self):
        del self.bondtype_i 
        self.bondtype_i = None

class Testangletype(unittest.TestCase):
    def setUp(self):
	self.angletype_i = parameters.Angletype("HC","CH","HC")
        self.angletype_i.theta0 = 120.0
        self.angletype_i.kb = 4.56

    def test_anglestr(self):
        angle_str = ' angle  HC - CH - HC type harmonic \n  harmonic theta_0 = 120.000000 K = 4.560000 lammps index 0  gromcas index 0  '
        self.assertEqual(str(self.angletype_i),angle_str)
        
    def tearDown(self):
        del self.angletype_i 
        self.angletype_i = None


class Testdihtypeharmonic(unittest.TestCase):
    
    def setUp(self):
	self.dihtype_i = parameters.Dihtype("HC","CH","CH","HC",type="harmonic")
        self.dihtype_i.d = 4.0
        self.dihtype_i.mult = 3.0
        self.dihtype_i.theat_s = 45.0
        self.dihtype_i.kb = 80.6
        

    def test_dihstr(self):
        dih_str = ' dihedral  HC - CH - CH - HC type harmonic \n  harmonic d = 4.000000 mult = 3.000000 K = 80.600000 theat_s = 45.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None


class Testdihtypemultiharmonic(unittest.TestCase):
    def setUp(self):
	self.dihtype_i = parameters.Dihtype("HC","CH","CH","HC",type="multiharmonic")
        self.dihtype_i.d = 4.0
        self.dihtype_i.mult = 3.0
        self.dihtype_i.theat_s = 45.0
        self.dihtype_i.kb = 80.6
        

    def test_dihstr(self):
        dih_str = ' dihedral  HC - CH - CH - HC type multiharmonic \n  harmonic d = 4.000000 mult = 3.000000 K = 80.600000 theat_s = 45.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None


class Testdihtypeopls(unittest.TestCase):
    def setUp(self):
	self.dihtype_i = parameters.Dihtype("HC","CH","CH","HC",type="opls")
        self.dihtype_i.setopls(14.0,1.0,45.0,100.0)

    def test_dihstropls(self):
        dih_str = ' dihedral  HC - CH - CH - HC type opls \n  k1 = 14.000000 k2 = 1.000000 k3 = 45.000000 k4 = 100.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)

    def test_dihstrrb(self):
        self.dihtype_i.type = "rb"
        dih_str = ' dihedral  HC - CH - CH - HC type rb \n  C0 = 30.500000  C1 = 60.500000 C2 = 179.000000 C3 = -90.000000 C4 = -400.000000  C5 = 0.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
                
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None

class Testdihtyperb(unittest.TestCase):
    def setUp(self):
	self.dihtype_i = parameters.Dihtype("HC","CH","CH","HC",type="rb")        
        self.dihtype_i.setrb(0.1,23.4,73.1,32.5,66.7,55.0)        

        
    def test_dihstrrb(self):
        dih_str = ' dihedral  HC - CH - CH - HC type rb \n  C0 = 0.100000  C1 = 23.400000 C2 = 73.100000 C3 = 32.500000 C4 = 66.700000  C5 = 55.000000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    def test_dihstropls(self):
        self.dihtype_i.type = "opls"
        dih_str = ' dihedral  HC - CH - CH - HC type opls \n  k1 = -95.550000 k2 = -139.800000 k3 = -16.250000 k4 = -16.675000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.dihtype_i),dih_str)
        
    def tearDown(self):
        del self.dihtype_i 
        self.dihtype_i = None
                         

class Testimproper(unittest.TestCase):
    
    def setUp(self):
	self.imptype_i = parameters.Imptype("C1","C2","C3","C4",type="improper")
        self.imptype_i.e0 = 180.0
        self.imptype_i.ke = 67.3
        self.imptype_i.pn = 4.0
        

    def test_impstr(self):
        imp_str = ' improper  C1 - C2 - C3 - C4 type improper \n  imp e0 = 180.000000 ke = 67.300000 lammps index 0  gromcas index 0 '
        self.assertEqual(str(self.imptype_i),imp_str)
        
    def tearDown(self):
        del self.imptype_i 
        self.imptype_i = None
                                             
class TestParameter(unittest.TestCase):
    def setUp(self):
	self.paramC = parameters.Container()

    def test_str(self):
        empty_paramC_str = '\n---------------------------------------------------------------------\n    Parameters \n---------------------------------------------------------------------\n      LJ parameters 0 \n      Bond parameters 0 \n      Angle parameters 0 \n      Dihedral parameters 0 \n      Imporper Dihedral parameters 0 \n'
        self.assertEqual(str(self.paramC) ,empty_paramC_str)


    def test_LJ(self):
        lj_str = []
	self.ljtype_i = parameters.LJtype("Ir")
        self.ljtype_i.epsilon = 2.35
        self.ljtype_i.sigma = 4.15
        self.paramC.add_LJtype(self.ljtype_i)
        lj_str.append(' LJ  Ir  epsilon 2.350000 sigma 4.150000 ')
	self.ljtype_i = parameters.LJtype("C")
        self.ljtype_i.epsilon = 1.05
        self.ljtype_i.sigma = 3.25
        self.paramC.add_LJtype(self.ljtype_i)
        lj_str.append(' LJ  C  epsilon 1.050000 sigma 3.250000 ')
	self.ljtype_i = parameters.LJtype("H")
        self.ljtype_i.epsilon = 0.75
        self.ljtype_i.sigma = 3.15
        self.paramC.add_LJtype(self.ljtype_i)
        lj_str.append(' LJ  H  epsilon 0.750000 sigma 3.150000 ')
        for ljtkey_i, ljtype_i  in self.paramC.ljtypes.iteritems():
            self.assertEqual(str(ljtype_i),lj_str[ljtkey_i])
        


    def test_bond(self):
        bond_str = []
	self.bondtype_i = parameters.Bondtype("Ir","C")
        self.bondtype_i.r0 = 1.02
        self.bondtype_i.kb = 13.563
        self.paramC.add_bondtype(self.bondtype_i)
        bond_str.append(' bond  Ir - C type harmonic \n  harmonic r_0 = 1.020000 K = 13.563000 lammps index 0  gromcas index 0  ')

	self.bondtype_i = parameters.Bondtype("C","C")
        self.bondtype_i.r0 = 0.56
        self.bondtype_i.kb = 24.023
        self.paramC.add_bondtype(self.bondtype_i)
        bond_str.append(' bond  C - C type harmonic \n  harmonic r_0 = 0.560000 K = 24.023000 lammps index 0  gromcas index 0  ')

	self.bondtype_i = parameters.Bondtype("C","H")
        self.bondtype_i.r0 = 0.43
        self.bondtype_i.kb = 65.123
        self.paramC.add_bondtype(self.bondtype_i)
        bond_str.append(' bond  C - H type harmonic \n  harmonic r_0 = 0.430000 K = 65.123000 lammps index 0  gromcas index 0  ')
        for btkey_i,bondtype_i  in self.paramC.bondtypes.iteritems():
            self.assertEqual(str(bondtype_i),bond_str[btkey_i])


    def test_angle(self):
        angle_str = []
	self.angletype_i = parameters.Angletype("H","C","H")
        self.angletype_i.theta0 = 120.0
        self.angletype_i.kb = 4.56
        angle_str.append(' angle  H - C - H type harmonic \n  harmonic theta_0 = 120.000000 K = 4.560000 lammps index 0  gromcas index 0  ')
        
	self.angletype_i = parameters.Angletype("C","Ir","C")
        self.angletype_i.theta0 = 90.0
        self.angletype_i.kb = 2.86
        angle_str.append(' angle  C - Ir - C type harmonic \n  harmonic theta_0 = 90.000000 K = 2.860000 lammps index 0  gromcas index 0  ')

	self.angletype_i = parameters.Angletype("Ir","C","H")
        self.angletype_i.theta0 = 120.0
        self.angletype_i.kb = 1.73
        angle_str.append(' angle  Ir - C - H type harmonic \n  harmonic theta_0 = 120.000000 K = 1.730000 lammps index 0  gromcas index 0  ')

        for atkey_i,angletype_i  in self.paramC.angletypes.iteritems():
            self.assertEqual(str(angletype_i),angle_str[atkey_i])


    def test_dih(self):
        dih_str = []
	self.dihtype_i = parameters.Dihtype("H","C","Ir","H",type="harmonic")
        self.dihtype_i.d = 4.0
        self.dihtype_i.mult = 3.0
        self.dihtype_i.theat_s = 45.0
        self.dihtype_i.kb = 80.6
        dih_str.append(' dihedral  HC - CH - CH - HC type harmonic \n  harmonic d = 4.000000 mult = 3.000000 K = 80.600000 theat_s = 45.000000 lammps index 0  gromcas index 0 ')
        for dtkey_i, dihtype_i  in self.paramC.dihtypes.iteritems():
            self.assertEqual(str(dihtype_i),dih_str[dtkey_i])
            
        
    def test_imp(self):
        imp_str = []
	self.imptype_i = parameters.Imptype("Ir","C","C","C",type="harmonic")
        self.imptype_i.e0 = 180.0
        self.imptype_i.ke = 67.3
        self.imptype_i.pn = 4.0
        imp_str.append(' improper  Ir - C - C - C type harmonic \n  imp e0 = 180.000000 ke = 67.300000 lammps index 0  gromcas index 0 ')
        for itkey_i, imptype_i  in self.paramC.imptypes.iteritems():    
            self.assertEqual(str(imptype_i),imp_str[itkey_i])
        
    def tearDown(self):
        del self.paramC 
        self.paramC = None
                        
if __name__ == '__main__':
    unittest.main()


