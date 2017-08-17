
class TestDihedral(unittest.TestCase):
    def setUp(self):
	self.dih_i = structure.Dihedral(2,0,1,5)
        
    def test_str(self):
        dih_str = ' 2 - 0 - 1 - 5 '
        self.assertEqual(str(self.dih_i),dih_str)

    def tearDown(self):
        del self.dih_i 
        self.dih_i = None


class TestImproper(unittest.TestCase):
    def setUp(self):
	self.imp_i = structure.Improper(0,1,2,3)
        
    def test_str(self):
        imp_str = ' 0 - 1 - 2 - 3 '
        self.assertEqual(str(self.imp_i),imp_str)

    def tearDown(self):
        del self.imp_i 
        self.imp_i = None
