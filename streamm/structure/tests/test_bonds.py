
class TestBond(unittest.TestCase):
    def setUp(self):
	self.bond_i = structure.Bond(0,1)
        
    def test_str(self):
        bond_str = ' 0 - 1  '
        self.assertEqual(str(self.bond_i),bond_str)

    def tearDown(self):
        del self.bond_i 
        self.bond_i = None
