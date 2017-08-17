
class TestLattice(unittest.TestCase):
    def setUp(self):
        # fixture 
        matrix = [ 100,0,0,0,100,0,0,0,100 ]
        self.lat = structure.Lattice()
        self.lat.set_matrix(matrix)
        
        self.assertEqual(self.lat._lengths[0],100.0)
        self.assertEqual(self.lat._lengths[1],100.0)
        self.assertEqual(self.lat._lengths[2],100.0)

        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
        
    def test_changematrix(self):
        matrix = [ 132,0,0,0,127,0,0,0,150 ]
        self.lat.set_matrix(matrix)
        self.assertEqual(self.lat._lengths[0],132.0)
        self.assertEqual(self.lat._lengths[1],127.0)
        self.assertEqual(self.lat._lengths[2],150.0)

        self.assertEqual(self.lat._angles[0],90.0)
        self.assertEqual(self.lat._angles[1],90.0)
        self.assertEqual(self.lat._angles[2],90.0)
                
    def tearDown(self):
        del self.lat 
        self.lat = None
