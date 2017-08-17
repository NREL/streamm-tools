
class TestAngle(unittest.TestCase):
    def setUp(self):
	self.angle_i = structure.Angle(2,0,1)
        
    def test_str(self):
        angle_str = ' 2 - 0 - 1  '
        self.assertEqual(str(self.angle_i),angle_str)

    def tearDown(self):
        del self.angle_i 
        self.angle_i = None

