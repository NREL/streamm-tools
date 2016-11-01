import unittest

from streamm import periodictable 

class TestElement(unittest.TestCase):

    def test_n_elements(self):
        self.assertEqual(periodictable.n_elements,109)

    def test_elements(self):
        el_i = periodictable.elements['H']
        self.assertEqual(el_i['mass'],1.0080,"Hyrdogen mass set incorrectly")
        self.assertEqual(el_i['cov_radii'],0.230,"Hyrdogen cov_radii set incorrectly")
        self.assertEqual(el_i['vdw_radii'],1.090,"Hyrdogen vdw_radii set incorrectly")
        self.assertEqual(el_i['number'],1,"Hyrdogen atomic number set incorrectly")

        el_i = periodictable.elements['C']
        self.assertEqual(el_i['mass'],12.0110,"Carbon mass set incorrectly")
        self.assertEqual(el_i['cov_radii'],0.680,"Carbon cov_radii set incorrectly")
        self.assertEqual(el_i['vdw_radii'],1.700,"Carbon vdw_radii set incorrectly")
        self.assertEqual(el_i['number'],6,"Carbon atomic number set incorrectly")
        
    def test_element_symbols(self):
        el_i =  periodictable.element_symbol("S")
        self.assertEqual(el_i['number'],16)
        self.assertEqual(el_i['mass'],32.0660,"Sulfur mass set incorrectly")
        self.assertEqual(el_i['cov_radii'],1.020,"Sulfur cov_radii set incorrectly")
        self.assertEqual(el_i['vdw_radii'],1.800,"Sulfur vdw_radii set incorrectly")
                  
    def test_element_mass(self):
        el_i =  periodictable.element_mass(15.999)
        self.assertEqual(el_i['number'],8)
        self.assertEqual(el_i['mass'],15.999,"Oxygen mass set incorrectly")
        self.assertEqual(el_i['cov_radii'],0.68,"Oxygen cov_radii set incorrectly")
        self.assertEqual(el_i['vdw_radii'],1.52,"Oxygen vdw_radii set incorrectly")
                              
    def test_element_number(self):
        el_i =  periodictable.element_number(7)
        self.assertEqual(el_i['number'],7)
        self.assertEqual(el_i['mass'],14.007,"Nitrogen mass set incorrectly")
        self.assertEqual(el_i['cov_radii'],0.68,"Nitrogen cov_radii set incorrectly")
        self.assertEqual(el_i['vdw_radii'],1.55,"Nitrogen vdw_radii set incorrectly")
                                       
if __name__ == '__main__':
    unittest.main()

