import unittest, os , random 
import numpy as np

from streamm import resource

class TestInt(unittest.TestCase):
    def setUp(self):
        os.chdir(os.path.dirname(__file__))
        self.res_i = resource.Resource('res1',res_type = 'local')
        
    def test_makedir(self):
        self.res_i.make_dir()
        #self.assertEqual(,)

    def tearDown(self):
        del self.res_i
        

class TestInt2(unittest.TestCase):
    def setUp(self):
        self.res_i = resource.Resource('res1',home_dir=os.path.dirname(__file__),res_type = 'local')
        
    def test_makedir(self):
        self.res_i.make_dir()
        #self.assertEqual(,)

    def tearDown(self):
        del self.res_i
        

class TestReadWriteJSON(unittest.TestCase):
    def setUp(self):
        self.res_i = resource.Resource('res1',home_dir=os.path.dirname(__file__),res_type = 'local')

    def test_writejson(self):
        os.chdir(os.path.dirname(__file__))
        self.res_i.dump_json()

        self.res_j = resource.Resource("res2")
        self.res_j.load_json('res_res1.json')
        self.assertEqual(self.res_i.meta['type'],self.res_j.meta['type'])
        self.assertEqual(self.res_i.dir['home'],self.res_j.dir['home'])
        self.assertEqual(self.res_i.dir['scratch'],self.res_j.dir['scratch'])
        self.assertEqual(self.res_i.dir['scripts'],self.res_j.dir['scripts'])
        self.assertEqual(self.res_i.dir['templates'],self.res_j.dir['templates'])
        

    def tearDown(self):
        del self.res_i
                
