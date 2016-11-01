import unittest, os , random 
import numpy as np
import logging
logger = logging.getLogger(__name__)


from streamm import gaussian
from streamm import structure
from streamm import buildingblock
from streamm import parameters
from streamm import periodictable
from streamm import resource

n_dec = 4


class TestAddStruc(unittest.TestCase):
    def setUp(self):

        self.bbTh = buildingblock.Container('thiophene')
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
            pt_i = buildingblock.BBatom(symbols[i])
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.bbTh.add_partpos(pt_i,pos_i)
        self.bbTh.particles[0].properties['cplytag'] = 'term_'
        self.bbTh.particles[1].properties['cplytag'] = 'func_'
        self.bbTh.particles[3].properties['cplytag'] = 'term_'
        self.bbTh.particles[5].properties['cplytag'] = 'termcap_'
        self.bbTh.particles[6].properties['cplytag'] = 'funccap_'
        self.bbTh.particles[8].properties['cplytag'] = 'termcap_'
        self.bbTh.parse_cplytag()
        self.bbTh.bonded_nblist.guess_nblist(self.bbTh.lat,self.bbTh.particles,self.bbTh.positions,"cov_radii",radii_buffer=1.25)        
        self.bbTh.lat_cubic(10.0)
        os.chdir(os.path.dirname(__file__))
        self.bbTh.write_cply()
        self.bbTh.write_xyz()
        self.calc = gaussian.Gaussian("g001")
        
    def test_addstrucC(self):
        os.chdir(os.path.dirname(__file__))
        self.calc.add_strucC(self.bbTh)
        
    def tearDown(self):
        del self.calc
        del self.bbTh

        
class TestAddRes(unittest.TestCase):
    def setUp(self):
        os.chdir(os.path.dirname(__file__))
        self.calc_i = gaussian.Gaussian("g001")
        self.res_loc = resource.Resource('local',res_type = 'local')
        
    def test_addres(self):
        self.calc_i.set_resource(self.res_loc)
        self.res_loc.dump_json()
        self.calc_i.dump_json()
        
    def tearDown(self):
        del self.calc_i
        del self.res_loc
        
        
class TestcpTemplates(unittest.TestCase):
    
    def setUp(self):
        os.chdir(os.path.dirname(__file__))
        self.calc_i = gaussian.Gaussian("gTh")
        self.res_loc = resource.Resource('local',res_type = 'local')
        self.res_loc.load_json()
        self.res_loc.dir['root_templates'] = '../../templates/'
        self.calc_i.set_resource(self.res_loc)
        
        file_type = 'input'
        file_key = 'cply'
        file_name = "thiophene.cply"
        from_dirkey = 'home'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
    def test_cptemplate(self):

        file_type = 'input'
        file_key = 'com_template'
        file_name = "gaussian.com"
        from_dirkey = 'root_templates'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

    def test_cptemplate2(self):

        file_type = 'input'
        file_key = 'run_template'
        file_name = "gaussian.pbs"
        from_dirkey = 'root_templates'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
                                
    def test_check(self):
        print "check"


    def tearDown(self):
        del self.calc_i
        del self.res_loc
        

class TestWriteTemplates(unittest.TestCase):
    
    def setUp(self):
        self.calc_i = gaussian.Gaussian("gTh")
        self.res_loc = resource.Resource('local',res_type = 'local')
        self.res_loc.load_json()
        self.res_loc.dir['root_templates'] = '../../templates/'
        self.calc_i.set_resource(self.res_loc)
        
        file_type = 'input'
        file_key = 'cply'
        file_name = "thiophene.cply"
        from_dirkey = 'home'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
        file_type = 'templates'
        file_key = 'com'
        file_name = "gaussian.com"
        from_dirkey = 'root_templates'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
        file_type = 'templates'
        file_key = 'run'
        file_name = "gaussian.pbs"
        from_dirkey = 'root_templates'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
    def test_writecom(self):
        os.chdir(self.calc_i.dir['scratch'])
        input_file = self.calc_i.files['input']['cply']
        self.calc_i.strucC.read_cply(input_file)

        self.calc_i.load_str('templates','com')
        self.calc_i.load_str('templates','run')
        # Set properties
        self.calc_i.properties['commands'] = 'SP'
        self.calc_i.properties['charge'] = 0
        self.calc_i.properties['spin_mult'] = 1
        self.calc_i.properties['coord'] = self.calc_i.strucC.write_coord()

        self.calc_i.properties['charge'] = 1
        self.calc_i.properties['spin_mult'] = 0
        com_str = self.calc_i.replace_prop('com')
        self.calc_i.replacewrite_prop('com','input','com','%s.com'%(self.calc_i.tag))
        self.calc_i.replacewrite_prop('com','input','com2','%sv2.com'%(self.calc_i.tag))

        self.calc_i.properties['scratch_dir'] = self.calc_i.dir['scratch']
        self.calc_i.properties['input_com'] = self.calc_i.files['input']['com']
        self.calc_i.properties['input_com'] = self.calc_i.files['input']['com']
        
        self.calc_i.replacewrite_prop('run','scripts','run','%s.pbs'%(self.calc_i.tag))

        
        from_dirkey = 'scratch'
        file_type = 'input'
        self.calc_i.compress_dirfiles(file_type,from_dirkey)
        
        os.chdir(self.calc_i.dir['home'])
        self.calc_i.dump_json()
        

    def tearDown(self):
        del self.calc_i
        del self.res_loc
        
