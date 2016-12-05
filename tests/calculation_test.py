import unittest, os , random 
import numpy as np
import logging
logger = logging.getLogger(__name__)


from streamm import calculation
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
        self.calc = calculation.Calculation("calc001")
        
    def test_addstrucC(self):
        os.chdir(os.path.dirname(__file__))
        self.calc.add_strucC(self.bbTh)
        
    def tearDown(self):
        del self.calc
        del self.bbTh 
        

class TestFiles(unittest.TestCase):
    def setUp(self):

        self.calc = calculation.Calculation("calc001")
        
    def test_tag(self):
        
        self.assertEqual(str(self.calc.tag),"calc001")
        
    def test_addfile(self):
        file_type = 'input'
        file_key = 'cply'
        file_name = "thiophene.cply"
        self.calc.add_file(file_type,file_key,file_name)
        self.assertEqual(self.calc.files[file_type][file_key],file_name)
        file_type = 'output'
        file_key = 'log'
        file_name = "%s.log"%(self.calc.tag)
        self.calc.add_file(file_type,file_key,file_name)
        self.assertEqual(self.calc.files[file_type][file_key],file_name)
        file_type = 'data'
        file_key = 'groups'
        file_name = "groups_mol.csv"
        self.calc.add_file(file_type,file_key,file_name)
        self.assertEqual(self.calc.files[file_type][file_key],file_name)

        
    def tearDown(self):
        del self.calc
        

class TestReadWriteJSON(unittest.TestCase):
    def setUp(self):

        self.calc = calculation.Calculation("calc001")
        file_type = 'input'
        file_key = 'cply'
        file_name = "thiophene.cply"
        self.calc.add_file(file_type,file_key,file_name)
        self.assertEqual(self.calc.files[file_type][file_key],file_name)
        file_type = 'output'
        file_key = 'log'
        file_name = "%s.log"%(self.calc.tag)
        self.calc.add_file(file_type,file_key,file_name)
        self.assertEqual(self.calc.files[file_type][file_key],file_name)
        file_type = 'data'
        file_key = 'groups'
        file_name = "groups_mol.csv"
        self.calc.add_file(file_type,file_key,file_name)
        self.assertEqual(self.calc.files[file_type][file_key],file_name)
        
    def test_writejson(self):
        os.chdir(os.path.dirname(__file__))
        self.calc.dump_json()
        del self.calc
        self.calc = calculation.Calculation("calc001")
        self.calc.load_json()
      
    def tearDown(self):
        del self.calc

class TestThBL(unittest.TestCase):
    
    def setUp(self):
        self.calc_i = calculation.Calculation("ThBL")
        file_type = 'input'
        file_key = 'xyz'
        file_name = "thiophene.xyz"
        self.calc_i.add_file(file_type,file_key,file_name)
        os.chdir(os.path.dirname(__file__))
        self.struc_i = buildingblock.Container()
        input_file = self.calc_i.files['input']['xyz']
        self.struc_i.read_xyz(input_file)
        self.struc_i.lat_cubic(10.0)
        self.struc_i.bonded_nblist.guess_nblist(self.struc_i.lat,self.struc_i.particles,self.struc_i.positions,"cov_radii",radii_buffer=1.25)

    def test_bldist(self):
        # Set up log file 
        log_file = "%s.log"%self.calc_i.tag
        self.calc_i.add_file('output','log',log_file)
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        hdlr = logging.FileHandler(log_file,mode='w')
        hdlr.setLevel(logging.DEBUG)
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        prop = 'length'
        units = 'Angstroms'
        self.calc_i.units[prop] = units
        logger.info(" Setting prop %s to units %s "%(prop,units))
        
        self.struc_i.bonded_bonds()
        # Get C-C bond lengths
        list_i = []
        list_j = []
        for k,p in self.struc_i.particles.iteritems():
            if( p.properties['symbol'] == 'C' ):
                list_i.append(k)
                list_j.append(k)
        CC_bondkeys = self.struc_i.find_bonds(list_i,list_j)
        self.struc_i.calc_bonds(CC_bondkeys)
        CC_bondlengths = structure.prop_list('length',CC_bondkeys,self.struc_i.bonds)
        self.assertEqual(round(CC_bondlengths[0],n_dec),1.3772)
        self.assertEqual(round(CC_bondlengths[1],n_dec),1.4321)
        self.assertEqual(round(CC_bondlengths[2],n_dec),1.3772)
        outf = 'CCbonds.csv'
        self.calc_i.add_file('data','CC_bondlengths',outf)
        self.struc_i.write_bonds(outf,CC_bondkeys)
        # Get C-S bond lengths
        list_i = []
        list_j = []
        for k,p in self.struc_i.particles.iteritems():
                if( p.properties['symbol'] == 'C' ):
                        list_i.append(k)
                if( p.properties['symbol'] == 'S' ):
                        list_j.append(k)
        CS_bondkeys = self.struc_i.find_bonds(list_i,list_j)
        self.struc_i.calc_bonds(CS_bondkeys)
        CS_bondlengths = structure.prop_list('length',CS_bondkeys,self.struc_i.bonds)
        self.assertEqual(round(CS_bondlengths[0],n_dec),1.6722)
        self.assertEqual(round(CS_bondlengths[0],n_dec),1.6722)
        outf = 'CSbonds.csv'
        self.calc_i.add_file('data','CS_bondlengths',outf)
        self.struc_i.write_bonds(outf,CS_bondkeys)
        self.calc_i.dump_json()


    def tearDown(self):
        del self.struc_i
        del self.calc_i


class TestThCosAng(unittest.TestCase):
    
    def setUp(self):
        self.calc_i = calculation.Calculation("ThCosAng")
        file_type = 'input'
        file_key = 'xyz'
        file_name = "thiophene.xyz"
        self.calc_i.add_file(file_type,file_key,file_name)
        os.chdir(os.path.dirname(__file__))
        self.struc_i = buildingblock.Container()
        input_file = self.calc_i.files['input']['xyz']
        self.struc_i.read_xyz(input_file)
        self.struc_i.lat_cubic(10.0)
        self.struc_i.bonded_nblist.guess_nblist(self.struc_i.lat,self.struc_i.particles,self.struc_i.positions,"cov_radii",radii_buffer=1.25)

    def test_CosAngDist(self):
        # Set up log file 
        log_file = "%s.log"%self.calc_i.tag
        self.calc_i.add_file('output','log',log_file)
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        hdlr = logging.FileHandler(log_file,mode='w')
        hdlr.setLevel(logging.DEBUG)
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        prop = 'cosine'
        units = None
        self.calc_i.units[prop] = units
        logger.info(" Setting prop %s to units %s "%(prop,str(units)))
        
        
        self.struc_i.bonded_angles()
        # Get C-C-C angles
        list_k = []
        list_i = []
        list_j = []
        for k,p in self.struc_i.particles.iteritems():
                if( p.properties['symbol'] == 'C' ):
                        list_k.append(k)
                        list_i.append(k)
                        list_j.append(k)
        CCC_keys = self.struc_i.find_angles(list_k,list_i,list_j)
        self.struc_i.calc_angles(CCC_keys)
        CCC_cos = structure.prop_list('cosine',CCC_keys,self.struc_i.angles)
        self.assertEqual(round(CCC_cos[0],n_dec),0.3669)
        self.assertEqual(round(CCC_cos[1],n_dec),0.3669)
        outf = 'CCC_cos.csv'
        self.calc_i.add_file('data','CCC_cos',outf)
        self.struc_i.write_angles(outf,CCC_keys)
        # Get C-S bond lengths
        list_k = []
        list_i = []
        list_j = []
        for k,p in self.struc_i.particles.iteritems():
                if( p.properties['symbol'] == 'C' ):
                        list_k.append(k)
                        list_j.append(k)
                if( p.properties['symbol'] == 'S' ):
                        list_i.append(k)
        CSC_keys = self.struc_i.find_angles(list_k,list_i,list_j)
        self.struc_i.calc_angles(CSC_keys)
        CSC_cos = structure.prop_list('cosine',CSC_keys,self.struc_i.angles)
        self.assertEqual(round(CSC_cos[0],n_dec),0.0669)
        outf = 'CSC_cos.csv'
        self.calc_i.add_file('data','CSC_cos',outf)
        self.struc_i.write_angles(outf,CSC_keys)
        self.calc_i.dump_json()

    def tearDown(self):
        del self.struc_i
        del self.calc_i


class TestThDihAng(unittest.TestCase):
    
    def setUp(self):
        self.calc_i = calculation.Calculation("ThCosAng")
        file_type = 'input'
        file_key = 'xyz'
        file_name = "thiophene.xyz"
        self.calc_i.add_file(file_type,file_key,file_name)
        os.chdir(os.path.dirname(__file__))
        self.struc_i = buildingblock.Container()
        input_file = self.calc_i.files['input']['xyz']
        self.struc_i.read_xyz(input_file)
        self.struc_i.lat_cubic(10.0)
        self.struc_i.bonded_nblist.guess_nblist(self.struc_i.lat,self.struc_i.particles,self.struc_i.positions,"cov_radii",radii_buffer=1.25)

    def test_CosAngDist(self):
        # Set up log file 
        log_file = "%s.log"%self.calc_i.tag
        self.calc_i.add_file('output','log',log_file)
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        hdlr = logging.FileHandler(log_file,mode='w')
        hdlr.setLevel(logging.DEBUG)
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        prop = 'cosine'
        units = None
        self.calc_i.units[prop] = units
        logger.info(" Setting prop %s to units %s "%(prop,str(units)))
        
        
        self.struc_i.bonded_dih()
                     
        # Select some particles 
        list_k = []
        list_i = []
        list_j = []
        list_l = []

        for k,p in self.struc_i.particles.iteritems():
                if( p.properties['symbol'] == 'C' ):
                        list_k.append(k)
                        list_i.append(k)
                        list_j.append(k)
                        list_l.append(k)

        CCCC_keys = self.struc_i.find_dihedrals(list_k,list_i,list_j,list_l)
        self.struc_i.calc_dihedrals(CCCC_keys)
        CCCC_cos = structure.prop_list('cosine',CCCC_keys,self.struc_i.dihedrals)
        self.assertEqual(round(CCCC_cos[0],n_dec),1.0)
        outf = 'CCCC_cos.csv'
        self.calc_i.add_file('data','CCCC_cos',outf)
        self.struc_i.write_dihedrals(outf,CCCC_keys)
        self.calc_i.dump_json()

    def tearDown(self):
        del self.struc_i
        del self.calc_i

class TestThDist(unittest.TestCase):
    
    def setUp(self):
        os.chdir(os.path.dirname(__file__))
        self.calc_i = calculation.CalculationRes("ThDist")
        file_type = 'input'
        file_key = 'xyz'
        file_name = "thiophene.xyz"
        self.calc_i.add_file(file_type,file_key,file_name)
        self.res_i = resource.Resource('res1',res_type = 'local')
        self.res_i.make_dir()
        self.struc_i = buildingblock.Container()
        self.calc_i.set_resource(self.res_i)
        self.calc_i.make_dir()
        
        file_dir = os.path.dirname(__file__)
        file_type = 'input'
        file_key = 'xyz'
        file_name = "thiophene.xyz"
        from_dirkey = 'home'
        to_dirkey = 'scratch'
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        
    def test_compress(self):
        self.calc_i.compress_dirfiles('input','scratch')
        
    def test_calcDist(self):
        os.chdir(self.calc_i.dir['scratch'])
        input_file = self.calc_i.files['input']['xyz']
        self.struc_i.read_xyz(input_file)
        self.struc_i.lat_cubic(10.0)
        self.struc_i.bonded_nblist.guess_nblist(self.struc_i.lat,self.struc_i.particles,self.struc_i.positions,"cov_radii",radii_buffer=1.25)

        # Set up log file 
        log_file = "%s.log"%self.calc_i.tag
        self.calc_i.add_file('output','log',log_file)
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        hdlr = logging.FileHandler(log_file,mode='w')
        hdlr.setLevel(logging.DEBUG)
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        prop = 'cosine'
        units = None
        self.calc_i.units[prop] = units
        logger.info(" Setting prop %s to units %s "%(prop,str(units)))
        dihcalc_tag = 'all_cos'
        self.struc_i.bonded_dih()
        self.struc_i.calc_dihedrals()
        outf = '%s.csv'%(dihcalc_tag)
        self.calc_i.add_file('data',dihcalc_tag,outf)
        self.struc_i.write_dihedrals(outf)
        # Get output
        file_key = self.calc_i.properties['comp_key']
        from_dirkey = 'scratch'
        to_dirkey = 'storage'
        file_type = 'output'
        self.calc_i.compress_dirfiles(file_type,from_dirkey)
        file_name = self.calc_i.files[file_type][self.calc_i.properties['comp_key']]
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        file_type = 'data'
        self.calc_i.compress_dirfiles(file_type,from_dirkey)
        file_name = self.calc_i.files[file_type][self.calc_i.properties['comp_key']]
        self.calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
        # 
    def tearDown(self):
        del self.struc_i
        del self.calc_i
