"""
Project is a set of calcultions

Projects track:
  file IO
  simulations ran
  computational resources used for each simulation

Data for each action is stored in a json file

"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

import copy , os , json , sys
import time, datetime
from string import replace

import logging
logger = logging.getLogger(__name__)

# streamm
from resource import Resource 
import resource, calculation, gaussian,lammps,nwchem,gromacs
from calculation import CalculationRes



class Project(CalculationRes):
    '''
    Data structure for a project 
    '''
    def __init__(self,tag):
        """
        Constructor for a general project object.
        
        Args:
            tag (str): String identifier for object
        """
        # Base class constructor is called
        CalculationRes.__init__(self,tag)
        self.meta['software'] = 'streamm_proj'
        self.prefix = 'proj'        
        self.calculations = dict()
        self.resources = dict()
        
    def __del__(self):
        """
        Delete Calculation object
        """
        # Call base class destructor
        CalculationRes.__del__(self)        
        del self.calculations
        del self.resources


    def dump_json(self):
        '''
        Dump json file for reference 
        '''
        json_data = dict()
        json_data['meta'] = self.meta
        json_data['units'] = self.units
        json_data['files'] = self.files
        json_data['data'] = self.data
        json_data['properties'] = self.properties
        json_data['references'] = self.references
        json_data['dir'] = self.dir
        json_data['calculations'] =  dict()
        for calc_key,calc_i in self.calculations.iteritems():
            json_data['calculations'][calc_key] = calc_i.meta['software']
        json_data['resources'] = self.resources.keys()
        
        json_file = "%s_%s.json"%(self.prefix,self.tag)
        f = open(json_file, 'w')
        json.dump(json_data,f, indent=2)
        f.close()


    def load_json(self):
        '''
        Load json file for reference 
        '''        
        json_file = "%s_%s.json"%(self.prefix,self.tag)
        try:
            with open(json_file) as f:            
                json_data = json.load(f)
                f.close()

                self.meta = json_data['meta']
                self.units = json_data['units']
                self.files = json_data['files']
                self.data = json_data['data'] 
                self.properties = json_data['properties'] 
                self.references = json_data['references'] 
                self.dir = json_data['dir']
                
                for calc_key,software_i in json_data['calculations'].iteritems():
                  calc_key = str(calc_key)
                  logger.debug("Loading calculation %s using %s module "%(calc_key,software_i))
                  if( software_i == 'nwchem' ):
                    calc_i = nwchem.NWChem(calc_key)
                  elif( software_i == 'lammps' ):
                    calc_i = lammps.LAMMPS(calc_key)
                  elif( software_i == 'gaussian' ):
                    calc_i = gaussian.Gaussian(calc_key)
                  elif( software_i == 'gromacs' ):
                    calc_i = gromacs.GROMACS(calc_key)
                  elif( software_i == 'streamm_proj' ):
                    calc_i = Project(calc_key)
                  elif( software_i == 'streamm_calc' ):
                    calc_i = CalculationRes(calc_key)
                  else:
                    print "Unknow software %s will set as general calculation object "%(software_i)
                    calc_i = calculation.CalculationRes(calc_key)
                  # Load calculation
                  calc_i.load_json()
                  self.calculations[calc_key] = calc_i

        except IOError:
            logger.warning(" File not found %s in %s "%(json_file,os.getcwd()))


    def check(self):
        '''
        Check if calculations in the project have finished
        '''
        
        for calc_key,calc_i in self.calculations.iteritems():
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['scratch'])
            calc_i.check()
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['home'])
            print "Calculation %s has status %s"%(calc_i.tag,calc_i.meta['status'])

    def run(self):
        '''
        Run calculations in the project 
        '''
        
        for calc_key,calc_i in self.calculations.iteritems():
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['scratch'])
            print os.getcwd()
            calc_i.run()
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['home'])
                        

    def store(self):
        '''
        Store calculations in the project 
        '''
        
        for calc_key,calc_i in self.calculations.iteritems():
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['scratch'])
            calc_i.store()
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['home'])
                                
                 
    def set_resource(self,resource_i):
        '''
        Set resource for simulation 
        '''
        self.resource = resource_i
        self.meta['resource'] = resource_i.tag
        # Add resource properties to calculation properties
        self.properties.update(resource_i.properties)
        # Set simulation directories based on resource
        self.dir = copy.deepcopy(resource_i.dir)
        self.dir['scratch'] = resource_i.dir['home']
        self.dir['launch'] = resource_i.dir['home']
        self.properties['scratch'] = resource_i.dir['scratch'] 
