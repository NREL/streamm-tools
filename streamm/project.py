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
        
        self.prefix = 'proj'        
        self.calculations = dict()
        self.resources = dict()
        
    def __del__(self):
        """
        Delete Calculation object
        """
        # Call base class destructor
        Calculation.__del__(self)        
        del self.calculations
        del self.resources


    def dump_json(self):
        '''
        Dump json file for reference 
        '''
        json_data = dict()
        json_data['meta'] =  self.meta
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
            if( calc_i.resource.meta['type'] == "local" )
                os.chdir(calc_i.dir['scratch'])
            calc_i.check()
            if( calc_i.resource.meta['type'] == "local" )
                os.chdir(calc_i.dir['home'])

    def run(self):
        '''
        Run calculations in the project 
        '''
        
        for calc_key,calc_i in self.calculations.iteritems():
            if( calc_i.resource.meta['type'] == "local" )
                os.chdir(calc_i.dir['scratch'])
            calc_i.run()
            if( calc_i.resource.meta['type'] == "local" )
                os.chdir(calc_i.dir['home'])
                        

    def store(self):
        '''
        Store calculations in the project 
        '''
        
        for calc_key,calc_i in self.calculations.iteritems():
            if( calc_i.resource.meta['type'] == "local" )
                os.chdir(calc_i.dir['scratch'])
            calc_i.store()
            if( calc_i.resource.meta['type'] == "local" )
                os.chdir(calc_i.dir['home'])
                                
