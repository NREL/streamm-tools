# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"


import copy , os , json , sys
import time, datetime
from string import replace

import logging
logger = logging.getLogger(__name__)



try:
    # Import pymatgen Class 
    import pymatgen_core.core.units as units 
except:
    raise ImportError("pymatgen import error for units object")

# Import streamm dependencies 
from streamm.calculations.resource import Resource
from streamm.calculations.nwchem import NWChem
from streamm.calculations.gaussian import Gaussian
from streamm.calculations.lammps import LAMMPS

from streamm.calculations.resource import CalculationRes


class Project(CalculationRes):
    '''
    Data structure for a project
    
    
    Project is a set of calculations
    
    Projects track:
      file IO
      simulations ran
      computational resources used for each simulation
    
    Data for each action is stored in a json file

    Args:
        tag (str): String identifier for object
        
    Kwargs:
        units_conf (dict): Dictionary of units for each attribute type
        
    '''
    def __init__(self,tag,unit_conf=units.unit_conf ):
        # Base class constructor is called
                # Store the units of each attribute type 
        self.unit_conf = unit_conf

        CalculationRes.__init__(self,tag,unit_conf = unit_conf)
        self.meta['software'] = 'streamm_proj'
        self.prefix = 'proj'        
        self.calculations = dict()
        self.resources = dict()
        
    def __del__(self):
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
        
        self.files['data']['json'] = "%s_%s.json"%(self.prefix,self.tag)
        f = open(self.files['data']['json'], 'w')
        json.dump(json_data,f, indent=2)
        f.close()

        for calc_key,calc_i in self.calculations.iteritems():
            calc_i.dump_json()
            
    def load_json(self):
        '''
        Load json file for reference 
        '''        
        self.files['data']['json'] = "%s_%s.json"%(self.prefix,self.tag)
        
        try:
            with open(self.files['data']['json']) as f:            
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
                    calc_i = NWChem(calc_key)
                  elif( software_i == 'lammps' ):
                    calc_i = LAMMPS(calc_key)
                  elif( software_i == 'gaussian' ):
                    calc_i = Gaussian(calc_key)
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
            logger.warning(" File not found %s in %s "%(self.files['data']['json'],os.getcwd()))

    def add_calc(self,calc_i,deepcopy = False ):
        '''
        Add Calculation to the project
        
        Args:
            calc_i (resource.CalculationRes): calculation object
            
        Kwargs:
            deepcopy (boolean): whether to make a deepcopy of the calculation
            
        '''
        if( deepcopy ):
            self.calculations[calc_i.tag] = copy.deepcopy(calc_i)
        else:
            self.calculations[calc_i.tag] = calc_i
            
        return 

        
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

    def make_dir(self):
        '''
        Run calculations in the project 
        '''
        
        for calc_key,calc_i in self.calculations.iteritems():
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['home'])
            calc_i.make_dir()
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['home'])
                        
                  
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
                    
              
    def pull(self):
        '''
        Run calculations in the project 
        '''
        
        for calc_key,calc_i in self.calculations.iteritems():
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['scratch'])
            print os.getcwd()
            calc_i.pull()
            if( calc_i.resource.meta['type'] == "local" ):
                os.chdir(calc_i.dir['home'])
            print "Calculation %s has status %s"%(calc_i.tag,calc_i.meta['status'])
                                            

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
