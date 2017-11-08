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
from datetime import datetime
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
from streamm.calculations.calculation import Calculation

class Project():
    '''
    Data structure for a project
    
    Project is a set of calculations
    
    Projects track:
      * file IO
      * simulations ran
      * computational resources used for each simulation
    
    Data for each action is stored in a json file

    Args:
        * tag (str): String identifier for object
        
    Kwargs:
        * units_conf (dict): Dictionary of units for each attribute type
        
    '''
    def __init__(self,tag):
        
        self.suffix = 'proj'        
        self.tag = str(tag)
        
        self.meta = {}
        
        self.dir = {}
        self.calculations = dict()
        self.resources = dict()
                
        self.meta['software'] = 'streamm_proj'
        dt = datetime.fromtimestamp(time.time())        
        self.meta['date'] = dt.isoformat()
        
    def __del__(self):
        del self.tag
        del self.meta
        del self.calculations
        del self.resources
        del self.suffix 


    def add_calc(self,calc_i,deepcopy = False ):
        '''
        Add Calculation to the project
        
        Args:
            * calc_i (resource.Calculation): calculation object
            
        Kwargs:
            * deepcopy (boolean): whether to make a deepcopy of the calculation
            
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
        self.resources[resource_i.tag] = resource_i
        # Set simulation directories based on resource
        self.dir = copy.deepcopy(resource_i.dir)
        self.dir['scratch'] = resource_i.dir['home']
        self.dir['launch'] = resource_i.dir['home']



    def export_json(self,write_file=True):
        '''    
        Export particles to json
        
        Kwargs:
            * write_file (boolean) to dump json to a file
            
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        #
        json_data = dict()
        # 
        json_data['meta'] = self.meta
        # Save tags of reference calculation
        json_data['calculations'] = {}
        for ck,calc in self.calculations.iteritems():
            json_data['calculations'][ck] = calc.meta['software'] 
            calc.export_json()
            
        # Save tags of resouces calculation
        json_data['resources'] = self.resources.keys()
        for rk,res in self.resources.iteritems():
            res.export_json()
        #
        # Write file 
        if( write_file ):
            file_name = "{}_{}.json".format(self.tag,self.suffix)
            logger.debug("Writting {}".format(file_name))
            with open(file_name,'wb') as fl:
                json.dump(json_data,fl,indent = 2)
        #
        return json_data

    def import_json(self,json_data={},read_file=True):
        '''    
        Export object to json
        
        Kwargs:
            * json_lattice (dict) json representation of the object
            * read_file (boolean) to read json from a file
            
        '''
        
        # 
        if( read_file ):
            file_name = "{}_{}.json".format(self.tag,self.suffix)
            logger.debug("Reading {}".format(file_name))
            with open(file_name,'rb') as fl:
                json_data = json.load(fl)
        # 
        logger.debug("Set object properties based on json")
        #      
        if( 'meta' in json_data.keys() ):               
            self.meta = json_data['meta']
        else:
            logger.warning('meta not in json ')
        #
        if( 'dir' in json_data.keys() ):               
            self.dir = json_data['dir']
        else:
            logger.warning('dir not in json ')

        if( 'calculations' in json_data.keys() ):               
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
                calc_i = Calculation(calc_key)
              else:
                print "Unknow software %s will set as general calculation object "%(software_i)
                calc_i = calculation.Calculation(calc_key)
              # Load calculation 
              calc_i.import_json()
              self.calculations[calc_key] = calc_i
        else:
            logger.warning('calculations not in json ')
            
        if( 'resources' in json_data.keys() ):               
            for res_tag in json_data['resources']:
                res_i = Resource(res_tag)
                res_i.import_json()
                self.resources[res_tag] = res_i
        else:
            logger.warning('resources not in json ')
            
            
