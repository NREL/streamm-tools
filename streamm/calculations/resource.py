# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Ph.D."
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"


"""
Class data structures for resource data
"""


import copy
import os
import json
import sys
import shutil

import logging
logger = logging.getLogger(__name__)

try:
    # Import pymatgen Class 
    import pymatgen_core.core.units as units 
except:
    raise ImportError("pymatgen import error for units object")

class Resource(object):
    '''
    Data structure for a compute resource


    Args:
        * tag (str) unique id for project
            
    '''

    def set_home(self,home_dir):
        '''
        Set the home/root directory of a resource 
        '''
        self.dir['home'] = home_dir
        self.dir['templates'] = '%s/templates'%(self.dir['home'])
        self.dir['scripts'] = '%s/scripts'%(self.dir['home'])
        self.dir['launch'] = '%s/scratch'%(self.dir['home'])
        self.dir['scratch'] = '%s/scratch'%(self.dir['home'])
        self.dir['storage'] = '%s/storage'%(self.dir['home'])
        self.dir['materials'] = '%s/materials'%(self.dir['home'])
        
    def __init__(self,tag="blank",home_dir='',res_type = 'local'):
        '''
        Constructor for a general Resource object.
        
        '''
        
        self.tag = str(tag)
        
        self.sufix = 'res'
        self.meta = dict()
        self.meta['type'] = res_type
        self.ssh = dict()
        self.ssh['username'] = 'user'
        self.ssh['address'] = ''
        # This is not really meta data
        if( len(home_dir) == 0 ):
            home_dir = os.getcwd()            
        self.dir = dict()
        self.set_home(home_dir)
        # These will be used for defaults for the simulation specs
        self.properties = dict()
        self.properties['walltime'] = 24
        self.properties['allocation'] = ''
        self.properties['nodes'] = int(1)
        self.properties['ppn'] = int(1)
        self.properties['nproc'] = self.properties['nodes']*self.properties['ppn'] 
        self.properties['pmem'] = 1500
        self.properties['queue'] = 'batch'
        self.properties['feature'] = '24core'
        self.properties['exe_command'] = './'
        self.properties['allocation'] = ''

    def __del__(self):
        '''
        Free memory 
        '''
        del self.sufix
        del self.tag 
        # Dictionaries:
        del self.meta
        del self.ssh
        del self.dir
        del self.properties

    def __str__(self):
        """
        Print resource information 
        """
        return str(self.meta)
    
    def make_dir(self):
        '''
        Check that needed directories exist 
        '''
        logger.debug("Creating directories for resource %s "%(self.tag))
        if( self.meta['type'] == "local" ):
            if ( not os.path.isdir(self.dir['home']) ):
                os.mkdir(self.dir['home'])                
            os.chdir(self.dir['home'])
            for dkey,dir_i in self.dir.iteritems():
                if ( not os.path.isdir(dir_i) ):
                    os.mkdir(dir_i)
                    
    def export_json(self,write_file=True):
        '''    
        Export particles to json
        
        Kwargs:
            * write_file (boolean) to dump json to a file
            
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        json_data = dict()
        json_data['meta'] = self.meta
        json_data['ssh'] = self.ssh
        json_data['dir'] = self.dir
        json_data['properties'] = self.properties
        
        # Write file 
        if( write_file ):
            file_name = "{}_{}.json".format(self.tag,self.sufix)
            logger.debug("Writting {}".format(file_name))
            with open(file_name,'wb') as fl:
                json.dump(json_data,fl)

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
            file_name = "{}_{}.json".format(self.tag,self.sufix)
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
        if( 'ssh' in json_data.keys() ):
            self.ssh = json_data['ssh']
        else:
            logger.warning('ssh not in json ')
        # 
        if( 'dir' in json_data.keys() ):
            self.dir = json_data['dir']
        else:
            logger.warning('meta not in json ')
        #dir
        if( 'properties' in json_data.keys() ):
            self.properties = json_data['properties']
        else:
            logger.warning('properties not in json ')
            
