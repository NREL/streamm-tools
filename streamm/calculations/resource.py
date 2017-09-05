# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

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



from streamm.calculations.calculation import Calculation

class Resource(object):
    '''
    Data structure for a compute resource  
    '''
    def __init__(self,tag="blank",home_dir='',res_type = 'local'):
        '''
        Constructor for a general Resource object.
        
        Input:
            tag (str) unique id for project
            
        '''
        
        self.tag = str(tag)
        
        self.prefix = 'res'
        self.meta = dict()
        self.meta['type'] = res_type
        self.ssh = dict()
        self.ssh['username'] = 'user'
        self.ssh['address'] = ''
        # This is not really meta data
        if( len(home_dir) == 0 ):
            home_dir = os.getcwd()            
        self.dir = dict()
        self.dir['home'] = home_dir
        self.dir['templates'] = '%s/templates'%(self.dir['home'])
        self.dir['scripts'] = '%s/scripts'%(self.dir['home'])
        self.dir['launch'] = '%s/scratch'%(self.dir['home'])
        self.dir['scratch'] = '%s/scratch'%(self.dir['home'])
        self.dir['storage'] = '%s/storage'%(self.dir['home'])
        self.dir['materials'] = '%s/materials'%(self.dir['home'])
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
        del self.prefix
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
            os.chdir(self.dir['home'])
            for dkey,dir_i in self.dir.iteritems():
                if ( not os.path.isdir(dir_i) ):
                    os.mkdir(dir_i)
                    
    def dump_json(self):
        '''
        Dump json file for reference 
        '''
        res_data = dict()
        res_data['meta'] = self.meta
        res_data['ssh'] = self.ssh
        res_data['dir'] = self.dir
        res_data['properties'] = self.properties
        
        json_file = "%s_%s.json"%(self.prefix,self.tag)
        f = open(json_file, 'w')
        json.dump(res_data,f, indent=2)
        f.close()


    def load_json(self,json_file=''):
        '''
        Load json file for reference 
        '''
        if( len(json_file) == 0 ):
            json_file = "%s_%s.json"%(self.prefix,self.tag)

        logger.info("Reading in resource json file %s "%(json_file))
        
        try:
            with open(json_file) as f:            
                json_data = json.load(f)
                f.close()
                
                self.meta = json_data['meta']
                self.ssh = json_data['ssh']
                self.dir = json_data['dir']
                self.properties = json_data['properties']

        except IOError:
            logger.warning(" File not found %s in %s "%(json_file,os.getcwd()))



class CalculationRes(Calculation):
    '''
    Derived type of Calculation for running a calculation on resource.
    In that files will be compressed and moved to scratch directories to run,
    then output files and output data will be compressed and moved to storage.

    '''
    def __init__(self, tag ):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """
        # Base class constructor is called
        Calculation.__init__(self, tag)
        #super(Calculation,self).__init__()
        # Computational Resource used for simulation/calculation  
        self.resource = Resource()
        self.dir = dict()
        #
    def __del__(self):
        """
        Destructor, clears object memory
        """
        # Call base class destructor
        Calculation.__del__(self)
        del self.resource
        del self.dir 
        #


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
        json_data['properties']['run_list'] = ''
        
        json_data['references'] = dict()
        for ref_key,ref_calc in self.references.iteritems(): 
            json_data['references'][ref_key] = ref_calc.tag
            
        json_data['dir'] = self.dir

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
                self.dir = json_data['dir']

                # Load resource 
                try:
                    res_tag = json_data['meta']['resource']
                except:
                    res_tag = ''
                    logger.warning("No resource found ")
                if( len(res_tag) > 0 ):
                    logger.debug("Resource tag found %s "%(res_tag))
                    resource_i = resource.Resource(str(res_tag))
                    resource_i.load_json()
                    self.resource = resource_i
                                    
                # Load references 
                try:
                    ref_tags = json_data['references']
                    for rekey,ref_tag in ref_tags:
                        ref_i = Calculation(ref_tag)
                        ref_i.load_json()
                        self.add_refcalc(ref_i)
                        logger.debug(" Need to set reference calculation type ")
                except:
                    logger.warning("No references found ")
                    
        except IOError:
            error_msg = " File not found %s in %s "%(json_file,os.getcwd())
            error_msg +=  " File not found %s in %s "%(json_file,os.getcwd())
            logger.warning(error_msg)

        
    def make_dir(self):
        '''
        Check that needed directories exist 
        '''
        logger.debug("Creating directories for resource %s "%(self.tag))
        
        if( self.resource.meta['type'] == "local" ):
            os.chdir(self.dir['home'])
            for dkey,dir_i in self.dir.iteritems():
                if ( not os.path.isdir(dir_i) ):
                    logger.info("Making %s "%(dir_i))
                    os.mkdir(dir_i)
            os.chdir(self.dir['home'])
        elif( self.resource.meta['type'] == "ssh" ):
            # Create local directories 
            os.chdir(self.dir['home'])
            dkey = 'launch'
            try:
                dir_i = self.dir[dkey]
                if ( not os.path.isdir(dir_i) ):
                    logger.info("Making %s "%(dir_i))
                    os.mkdir(dir_i)                
            except:
                logger.warning("%s directory not set "%(dkey))
            # Create remote directories 
            ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])#
            for dkey in ['storage','scratch']:
                try:
                    dir_i = self.dir[dkey]
                    # Create  directory
                    bash_command = ' mkdir -p %s  '%(dir_i)
                    bash_command = 'ssh %s \'  %s \' '%(ssh_id,bash_command)
                    logger.debug("MAKEDIR:{}".format(bash_command))
                    os.system(bash_command)                       
                except:
                    logger.warning("%s directory not set "%(dkey))
                    
            
                            
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
        # Set storage and scratch to have calculation directories 
        self.dir['storage'] = '%s/%s/'%(resource_i.dir['storage'] ,self.tag)
        self.dir['scratch'] = '%s/%s/'%(resource_i.dir['scratch'],self.tag)
        self.dir['launch'] = '%s/%s/'%(resource_i.dir['launch'],self.tag)
        self.properties['scratch'] = resource_i.dir['scratch'] 


    def add_refcalc(self,ref_calc):
        '''
        Add reference calculation to current calculation
        this will copy the output of the reference calculation
        to the current calculation's scratch location to be used
        by the current calculation
        '''
        self.references[ref_calc.tag] = ref_calc
        
    def get_cp_str(self,file_type,file_key,file_name,from_dirkey,to_dirkey):
        '''
        Return bash command as string to copy a file
        '''
        cpfile = True
        try:
            fromdir = self.dir[from_dirkey]
        except:
            cpfile = False
            logger.warning(" dir dictionary key %s not set for calculation %s files "%(from_dirkey,self.tag))
        try:
            todir = self.dir[to_dirkey]
        except:
            cpfile = False
            logger.warning(" dir dictionary key %s not set for calculation %s files "%(to_dirkey,self.tag))
        if( cpfile ):
            from_pathfile = os.path.join(fromdir,file_name)
            to_pathfile = os.path.join(todir,file_name)
            return "cp %s %s "%(from_pathfile,to_pathfile)
        else:
            return ''
        
    def cp_file(self,file_type,file_key,file_name,from_dirkey,to_dirkey):
        '''
        Add file to calculation with add_file and copy the file to a directory
        '''
        logger.info(" in cp_file {} {} {} {} {}".format(file_type,file_key,file_name,from_dirkey,to_dirkey))
        
        cpfile = True
        try:
            fromdir = self.dir[from_dirkey]
        except:
            cpfile = False
            logger.info(" dir dictionary key %s not set for calculation %s files "%(from_dirkey,self.tag))
        try:
            todir = self.dir[to_dirkey]
        except:
            cpfile = False
            logger.info(" dir dictionary key %s not set for calculation %s files "%(to_dirkey,self.tag))
        if( cpfile ):
            from_pathfile = os.path.join(fromdir,file_name)
            to_pathfile = os.path.join(todir,file_name)
            logger.info("copying %s to %s "%(from_pathfile,to_pathfile))
            shutil.copyfile(from_pathfile,to_pathfile)
             
        self.add_file(file_type,file_key,file_name)
        #
    def compress_dirfiles(self,file_type,in_dirkey):
        '''
        Compress  files in directory
        '''
        os.chdir(self.dir[in_dirkey])
        self.compress_files(file_type)
        os.chdir(self.dir['home'])


    def push(self,file_type_list=[ 'output', 'data']):
        '''
        Push input files to resource 
        '''
        #
        logger.info(" Resource type %s "%(self.resource.meta['type']))

        if( self.meta['status'] == 'written' ):
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
                logger.info("runnning push function in %s "%(os.getcwd()))
                #
                # Scp compressed files over to resource 
                #
                from_dirkey = 'launch'
                fromdir = self.dir[from_dirkey]
                to_dirkey = 'scratch'
                todir = self.dir[to_dirkey]
                file_key = self.properties['comp_key']
                for file_type in ['input','templates','scripts']:
                    if( len(self.files[file_type]) ):
                        
                        logger.info("Compressing and copying %s files to scratch directory "%(file_type))
                        self.compress_files(file_type)
                        file_name = self.files[file_type][file_key]
                        # self.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
                        from_pathfile = os.path.join(fromdir,file_name)
                        to_pathfile = os.path.join(todir,file_name)
                        bash_command = "scp %s %s:%s"%(from_pathfile,ssh_id,todir)
                        logger.info("COPY:{}".format(bash_command))
                        os.system(bash_command)
                        # Uncompress file 
                        bash_command = "%s %s "%(self.properties['uncompress'],file_name)
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,todir,bash_command)
                        os.system(bash_command)
                    else:
                        logger.info("No files of type %s set"%(file_type))
            else:
                logger.info(" Resource type %s does not need to push files created locally "%(self.resource.meta['type']))
            #
            # Copy reference simulation output over to scratch directory
            #
            file_key = self.properties['comp_key']
            if( self.resource.meta['type'] == "ssh"  ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
            for ref_key,ref_calc in self.references.iteritems():
                logger.info("Copying output of reference calculations %s"%(ref_key))
                # file_type = 'output'
                for file_type in file_type_list:
                    
                    try:
                        file_name = ref_calc.files[file_type][file_key]
                    except:
                        logger.info("Fine type %s has no compressed files "%(file_type))
                        file_name = ''
                    if( len(file_name) > 0 ):
                        # Copy compressed reference output to scratch directory 
                        if( ref_calc.resource.meta['type'] == "local" and self.resource.meta['type'] == "local"):
                            bash_command = "cp %s/%s %s"%(ref_calc.dir['storage'],file_name,self.dir['scratch'])
                            os.system(bash_command)
                        elif( ref_calc.resource.meta['type'] == "ssh" and self.resource.meta['type'] == "ssh" ):
                            if( ref_calc.resource.ssh['address'] == self.resource.ssh['address'] ):
                                # Copy on external resource 
                                bash_command = "cp %s/%s %s"%(ref_calc.dir['storage'],file_name,self.dir['scratch'])
                                bash_command = 'ssh %s \'  %s \' '%(ssh_id,bash_command)
                                os.system(bash_command)                
                            else:
                                # Copy between resources 
                                ref_ssh_id = "%s@%s"%(ref_calc.resource.ssh['username'],ref_calc.resource.ssh['address'])
                                bash_command = "scp %s:%s/%s %s:%s"%(ref_ssh_id,ref_calc.dir['storage'],file_name,ssh_id,self.dir['scratch'])
                                os.system(bash_command)

                        elif( ref_calc.resource.meta['type'] == "ssh" and self.resource.meta['type'] == "local" ):
                                # Copy between resources 
                                ref_ssh_id = "%s@%s"%(ref_calc.resource.ssh['username'],ref_calc.resource.ssh['address'])
                                bash_command = "scp %s:%s/%s %s"%(ref_ssh_id,ref_calc.dir['storage'],file_name,self.dir['scratch'])
                                os.system(bash_command)                                            
                        else:
                            logger.info(" Copy from  type %s to type %s not set  "%(ref_calc.resource.meta['type'], self.resource.meta['type']))
                        # Uncompress reference output 
                        bash_command = "%s %s "%(self.properties['uncompress'],file_name)
                        if( self.resource.meta['type'] == "ssh" ):
                            bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.dir['scratch'],bash_command)
                        # Run bash command
                        os.system(bash_command)

                
                
            '''
            # Copy over reference output
            #for ref_tag in self.files['reference']:
            for ref_sim in self.references:
                #ref_sim = load_json(ref_tag)
                ref_compressed_output = "%s_output.%s"%(ref_sim.tag,ref_sim.compress_sufix)
                if( ref_sim.resource.meta['type'] == "local" and self.resource.meta['type'] == "local"):
                    bash_command = "cp %s/%s %s"%(ref_sim.meta['storage_dir'],ref_compressed_output,self.meta['scratch_dir'])
                    os.system(bash_command)                
                if( ref_sim.resource.meta['type'] == "ssh" and self.resource.meta['type'] == "ssh"):
                    if( ref_sim.resource.meta['address'] == self.resource.meta['address'] ):
                        bash_command = "cp %s/%s %s"%(ref_sim.meta['storage_dir'],ref_compressed_output,self.meta['scratch_dir'])
                        bash_command = 'ssh %s \'  %s \' '%(ssh_id,bash_command)
                        os.system(bash_command)                
                    else:
                        logger.info("Can't transfer simulations between two extrernal resources :( ")
                else:
                    logger.info(" Resource type %s unknown "%(ref_sim.resource.meta['type']))
                #
                # Record compressedinputs
                self.files['compressedinputs'].append(ref_compressed_output)            
                         
            os.chdir(self.meta['launch_dir'])
            #
            # Compress input files 
            self.compress_input()
            #
            # Copy over input files 
            if( self.resource.meta['type'] == "local" ):
                bash_command = "cp %s %s"%(self.compressed_input,self.meta['scratch_dir'])
                os.system(bash_command)                
            elif( self.resource.meta['type'] == "ssh" ):
                bash_command = "scp %s %s:%s"%(self.compressed_input,ssh_id,self.meta['scratch_dir'])
                os.system(bash_command)
            else:
                logger.info(" Resource type %s unknown "%(self.resource.meta['type']))
            #
            # Record compressedinputs
            self.files['compressedinputs'].append(self.compressed_input)


            if( self.resource.meta['type'] == "local" ):
                # Change to scratch directory 
                os.chdir(self.meta['scratch_dir'])
                logger.debug("scratch_dir %s "%(self.meta['scratch_dir']))
                # Uncompress input files

                
                
                for file_i in self.files['compressedinputs']:
                    bash_command = "%s %s "%(self.uncompress_method,file_i)
                    # Run bash command
                    os.system(bash_command)

            elif( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.meta['username'],self.resource.meta['address'])                              
                for file_i in self.files['compressedinputs']:
                    bash_command = "%s %s "%(self.uncompress_method,file_i)
                    bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.meta['scratch_dir'],bash_command)
                    # Run bash command
                    os.system(bash_command)
            '''
        
        
        
            
        
        
