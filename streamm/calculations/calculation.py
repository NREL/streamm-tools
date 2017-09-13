# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

"""
Class data structures for calculation data
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"


import copy, sys, os, shutil, math
import time, datetime
import json
import numpy as np
from string import replace
from datetime import datetime

try:
    # Import pymatgen Class 
    import pymatgen.core.periodic_table as pymatgen_pt
except:
    raise ImportError("pymatgen import error for periodic_table object")


# Import streamm dependencies 
import streamm.util.units as units
from streamm.buildingblocks.container import Container as BBCont
from streamm.forcefields.container import Container as ParamCont


import logging
logger = logging.getLogger(__name__)



class Calculation(object):
    '''
    Data structure for a calculation where input files are read in
    and output files and data files are produced.

    These input, output and data files are tracked within a dictionary
    along with dictionaries of meta data and units.
        
        
    Args:
        tag (str): String identifier for object
    
    Kwargs:
        units_conf (dict): Dictionary of units for each attribute type
        

    .. attribute:: data (dict)
    
        Calculation data
        
    .. attribute:: units (dict)
    
        Track the units used in calculation
        
    .. attribute:: meta (dict)
        
        Track Calculation data
        
    .. attribute:: files (dict)
        
        Files needed and produced by calculation
    
    .. attribute:: files['input']  (dict)
    
        Input files for calculation
        
    .. attribute:: files['templates']  (dict)
        
        Templates files for calculation
        
    .. attribute:: files['scripts']  (dict)
        
        Scripts files for calculation
        
    .. attribute:: files['output'] (dict)
        
        Output files
        
    .. attribute:: files['data']   (dict)
    
        Data produced by calculation
    
    .. attribute:: str['templates'] (dict)
    
        Template files
        
    '''
    def __init__(self, tag,unit_conf=units.unit_conf ):
        
        self.tag = str(tag)
        # Store the units of each attribute type 
        self.unit_conf = unit_conf  
        
        self.prefix = 'calc'
        self.data = dict()

        self.strucC = BBCont(unit_conf=unit_conf )
        self.paramC = ParamCont()
        
        dt = datetime.fromtimestamp(time.time())
        self.meta = dict()
        self.meta['date'] = dt.isoformat()
        self.meta['status'] = 'written'
        self.meta['software'] = 'streamm_calc'
        
        self.units = dict()
        self.units['distance'] = 'angstroms'
        self.units['angle'] = 'radians'
        
        self.files = dict()
        self.files['input'] = dict()     
        self.files['templates'] = dict()
        self.files['scripts'] = dict()
        self.files['output'] = dict()
        self.files['data'] = dict()
        
        
        self.str = dict()

        self.properties = dict()        
        # Add compression properties
        self.properties['compress'] =  "tar -czf "
        self.properties['uncompress'] =  "tar -xzf "
        self.properties['compress_sufix'] =  "tgz"
        self.properties['comp_key'] = 'compressed'
        
        self.references = dict()        
        
        self.files['data']['json'] = "%s_%s.json"%(self.prefix,self.tag)
        
    def __del__(self):
        """
        Delete Calculation object
        """
        del self.prefix
        del self.tag
        del self.strucC
        del self.paramC
        del self.data
        del self.meta
        del self.units        
        del self.files        
        del self.str        
        del self.properties        

    def add_file(self,file_type,file_key,file_name):
        '''
        Add input file to simulation 
        '''
        # Check that file is not already in output dic
        for k_j,f_j in self.files[file_type].iteritems():
            if( f_j == file_name ):
                return
            if( k_j == file_key):
                logger.warning("File %s with key %s in type %s will be replaced with %s "%(f_j,k_j,file_type,file_name))
        # Not in dic add 
        self.files[file_type][file_key] = file_name
        logger.info("add_file file_type %s file_key %s file_name %s "%(file_type,file_key,file_name))
        #
        return    
        #
    def get_compress_str(self,file_type):
        '''
        get bash string to compress files
        '''
        
        def file_test(file_i):
            '''
            Read file into string 
            '''
            try:
                with open(file_i) as F:            
                    F.close()
                    return True  

            except IOError:
                logger.warning(" File not found %s "%(file_i))

            return False

        compress_files = ''
        for fkey,file_i in self.files[file_type].iteritems():
            # if( fkey != self.properties['comp_key'] and file_test(file_i) ):
            if( fkey != self.properties['comp_key']  ):
                logger.info("Adding %s "%(file_i))
                compress_files += ' %s'%(file_i)
        compressed_file = "%s_%s.%s"%(self.tag,file_type,self.properties['compress_sufix'] )
        #
        logger.info(" compressed_file:{} ".format(compressed_file))
        #
        # Add compressed file to file dictionary 
        self.add_file(file_type,self.properties['comp_key'],compressed_file)
        # Return bash command as string 
        return "%s %s %s "%(self.properties['compress'],compressed_file,compress_files)
        
    def compress_files(self,file_type):
        '''
        Compress  files
        '''
        logger.info("Compressing %s files in dir %s "%(file_type,os.getcwd()))
        bash_command = self.get_compress_str(file_type)
        logger.debug(" compressing with %s "%(bash_command))
        # Run bash command
        os.system(bash_command)


    def load_str(self,file_type,file_key):
        '''
        Read in a file as string and store in dictionary
        
        Args:
           file_type (str): dictionary key for files dictionary
           file_key (str): dictionary key for files[file_type] dictionary
        
        '''
        
        file_name = self.files[file_type][file_key]
        try:
            with open(file_name) as F:
                self.str[file_key] = ''
                lines = F.readlines()
                F.close()
                for line in lines:
                    self.str[file_key] += line

        except IOError:
            logger.info(" File not found %s "%(file_name))

    def replace_prop(self,str_key):
        '''
        Replace properties in string with values from properties dictionary 
        '''
        new_str = copy.deepcopy(self.str[str_key])
        # Replace some special values 
        new_str = replace(new_str,"<tag>",self.tag)
        
        for propkey,prop_i in self.properties.iteritems():
                prop_replace = "<%s>"%(propkey)
                new_str = replace(new_str,prop_replace,str(prop_i))

        return new_str

    def replacewrite_prop(self,str_key,file_type,file_key,file_name):
        '''
        Replace properties in string and write to file
        '''
        new_str = self.replace_prop(str_key)
        
        
        F = open(file_name,"w")
        F.write(new_str)
        F.close()

        self.add_file(file_type,file_key,file_name)

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
        
        # Update file in case tag has changed 
        self.files['data']['json'] = "%s_%s.json"%(self.prefix,self.tag)
        
        f = open(self.files['data']['json'], 'w')
        json.dump(json_data,f, indent=2)
        f.close()

    def load_json(self):
        '''
        Load json file for reference 
        '''
        #
        # Update file in case tag has changed
        #
        self.files['data']['json'] = "%s_%s.json"%(self.prefix,self.tag)
        #
        logger.info("Reading %s"%(self.files['data']['json']))
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

        except IOError:
            logger.warning(" File not found %s in %s "%(self.files['data']['json'],os.getcwd()))

    def set_strucC(self,strucC_i):
        '''
        Add a structure.Container to the simulation 
        '''
        self.strucC = copy.deepcopy(strucC_i)

    def add_strucC(self,strucC_i):
        '''
        Add a structure.Container to the simulation 
        '''
        self.strucC += strucC_i
        
    def add_paramC(self,paramC_i):
        '''
        Add force field parameters based on fftypes in structure 
        '''
        self.paramC += paramC_i

    def run(self):
        '''
        Run calculation script
        '''
        logger.info("Calculation with status %s "%(self.meta['status']))
        if( self.meta['status'] == 'written' ):
            logger.info("Resource type %s "%(self.resource.meta['type'] ))
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
            for fkey,file_i in self.files['scripts'].iteritems():
                if( fkey != self.properties['comp_key'] ):
                    
                    bash_command = "chmod 755   %s "%(file_i)
                    if( self.resource.meta['type'] == "ssh" ):
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.dir['scratch'],bash_command)
                    os.system(bash_command)
                        
                    bash_command = "%s%s"%(self.properties['exe_command'],file_i)
                    if( self.resource.meta['type'] == "ssh" ):
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.dir['scratch'],bash_command)
                    logger.info("Executing run command %s "%(bash_command) )
                    os.system(bash_command)

            self.meta['status'] == 'submitted'

    def check(self,output_key='log'):
        '''
        Check if calculation has finished        
        '''
        if( self.meta['status'] != 'stored' ):
        
            # Find output_key file 
            try:
                output_file = self.files['output'][output_key]
            except KeyError:
                logger.info("Calculation %s no output_file file  with key %s found"%(self.tag,output_key))
                return
            
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
                bash_command = "ssh  %s  grep \'%s\' %s%s >> ./%s "%(ssh_id,self.properties['finish_str'],self.dir['scratch'],output_file,output_file)
                os.system(bash_command)        
            elif( self.resource.meta['type'] == "local" ):
                # Change to scratch directory
                logger.debug("Running check in local directory %s "%(os.getcwd()))
            else:
                logger.info(" Resource type %s unknown "%(self.resource.meta['type']))

            # Check output file for key string  
            try:
                with open(output_file) as F:
                    output_lines = F.readlines()
                    F.close()
                    if( len(output_lines) > 0 ):
                        # If output file exists assume it's running 
                        self.meta['status'] = 'running'
                        key_cnts = 0 
                        for line in output_lines:
                            if( self.properties['finish_str'] in line ):
                                self.meta['status'] = 'finished'
                                #return 

            except IOError:
                logger.info("If ouput file missing leave status and return")
                return
        else:
            logger.info("Calculation %s has status %s and will not be checked "%(self.tag,self.meta['status']))
            


    def proc_log(self,log_file):
        """
        Read in new files created by script recorded in the log file.

        Args:
            log_file (str) script log file

        """
        f = open(log_file,'r')
        log_lines = f.readlines()
        f.close()
        for line in log_lines:
            # logger.info("line:",line
            llow = line.lower()
            if( 'file:' in  llow):
                col = line.split()
                file_type = col[8]
                file_key = col[9]
                file_name = col[10]
                logger.info("Adding file type %s with key %s and name %s "%(file_type,file_key,file_name))
                self.add_file(file_type,file_key,file_name)
        #
        #
    def analysis(self,output_key='log'):
        """
        Read in results from streamm script
        

        Args:
            output_key (str): dictionary key for files['output'][output_key]
        
        """
        # Find output_key file 
        try:
            output_file = self.files['output'][output_key]
        except KeyError:
            output_file = ''
            logger.info("Calculation %s No output_file file  with key %s found"%(self.tag,output_key))
        if( len(output_file) > 0 ):
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])                            
                bash_command = "scp  %s:%s%s  ./ "%(ssh_id,self.dir['scratch'],output_file)
                os.system(bash_command)                
            self.proc_log(output_file)            
            
        
        
    def store(self,file_type_list=['input','scripts','output','data']):
        '''
        Copy input, output and data of simulation to storage. 
        
        This assumes storage is accesable by cp from ssh resource
        
        '''
        if( self.meta['status'] == 'finished' ):
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
            logger.info("runnning store function in %s "%(os.getcwd()))
            from_dirkey = 'scratch'
            fromdir = self.dir[from_dirkey]
            to_dirkey = 'storage'
            todir = self.dir[to_dirkey]
            file_key = self.properties['comp_key']
            # for file_type in ['input','scripts','output','data']:
            for file_type in file_type_list:
                if( len(self.files[file_type]) ):
                    logger.info("Storing %s files "%(file_type))
                    bash_command = self.get_compress_str(file_type)
                    if( self.resource.meta['type'] == "ssh" ):
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.dir['scratch'],bash_command)
                    logger.info(bash_command)
                    os.system(bash_command)
                    file_name = self.files[file_type][file_key]
                    bash_command = self.get_cp_str(file_type,file_key,file_name,from_dirkey,to_dirkey)
                    if( self.resource.meta['type'] == "ssh" ):
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.dir['scratch'],bash_command)
                    logger.info(bash_command)
                    os.system(bash_command)
                    self.add_file(file_type,file_key,file_name)
                else:
                    logger.info("No files of type %s present"%(file_type))
            self.meta['status'] = 'stored'
            

    def pull(self,file_type_list=['output','data']):
        '''
        Copy output and data of simulation from storage
        
        '''
        if( self.meta['status'] == 'stored' ):
            if( self.resource.meta['type'] == "local" ):
                logger.info("runnning pull function in %s "%(os.getcwd()))
                from_dirkey = 'storage'
                to_dirkey = 'scratch'
                file_key = self.properties['comp_key']
                for file_type in file_type_list:
                    file_name = self.files[file_type][file_key]
                    self.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
                    os.chdir( self.dir[to_dirkey] )
                    bash_command = "%s %s "%(self.properties['uncompress'],file_name)
                    os.system(bash_command)
                
            elif( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
                from_dirkey = 'scratch'
                fromdir = self.dir[from_dirkey]
                to_dirkey = 'launch'
                todir = self.dir[to_dirkey]
                comp_key = self.properties['comp_key']
                #for file_type,files in self.files:
                for file_type in file_type_list:
                    #if( len(self.files[file_type]) ):
                    run_cp  = False 
                    try:
                        file_name = self.files[file_type][comp_key]
                        run_cp  = True 
                        logger.info(file_type,self.files[file_type][comp_key])
                    except:
                        logging.warning("No commpressed files for file type %s with key %s "%(file_type,comp_key))
                    if( run_cp ):
                        from_pathfile = os.path.join(fromdir,file_name)
                        to_pathfile = os.path.join(todir,file_name)
                        bash_command = "scp %s:%s%s %s "%(ssh_id,fromdir,file_name,todir)
                        logger.debug("bash_command:{} ".format(bash_command))
                        os.system(bash_command)
                        # Uncompress file
                        os.chdir(todir)
                        bash_command = "%s %s "%(self.properties['uncompress'],file_name)
                        os.system(bash_command)
            else:
                logger.info(" Resource type %s unknown "%(self.resource.meta['type']))
            

    def set_ffparam(self):
        '''
        Create new parameter container to hold each unique type.
        
        This is necessary for parameter outputs to no have
        redundent values.
        
           
        '''
        use_last = True 
        # Open log file 
        log_file = "param.log"
        param_out = open(log_file,"w")
        #
        # Create a container of possible parameters
        #
        self.paramC_o = copy.deepcopy( self.paramC)
        # Reset parameters container
        self.paramC = ParamCont()
        #
        # Set general parameters 
        #
        self.paramC.nbfunc = self.paramC_o.nbfunc
        self.paramC.combmixrule = self.paramC_o.combmixrule
        self.paramC.genpairs = self.paramC_o.genpairs
        self.paramC.fudgeJ = self.paramC_o.fudgeLJ
        self.paramC.fudgeQQ = self.paramC_o.fudgeQQ
        #
        # Examine atom types
        #
        for pkey_o, particle_o  in self.strucC.particles.iteritems():    
            new_type = True
            fftype_i = particle_o.ffkey
            #mass_i = particle_o.mass
            lmptype_p = 0

            logger.debug(" Checking type {}".format(fftype_i))

            for ptkey_i, ptype_i  in self.paramC.particletypes.iteritems():
                lmptype_p = ptkey_i + 1
                if( fftype_i == ptype_i.fftype1 ):
                    new_type = False
                    # particle_o.type = str(lj_p)
                    particle_o.ff = copy.deepcopy(ptype_i) 
                    particle_o.ff.lammps_index =  lmptype_p
                    logger.debug(' maches previous {} {} {} '.format(pkey_o, ptkey_i,ptype_i.fftype1))

            if( new_type ):
                lmptype_p += 1
                logger.debug(" New type {} {}".format(fftype_i,lmptype_p))
                # Find type in ljtypC_all
                cnt_check = 0
                type_found = False 
                # particle_o.type = str(lj_p+1)
                for lj_all, ljObj_all  in self.paramC_o.particletypes.iteritems():
                    all_i = ljObj_all.fftype1
                    if( fftype_i == all_i ):
                        type_found = True

                    if( type_found ):
                        cnt_check += 1
                        ljObj_all.lammps_index = lmptype_p
                        # set_obj.setpid( particle_o.properties["number"] )
                        #ljObj_all.setmass( mass_i )
                        self.paramC.add_particletype(ljObj_all,deepcopy = True)

                        particle_o.ff = copy.deepcopy(ljObj_all) 
                        particle_o.ff.lammps_index = lmptype_p
                        
                        type_found = False 

                if( cnt_check < 1 ):
                    self.paramC = self.paramC_o
                    raise TypeError(" No particle parameters were found for atom type %s "%fftype_i)
                     
                elif( cnt_check > 1 ):
                    logger.warning(" Multiple particle parameters (%d) were found for atom type %s "%(cnt_check,fftype_i))
                    if( not use_last  ):
                        self.paramC = self.paramC_o
                        raise TypeError("Last parameter will not be used")
                    
        #
        # Examine  bonds types
        #
        for bkey_o, bondObj_o  in self.strucC.bonds.iteritems():
            new_type = True
            lmptype_p = 0
            pid_i = bondObj_o.pkey1
            pid_j = bondObj_o.pkey2
            fftype_i =  self.strucC.particles[ pid_i ].ffkey
            fftype_j =  self.strucC.particles[ pid_j ].ffkey
            #r_i = np.array( self.strucC.particles[ bondObj_o.pkey1 ].position  )
            #r_j = np.array( self.strucC.particles[ bondObj_o.pkey2 ].position  )
            #bond_len = np.linalg.norm(delta_r_c(r_i,r_j,struc_o.getLatVec() ))

            for btyp_p, btypObj_p  in self.paramC.bondtypes.iteritems():
                lmptype_p = btypObj_p.lammps_index
                p_i = btypObj_p.fftype1 
                p_j = btypObj_p.fftype2
                match = False 
                if( fftype_i == p_i  and  fftype_j == p_j ):
                    match = True
                elif( fftype_i == p_j  and  fftype_j == p_i ):
                    match = True
                if( match ):
                    new_type = False
                    bondObj_o.lammps_index = btypObj_p.lammps_index
                    bondObj_o.gromacs_index = btypObj_p.gromacs_index
                    # log_line=" Setting bond atoms %s - %s numbers %d - %d wiht bond length %f to type %d with r_o %f  delta %f \n"%(fftype_i,fftype_j,pid_i,pid_j,bond_len,btyp_p,btypObj_p.get_r0(),bond_len-btypObj_p.get_r0() )
                    #log_line=" Setting to lmptyp %d bond atoms %s - %s numbers %d - %d "%(lmptype_p,fftype_i,fftype_j,pid_i,pid_j)
                    #param_out.write(log_line+'\n')
                    break 

            if( new_type ):
                # Find type in btypC_all
                lmptype_p += 1
                cnt_check = 0
                type_found = False 
                copy_type = False 
                
                for b_all, btypObj_all  in self.paramC_o.bondtypes.iteritems():
                    all_i = btypObj_all.fftype1 
                    all_j = btypObj_all.fftype2
                    if( fftype_i == all_i  and  fftype_j == all_j ):
                        copy_type = True 
                    elif( fftype_j == all_i  and  fftype_i == all_j ):
                        copy_type = True

                    if( copy_type ):
                        cnt_check += 1
                        bondObj_temp = copy.deepcopy(btypObj_all)
                        type_found = True
                        copy_type = False 

                if( cnt_check < 1 ):
                    self.paramC = self.paramC_o
                    raise TypeError(" No Bond parameters were found for bond type %s-%s "%(fftype_i,fftype_j))
                elif( cnt_check > 1 ):
                    logger.warning(" %d  Bond parameters were found for bond type %s-%s "%(cnt_check,fftype_i,fftype_j))
                    for btyp_p, btypObj_p  in self.paramC.bondtypes.iteritems():
                        logger.info(btyp_p ,btypObj_p.fftype1 ,btypObj_p.fftype2      )              
                    if( not use_last  ):
                        self.paramC = self.paramC_o
                        raise TypeError("Last parameter will not be used")
                        
                if( type_found ):
                        
                        bondObj_o.gromacs_index = bondObj_temp.gromacs_index
                        # Set lammps type for bond and bond type
                        # These have to match!!! 
                        bondObj_o.lammps_index = lmptype_p          # Set LAMMPS index for bond
                        bondObj_temp.lammps_index = lmptype_p
                        self.paramC.add_bondtype(bondObj_temp)
                        logger.debug(" %d  Bond parameters were found for bond type %s-%s "%(cnt_check,fftype_i,fftype_j))

                        # log_line=" Setting bond atoms %s - %s numbers %d - %d wiht bond length %f to type %d with r_o %f  delta %f \n"%(fftype_i,fftype_j,pid_i,pid_j,bond_len,btyp_p,btypObj_all.get_r0(),bond_len-btypObj_all.get_r0() )
                        log_line="Adding new lmptyp %d for bond atoms %s - %s numbers %d - %d "%(lmptype_p,fftype_i,fftype_j,pid_i,pid_j)
                        param_out.write(log_line+'\n')


        #
        # Examine  angles types
        #
        # for akey_i, angle_i in self.strucC.angles.iteritems():
        for a_o,angleObj_o in self.strucC.angles.iteritems():
            new_type = True
            lmptype_p = 0
            pid_k = angleObj_o.pkey1
            pid_i = angleObj_o.pkey2
            pid_j = angleObj_o.pkey3
            fftype_k =  self.strucC.particles[ angleObj_o.pkey1 ].ffkey
            fftype_i =  self.strucC.particles[ angleObj_o.pkey2 ].ffkey
            fftype_j =  self.strucC.particles[ angleObj_o.pkey3 ].ffkey
            #r_k = np.array( self.strucC.particles[ pid_k ].position  )
            #r_i = np.array( self.strucC.particles[ pid_i ].position  )
            #r_j = np.array( self.strucC.particles[ pid_j ].position  )
            #r_ik = delta_r_c(r_i,r_k,struc_o.getLatVec() )
            #r_ij = delta_r_c(r_i,r_j,struc_o.getLatVec() )
            #angle_kij = getAngle(r_ik,r_ij)
            for atyp_p, atypObj_p  in self.paramC.angletypes.iteritems():
                lmptype_p = atyp_p + 1
                p_k = atypObj_p.fftype1 
                p_i = atypObj_p.fftype2 
                p_j = atypObj_p.fftype3
                theta0_kij = atypObj_p.theta0
                # delta_theta = angle_kij - theta0_kij
                match = False 
                if( fftype_k == p_k  and  fftype_i == p_i  and  fftype_j == p_j ):
                    match = True
                elif( fftype_j == p_k  and  fftype_i == p_i  and  fftype_k == p_j ):
                    match = True
                if( match ):
                    new_type = False
                    angleObj_o.lammps_index = lmptype_p
                    angleObj_o.gromacs_index = atypObj_p.gromacs_index

                    #log_line=" Setting angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f to type %d with theta_o %f  delta %f \n"%(fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij,atyp_p,theta0_kij,delta_theta )
                    log_line=" Setting angle atoms %s - %s - %s numbers %d - %d - %d "%(fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j)
                    param_out.write(log_line+'\n')
                    break
                
            if( new_type ):
                # Find type in btypC_all
                lmptype_p += 1
                cnt_check = 0
                type_found = False 
                copy_type = False 
                for a_all, atypObj_all  in self.paramC_o.angletypes.iteritems():
                    all_k = atypObj_all.fftype1 
                    all_i = atypObj_all.fftype2 
                    all_j = atypObj_all.fftype3
                    theta0_kij =  atypObj_all.theta0
                    # delta_theta = angle_kij-theta0_kij
                    if( fftype_k == all_k  and fftype_i == all_i  and  fftype_j == all_j ):
                        copy_type = True 
                    elif( fftype_j == all_k  and  fftype_i == all_i  and  fftype_k == all_j ):
                        copy_type = True 

                    if( copy_type ):
                        cnt_check += 1
                        atypObj_temp = copy.deepcopy(atypObj_all)
                        type_found = True
                        copy_type = False 

                if( cnt_check < 1 ):
                    self.paramC = self.paramC_o
                    raise TypeError(" No Angles parameters were found for bond type %s-%s-%s "%(fftype_k,fftype_i,fftype_j))
                elif( cnt_check > 1 ):
                    # log_line=" %d Angles parameters were found for angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f  \n"%(cnt_check,fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij )
                    logger.warning(" %d Angles parameters were found for angle atoms %s - %s - %s numbers %d - %d - %d   \n"%(cnt_check,fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j ))
                    #param_out.write(log_line)
                    logger.info(log_line)
                    # atypC_p.findtype(fftype_k,fftype_i,fftype_j)

                    if( not use_last  ):
                        self.paramC = self.paramC_o
                        logger.warning("Last parameter will not be used")
                if( type_found ):
                        
                        angleObj_o.gromacs_index = atypObj_temp.gromacs_index
                        angleObj_o.lammps_index = lmptype_p
                        atypObj_temp.lammps_index = lmptype_p
                        self.paramC.add_angletype(atypObj_temp)
                        log_line=" Adding new lmptyp %d  angle atoms %s - %s - %s numbers %d - %d - %d "%(lmptype_p,fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j)
                        logger.debug(" %d Angles parameters were found for bond type %s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j))
                        param_out.write(log_line+'\n')


        #
        # Examine  dihedrals types
        #
        imp_cnt = 0 
        for d_o,dihObj_o in self.strucC.dihedrals.iteritems():
            new_type = True
            lmptype_p = 0
            pid_k = dihObj_o.pkey1
            pid_i = dihObj_o.pkey2 
            pid_j = dihObj_o.pkey3
            pid_l = dihObj_o.pkey4
            fftype_k =  self.strucC.particles[pid_k].ffkey
            fftype_i =  self.strucC.particles[pid_i].ffkey
            fftype_j =  self.strucC.particles[pid_j].ffkey
            fftype_l =  self.strucC.particles[pid_l].ffkey

            logger.debug(" checking ",fftype_k, fftype_i,  fftype_j , fftype_l)
            # Check to see if dihedral type is already in parameter set for the structure container
            for dtyp_p, dtypObj_p  in self.paramC.dihtypes.iteritems():
                lmptype_p = dtyp_p + 1
                p_k = dtypObj_p.fftype1 
                p_i = dtypObj_p.fftype2 
                p_j = dtypObj_p.fftype3
                p_l = dtypObj_p.fftype4
                match = False 
                if( fftype_k == p_k  and  fftype_i == p_i  and  fftype_j == p_j  and  fftype_l == p_l ):
                    match = True
                if( fftype_l == p_k  and  fftype_j == p_i  and  fftype_i == p_j  and  fftype_k == p_l ):
                    match = True
                if( match ):
                    new_type = False
                    dihObj_o.lammps_index = lmptype_p
                    dihObj_o.gromacs_index = dtypObj_p.gromacs_index

                    logger.debug(" dihObj_o.lammps_index  ",dihObj_o.lammps_index)
                    logger.debug(" dihObj_o.gromacs_index  ",dihObj_o.gromacs_index)
                    logger.debug("  previous type ",dtyp_p,p_k,p_i,p_j,p_l,dihObj_o.gromacs_index)
                    break
                    

            # If it is not in the parameter set for the struture container
            #  find it in the parameters from the reference parameter file 
            if( new_type ):
                # Find type in btypC_all
                lmptype_p += 1
                cnt_check = 0
                type_found = False 
                # Set type to new type = last type+1
                dihObj_o.lammps_index = lmptype_p

                logger.debug("  new type checking against %d read in parameters "%len(self.paramC_o.dihtypes))

                copy_type = False 
                for d_all, dtypObj_all  in self.paramC_o.dihtypes.iteritems():
                    all_k = dtypObj_all.fftype1 
                    all_i = dtypObj_all.fftype2 
                    all_j = dtypObj_all.fftype3
                    all_l = dtypObj_all.fftype4

                    if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                        copy_type = True    
                    elif ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l   ):
                        copy_type = True

                    if( copy_type ):
                        cnt_check += 1
                        dtypObj_temp = copy.deepcopy(dtypObj_all)
                        #dtypObj_temp.set_g_indx(dtypObj_all.gromacs_index)
                        type_found = True
                        copy_type = False 
                        logger.debug(" %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                        logger.debug("     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l)

                if( not type_found ):
                    logger.debug(" checking  X - FF - FF - FF ")
                    copy_type = False 
                    for d_all, dtypObj_all  in self.paramC_o.dihtypes.iteritems():
                        all_k = dtypObj_all.fftype1 
                        all_i = dtypObj_all.fftype2 
                        all_j = dtypObj_all.fftype3
                        all_l = dtypObj_all.fftype4

                        if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                            copy_type = True
                        if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                            copy_type = True
                        if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l  ):
                            copy_type = True
                        if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                            copy_type = True

                        if( copy_type ):
                            cnt_check += 1
                            dtypObj_temp = copy.deepcopy(dtypObj_all)
                            #dtypObj_temp.set_g_indx(dtypObj_all.gromacs_index)
                            type_found = True 
                            copy_type = False 
                            logger.debug(" %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                            logger.debug("     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l)

                if( not type_found ):
                    logger.debug(" checking  X - FF - FF - X ")
                    copy_type = False 
                    for d_all, dtypObj_all  in self.paramC_o.dihtypes.iteritems():
                        all_k = dtypObj_all.fftype1 
                        all_i = dtypObj_all.fftype2 
                        all_j = dtypObj_all.fftype3
                        all_l = dtypObj_all.fftype4

                        if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                            copy_type = True
                        if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                            copy_type = True

                        if( copy_type ):
                            cnt_check += 1
                            dtypObj_temp = copy.deepcopy(dtypObj_all)
                            #dtypObj_temp.set_g_indx(dtypObj_all.gromacs_index)
                            type_found = True 
                            copy_type = False 
                            logger.debug(" %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                            logger.debug("     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l)

                if( cnt_check < 1 ):
                    self.paramC = self.paramC_o
                    raise TypeError(" No Dih parameters were found for dih type %s-%s-%s-%s "%(fftype_k,fftype_i,fftype_j,fftype_l))
                elif( cnt_check > 1 ):
                    logger.info(" %d Dih parameters were found for dih type %s-%s-%s-%s please check parameter file  "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                    logger.info(dtypObj_temp)
                    #dtypObj_temp_list.findtype(fftype_k,fftype_i,fftype_j,fftype_l)
                    if( not use_last  ):
                        self.paramC = self.paramC_o
                        logger.warning('Last will not be used')

                if( type_found ):
                    logger.debug(" adding new type to dtypC_p  from ",dtypObj_temp.type, dtypObj_temp.fftype1,dtypObj_temp.fftype2,dtypObj_temp.fftype3,dtypObj_temp.fftype4)

                    # Set FF types to read in bond to remove X's 
                    dtypObj_temp.fftype1 = fftype_k
                    dtypObj_temp.fftype2 = fftype_i
                    dtypObj_temp.fftype3 = fftype_j
                    dtypObj_temp.fftype4 = fftype_l
                    
                    '''
                    Remove for now 
                    if( norm_dihparam ):
                        # normalize by number of nieghbors
                        dihen_norm = 1.0
                        if( debug):
                            logger.info(" Normalizing dihedral potential "
                            logger.info(" finding types for ",pid_i,pid_j
                        NNAB_i = calc_nnab(pid_i,cov_nbindx) - 1
                        NNAB_j = calc_nnab(pid_j,cov_nbindx) - 1

                        dihen_norm = float( NNAB_i + NNAB_j)/2.0

                        if(debug): logger.info(" dihen_norm ",dihen_norm

                        dtypObj_temp.normforceconstants(dihen_norm)
                    '''

                    dihObj_o.gromacs_index = dtypObj_temp.gromacs_index
                    dtypObj_temp.lammps_index = lmptype_p
                    self.paramC.add_dihtype(dtypObj_temp)                
                    logger.debug(" %d Dih parameters were found for dih type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                    logger.debug(" len(dtypC_p) ",len(self.paramC.dihtypes) )
                    logger.debug(" dtypObj_temp.lammps_index  ",dtypObj_temp.lammps_index)
                    logger.debug(" dtypObj_temp.gromacs_index  ",dtypObj_temp.gromacs_index)
        #
        # Examine improper dihedrals types
        #
        imp_cnt = 0 
        for imp_o,impObj_o in self.strucC.impropers.iteritems():
            new_type = True
            lmptype_p = 0
            pid_k = impObj_o.pkey1
            pid_i = impObj_o.pkey2 
            pid_j = impObj_o.pkey3
            pid_l = impObj_o.pkey4 
            fftype_k =  self.strucC.particles[pid_k].ffkey
            fftype_i =  self.strucC.particles[pid_i].ffkey
            fftype_j =  self.strucC.particles[pid_j].ffkey
            fftype_l =  self.strucC.particles[pid_l].ffkey

            logger.debug(" checking ",fftype_k, fftype_i,  fftype_j , fftype_l)
            # Check to see if impedral type is already in parameter set for the structure container
            for imptyp_p, imptypObj_p  in self.paramC.imptypes.iteritems():
                lmptype_p = imptyp_p + 1
                p_k = imptypObj_p.fftype1 
                p_i = imptypObj_p.fftype2 
                p_j = imptypObj_p.fftype3
                p_l = imptypObj_p.fftype4
                match = False 
                if( fftype_k == p_k  and  fftype_i == p_i  and  fftype_j == p_j  and  fftype_l == p_l ):
                    match = True
                if( fftype_l == p_k  and  fftype_j == p_i  and  fftype_i == p_j  and  fftype_k == p_l ):
                    match = True
                if( match ):
                    new_type = False
                    impObj_o.lammps_index = lmptype_p
                    impObj_o.gromacs_index = imptypObj_p.gromacs_index

                    logger.debug(" impObj_o.lammps_index  ",impObj_o.lammps_index)
                    logger.debug(" impObj_o.gromacs_index  ",impObj_o.gromacs_index)
                    logger.debug("  previous type ",imptyp_p,p_k,p_i,p_j,p_l,impObj_o.gromacs_index)
                    break
                    
            # If it is not in the parameter set for the struture container
            #  find it in the parameters from the reference parameter file 
            if( new_type ):
                # Find type in btypC_all
                lmptype_p += 1
                cnt_check = 0
                type_found = False 
                # Set type to new type = last type+1
                impObj_o.lammps_index = lmptype_p

                logger.debug("  new type checking against %d read in parameters "%len(imptypC_all))

                copy_type = False 
                for d_all, imptypObj_all  in self.paramC_o.imptypes.iteritems():
                    all_k = imptypObj_all.fftype1 
                    all_i = imptypObj_all.fftype2 
                    all_j = imptypObj_all.fftype3
                    all_l = imptypObj_all.fftype4

                    if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                        copy_type = True    
                    if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l   ):
                        copy_type = True

                    if( copy_type ):
                        cnt_check += 1
                        imptypObj_temp = copy.deepcopy(imptypObj_all)
                        #imptypObj_temp.set_g_indx(imptypObj_all.gromacs_index)
                        type_found = True
                        copy_type = False 
                        if( debug ):
                            logger.info(" %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                            logger.info("     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l)

                if( not type_found ):
                    logger.debug(" checking  X - FF - FF - FF ")
                    copy_type = False 
                    for d_all, imptypObj_all  in self.paramC_o.imptypes.iteritems():
                        all_k = imptypObj_all.fftype1 
                        all_i = imptypObj_all.fftype2 
                        all_j = imptypObj_all.fftype3
                        all_l = imptypObj_all.fftype4

                        if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                            copy_type = True
                        if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                            copy_type = True
                        if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l  ):
                            copy_type = True
                        if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                            copy_type = True

                        if( copy_type ):
                            cnt_check += 1
                            imptypObj_temp = copy.deepcopy(imptypObj_all)
                            #imptypObj_temp.set_g_indx(imptypObj_all.gromacs_index)
                            type_found = True 
                            copy_type = False 
                            logger.debug(" %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                            logger.debug("     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l)

                if( not type_found ):
                    logger.debug(" checking  X - FF - FF - X ")
                    copy_type = False 
                    for d_all, imptypObj_all  in self.paramC_o.imptypes.iteritems():
                        all_k = imptypObj_all.fftype1 
                        all_i = imptypObj_all.fftype2 
                        all_j = imptypObj_all.fftype3
                        all_l = imptypObj_all.fftype4

                        if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                            copy_type = True
                        if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                            copy_type = True

                        if( copy_type ):
                            cnt_check += 1
                            imptypObj_temp = copy.deepcopy(imptypObj_all)
                            #imptypObj_temp.set_g_indx(imptypObj_all.gromacs_index)
                            type_found = True 
                            copy_type = False 
                            logger.debug(" %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                            logger.debug("     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l)

                if( cnt_check < 1 ):
                    self.paramC = self.paramC_o
                    raise TypeError(" No Dih parameters were found for dih type %s-%s-%s-%s "%(fftype_k,fftype_i,fftype_j,fftype_l))
                elif( cnt_check > 1 ):
                    logger.info(" %d Dih parameters were found for dih type %s-%s-%s-%s please check parameter file  "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                    logger.info(imptypObj_temp)
                    #imptypObj_temp_list.findtype(fftype_k,fftype_i,fftype_j,fftype_l)
                    if( not use_last  ):
                        self.paramC = self.paramC_o
                        raise TypeError

                if( type_found ):
                    logger.debug(" adding new type to imptypC_p  from ",imptypObj_temp.type, imptypObj_temp.fftype1,imptypObj_temp.fftype2,imptypObj_temp.fftype3,imptypObj_temp.fftype4)

                    # Set FF types to read in bond to remove X's 
                    imptypObj_temp.fftype1 = fftype_k
                    imptypObj_temp.fftype2 = fftype_i
                    imptypObj_temp.fftype3 = fftype_j
                    imptypObj_temp.fftype4 = fftype_l

                    '''
                    norm_impdihparam = False 
                    if( norm_impdihparam ):
                        # normalize by number of nieghbors
                        dihen_norm = 1.0
                        if( debug):
                            logger.info(" Normalizing dihedral potential "
                            logger.info(" finding types for ",pid_i,pid_j
                        NNAB_i = calc_nnab(pid_i,cov_nbindx) - 1
                        NNAB_j = calc_nnab(pid_j,cov_nbindx) - 1

                        dihen_norm = float( NNAB_i + NNAB_j)/2.0

                        if(debug): logger.info(" dihen_norm ",dihen_norm
                        logger.info(" dihen_norm ",dihen_norm

                        imptypObj_temp.normforceconstants(dihen_norm)
                    '''
                    impObj_o.gromacs_index = imptypObj_temp.gromacs_index
                    imptypObj_temp.lammps_index = lmptype_p
                    self.paramC.add_imptype(imptypObj_temp,deepcopy = True)
                    logger.debug(" %d Dih parameters were found for dih type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l))
                    logger.debug(" len(imptypC_p) ",len(imptypC_p) )
                    logger.debug(" imptypObj_temp.lammps_index  ",imptypObj_temp.lammps_index)
                    logger.debug(" imptypObj_temp.gromacs_index  ",imptypObj_temp.gromacs_index)


        param_out.close()

        return           
        
        