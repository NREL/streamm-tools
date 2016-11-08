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

import periodictable, units, structure, parameters, buildingblock, resource

import logging
logger = logging.getLogger(__name__)



class Calculation:
    '''
    Data structure for a calculation where input files are read in
    and output files and data files are produced.

    These input, output and data files are tracked within a dictionary
    along with dictionaries of meta data and units.
        
    '''

    def __init__(self, tag):
        """
        Constructor for a general Calculation object.
        
        Args:
            tag (str): String identifier for object

        self.data (dict) Calculation data 
        self.units (dict) Track the units used in calculation
        self.meta (dict) Track Calculation data 
        self.files (dict) Files needed and produced by calculation
            self.files['input']  (dict) Input files for calculation
            self.files['templates']  (dict) Templates files for calculation
            self.files['scripts']  (dict) Scripts files for calculation
            self.files['output'] (dict) Output files 
            self.files['data']   (dict) Data produced by calculation
        self.str['templates'] (dict) of template strings
        """
        if isinstance(tag, str):
            self.tag = tag
        else:
            raise TypeError("1st arg (tag) in %s Container initialization should be string"%(__name__))
        
        self.prefix = 'calc'
        self.data = dict()

        self.strucC = buildingblock.Container()
        self.paramC = parameters.Container()
        
        dt = datetime.datetime.fromtimestamp(time.time())
        self.meta = dict()
        self.meta['date'] = dt.isoformat()
        self.meta['status'] = 'written'
        
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
        compress_files = ''
        for fkey,file_i in self.files[file_type].iteritems():
            print " file_i ",file_i
            if( fkey != self.properties['comp_key'] ):
                compress_files += ' %s'%(file_i)
        compressed_file = "%s_%s.%s"%(self.tag,file_type,self.properties['compress_sufix'] )
        #
        print "> compressed_file ",compressed_file
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
            print " File not found %s "%(file_name)

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

        except IOError:
            logger.warning(" File not found %s in %s "%(json_file,os.getcwd()))
            print " File not found %s in %s "%(json_file,os.getcwd())

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
        print "Calculation with status %s "%(self.meta['status'])
        if( self.meta['status'] == 'written' ):
            print "Resource type %s "%(self.resource.meta['type'] )
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
                    print "Executing run command %s "%(bash_command) 
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
                print "Calculation %s no output_file file  with key %s found"%(self.tag,output_key)

            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
                bash_command = "ssh  %s  grep \'%s\' %s%s >> ./%s "%(ssh_id,self.properties['finish_str'],self.dir['scratch'],output_file,output_file)
                os.system(bash_command)        
            elif( self.resource.meta['type'] == "local" ):
                # Change to scratch directory
                print "Running check in local directory %s "%(os.getcwd())
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
                print "If ouput file missing leave status and return"
                return
        else:
            print "Calculation %s has status %s and will not be checked "%(self.tag,self.meta['status'])

    def store(self):
        '''
        Copy input, output and data of simulation to storage 
        # This assumes storage is accesable by cp from ssh resource 
        '''
        if( self.meta['status'] == 'finished' ):
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
            print "runnning store function in %s "%(os.getcwd())
            from_dirkey = 'scratch'
            fromdir = self.dir[from_dirkey]
            to_dirkey = 'storage'
            todir = self.dir[to_dirkey]
            file_key = self.properties['comp_key']
            for file_type in ['input','scripts','output','data']:
                if( len(self.files[file_type]) ):
                    print "Storing %s files "%(file_type)
                    bash_command = self.get_compress_str(file_type)
                    if( self.resource.meta['type'] == "ssh" ):
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.dir['scratch'],bash_command)
                    os.system(bash_command)
                    file_name = self.files[file_type][file_key]
                    bash_command = self.get_cp_str(file_type,file_key,file_name,from_dirkey,to_dirkey)
                    if( self.resource.meta['type'] == "ssh" ):
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,self.dir['scratch'],bash_command)
                    os.system(bash_command)
                    self.add_file(file_type,file_key,file_name)
                else:
                    print "No files of type %s present"%(file_type)
            self.meta['status'] = 'stored'
            

    def pull(self,file_type='output'):
        '''
        Copy output and data of simulation from storage 
        '''
        print "394jr092"
        if( self.meta['status'] == 'stored' ):
            if( self.resource.meta['type'] == "local" ):
                print "runnning pull function in %s "%(os.getcwd())
                from_dirkey = 'storage'
                to_dirkey = 'scratch'
                file_key = self.properties['comp_key']
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
                for file_type in ['output','data']:
                    #if( len(self.files[file_type]) ):
                    run_cp  = False 
                    try:
                        file_name = self.files[file_type][comp_key]
                        run_cp  = True 
                        print file_type,self.files[file_type][comp_key]
                    except:
                        logging.warning("No commpressed files for file type %s with key %s "%(file_type,comp_key))
                        print "No commpressed files for file type %s with key %s "%(file_type,comp_key)
                    if( run_cp ):
                        from_pathfile = os.path.join(fromdir,file_name)
                        to_pathfile = os.path.join(todir,file_name)
                        bash_command = "scp %s:%s%s %s "%(ssh_id,fromdir,file_name,todir)
                        print "bash_command ",bash_command
                        os.system(bash_command)
                        # Uncompress file
                        os.chdir(todir)
                        bash_command = "%s %s "%(self.properties['uncompress'],file_name)
                        os.system(bash_command)
            else:
                print " Resource type %s unknown "%(self.resource.meta['type'])
            

    def set_ffparam(self):
        '''
         Create new parameter container to hold each unique type
           which exsists in the structure container strucC
           this is necessary for parameter outputs to no have
           redundent values
           
        '''
        use_last = True 
        # Open log file 
        log_file = "param.log"
        param_out = open(log_file,"w")
        #
        # Create a container of possible parameters
        #
        paramC_o = copy.deepcopy( self.paramC)
        # Reset parameters container
        self.paramC = parameters.Container()
        #
        # Set general parameters 
        #
        self.paramC.nbfunc = paramC_o.nbfunc
        self.paramC.combmixrule = paramC_o.combmixrule
        self.paramC.genpairs = paramC_o.genpairs
        self.paramC.fudgeJ = paramC_o.fudgeLJ
        self.paramC.fudgeQQ = paramC_o.fudgeQQ
        #
        # Examine atom types
        #
        debug_lj = False
        for pkey_o, particle_o  in self.strucC.particles.iteritems():    
            new_type = True
            fftype_i = particle_o.properties["fftype"]
            #mass_i = particle_o.mass
            lmptype_p = 0

            if( debug_lj ):
                print " Checking type ",fftype_i

            for lj_p, ljObj_p  in self.paramC.ljtypes.iteritems():
                lmptype_p = lj_p + 1
                if( fftype_i == ljObj_p.fftype1 ):
                    new_type = False
                    # particle_o.type = str(lj_p)
                    particle_o.properties["lmpindx"] = lmptype_p
                    if(debug_lj):
                        print ' maches previous ',pkey_o, lj_p,ljObj_p.fftype1

            if( new_type ):
                lmptype_p += 1
                if( debug_lj ):
                    print " New type ",fftype_i,lmptype_p
                # Find type in ljtypC_all
                cnt_check = 0
                type_found = False 
                # particle_o.type = str(lj_p+1)
                particle_o.properties["lmpindx"] = lmptype_p
                for lj_all, ljObj_all  in paramC_o.ljtypes.iteritems():
                    all_i = ljObj_all.fftype1
                    if( fftype_i == all_i ):
                        type_found = True

                    if( type_found ):
                        cnt_check += 1
                        ljObj_all.lmpindx = lmptype_p
                        # set_obj.setpid( particle_o.properties["number"] )
                        #ljObj_all.setmass( mass_i )
                        self.paramC.add_LJtype(ljObj_all,deepcopy = True)
                        type_found = False 


                if( cnt_check < 1 ):
                    raise TypeError(" No LJ parameters were found for atom type %s "%fftype_i)
                elif( cnt_check > 1 ):
                    logger.warning(" Multiple LJ parameters (%d) were found for atom type %s "%(cnt_check,fftype_i))
                    if( not use_last  ):
                        raise TypeError("Last parameter will not be used")
        if( debug_lj ):

            for lj_p, ljObj_p  in self.paramC.ljtypes.iteritems():
                print ljObj_p
            sys.exit("LJ debug 89798")

        #
        # Examine  bonds types
        #
        debug = 0
        for bkey_o, bondObj_o  in self.strucC.bonds.iteritems():
            new_type = True
            lmptype_p = 0
            pid_i = bondObj_o.pkey1
            pid_j = bondObj_o.pkey2
            fftype_i =  self.strucC.particles[ pid_i ].properties["fftype"]
            fftype_j =  self.strucC.particles[ pid_j ].properties["fftype"]
            #r_i = np.array( self.strucC.particles[ bondObj_o.pkey1 ].position  )
            #r_j = np.array( self.strucC.particles[ bondObj_o.pkey2 ].position  )
            #bond_len = np.linalg.norm(delta_r_c(r_i,r_j,struc_o.getLatVec() ))

            for btyp_p, btypObj_p  in self.paramC.bondtypes.iteritems():
                lmptype_p = btypObj_p.lmpindx
                p_i = btypObj_p.fftype1 
                p_j = btypObj_p.fftype2
                match = False 
                if( fftype_i == p_i  and  fftype_j == p_j ):
                    match = True
                elif( fftype_i == p_j  and  fftype_j == p_i ):
                    match = True
                if( match ):
                    new_type = False
                    bondObj_o.lmpindx = btypObj_p.lmpindx
                    bondObj_o.g_indx = btypObj_p.g_indx
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
                
                for b_all, btypObj_all  in paramC_o.bondtypes.iteritems():
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
                    raise TypeError(" No Bond parameters were found for bond type %s-%s "%(fftype_i,fftype_j))
                elif( cnt_check > 1 ):
                    logger.warning(" %d  Bond parameters were found for bond type %s-%s "%(cnt_check,fftype_i,fftype_j))
                    for btyp_p, btypObj_p  in self.paramC.bondtypes.iteritems():
                        print btyp_p ,btypObj_p.fftype1 ,btypObj_p.fftype2                    
                    if( not use_last  ):
                        raise TypeError("Last parameter will not be used")
                        
                if( type_found ):
                        
                        bondObj_o.g_indx = bondObj_temp.g_indx
                        # Set lammps type for bond and bond type
                        # These have to match!!! 
                        bondObj_o.lmpindx = lmptype_p          # Set LAMMPS index for bond
                        bondObj_temp.lmpindx = lmptype_p
                        self.paramC.add_bondtype(bondObj_temp)
                        if( debug ):
                            print " %d  Bond parameters were found for bond type %s-%s "%(cnt_check,fftype_i,fftype_j)

                        # log_line=" Setting bond atoms %s - %s numbers %d - %d wiht bond length %f to type %d with r_o %f  delta %f \n"%(fftype_i,fftype_j,pid_i,pid_j,bond_len,btyp_p,btypObj_all.get_r0(),bond_len-btypObj_all.get_r0() )
                        log_line="Adding new lmptyp %d for bond atoms %s - %s numbers %d - %d "%(lmptype_p,fftype_i,fftype_j,pid_i,pid_j)
                        param_out.write(log_line+'\n')


        #
        # Examine  angles types
        #
        debug = 0
        # for akey_i, angle_i in self.strucC.angles.iteritems():
        for a_o,angleObj_o in self.strucC.angles.iteritems():
            new_type = True
            lmptype_p = 0
            pid_k = angleObj_o.pkey1
            pid_i = angleObj_o.pkey2
            pid_j = angleObj_o.pkey3
            fftype_k =  self.strucC.particles[ angleObj_o.pkey1 ].properties["fftype"]
            fftype_i =  self.strucC.particles[ angleObj_o.pkey2 ].properties["fftype"]
            fftype_j =  self.strucC.particles[ angleObj_o.pkey3 ].properties["fftype"]
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
                    angleObj_o.lmpindx = lmptype_p
                    angleObj_o.g_indx = atypObj_p.g_indx

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
                for a_all, atypObj_all  in paramC_o.angletypes.iteritems():
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
                    raise TypeError(" No Angles parameters were found for bond type %s-%s-%s "%(fftype_k,fftype_i,fftype_j))
                elif( cnt_check > 1 ):
                    # log_line=" %d Angles parameters were found for angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f  \n"%(cnt_check,fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij )
                    logger.warning(" %d Angles parameters were found for angle atoms %s - %s - %s numbers %d - %d - %d   \n"%(cnt_check,fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j ))
                    #param_out.write(log_line)
                    print log_line
                    # atypC_p.findtype(fftype_k,fftype_i,fftype_j)

                    if( not use_last  ):
                        raise TypeError("Last parameter will not be used")
                if( type_found ):
                        
                        angleObj_o.g_indx = atypObj_temp.g_indx
                        angleObj_o.lmpindx = lmptype_p
                        atypObj_temp.lmpindx = lmptype_p
                        self.paramC.add_angletype(atypObj_temp)
                        if( debug ):
                            print " %d Angles parameters were found for bond type %s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j)
                            #log_line=" Setting angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f to type %d with theta_o %f  delta %f \n"%(fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij,atyp_p+1,theta0_kij,delta_theta )
                        log_line=" Adding new lmptyp %d  angle atoms %s - %s - %s numbers %d - %d - %d "%(lmptype_p,fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j)
                        param_out.write(log_line+'\n')


        #
        # Examine  dihedrals types
        #
        debug = False
        if( debug):
            for d_all, dtypObj_all  in paramC_o.dihtypes.iteritems():
                all_k = dtypObj_all.fftype1 
                all_i = dtypObj_all.fftype2 
                all_j = dtypObj_all.fftype3
                all_l = dtypObj_all.fftype4
                print " all types in parameters ",all_k,all_i,all_j,all_l,dtypObj_all.type
                #sys.exit("  type check debug ")

        imp_cnt = 0 
        for d_o,dihObj_o in self.strucC.dihedrals.iteritems():
            new_type = True
            lmptype_p = 0
            pid_k = dihObj_o.pkey1
            pid_i = dihObj_o.pkey2 
            pid_j = dihObj_o.pkey3
            pid_l = dihObj_o.pkey4
            fftype_k =  self.strucC.particles[pid_k].properties["fftype"]
            fftype_i =  self.strucC.particles[pid_i].properties["fftype"]
            fftype_j =  self.strucC.particles[pid_j].properties["fftype"]
            fftype_l =  self.strucC.particles[pid_l].properties["fftype"]

            if( debug):
                print " checking ",fftype_k, fftype_i,  fftype_j , fftype_l
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
                    dihObj_o.lmpindx = lmptype_p
                    dihObj_o.g_indx = dtypObj_p.g_indx

                    if( debug):
                        print " dihObj_o.lmpindx  ",dihObj_o.lmpindx
                        print " dihObj_o.g_indx  ",dihObj_o.g_indx
                        print "  previous type ",dtyp_p,p_k,p_i,p_j,p_l,dihObj_o.g_indx
                    break
                    

            # If it is not in the parameter set for the struture container
            #  find it in the parameters from the reference parameter file 
            if( new_type ):
                # Find type in btypC_all
                lmptype_p += 1
                cnt_check = 0
                type_found = False 
                # Set type to new type = last type+1
                dihObj_o.lmpindx = lmptype_p

                if( debug):
                    print "  new type checking against %d read in parameters "%len(paramC_o.dihtypes)

                copy_type = False 
                for d_all, dtypObj_all  in paramC_o.dihtypes.iteritems():
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
                        #dtypObj_temp.set_g_indx(dtypObj_all.g_indx)
                        type_found = True
                        copy_type = False 
                        if( debug ):
                            print " %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                            print "     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l

                if( not type_found ):
                    if(debug):
                        print " checking  X - FF - FF - FF "
                    copy_type = False 
                    for d_all, dtypObj_all  in paramC_o.dihtypes.iteritems():
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
                            #dtypObj_temp.set_g_indx(dtypObj_all.g_indx)
                            type_found = True 
                            copy_type = False 
                            if( debug ):
                                print " %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                                print "     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l

                if( not type_found ):
                    if(debug):
                        print " checking  X - FF - FF - X "
                    copy_type = False 
                    for d_all, dtypObj_all  in paramC_o.dihtypes.iteritems():
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
                            #dtypObj_temp.set_g_indx(dtypObj_all.g_indx)
                            type_found = True 
                            copy_type = False 
                            if( debug ):
                                print " %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                                print "     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l

                if( cnt_check < 1 ):
                    raise TypeError(" No Dih parameters were found for dih type %s-%s-%s-%s "%(fftype_k,fftype_i,fftype_j,fftype_l))
                elif( cnt_check > 1 ):
                    print " %d Dih parameters were found for dih type %s-%s-%s-%s please check parameter file  "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                    print dtypObj_temp
                    #dtypObj_temp_list.findtype(fftype_k,fftype_i,fftype_j,fftype_l)
                    if( not use_last  ):
                        raise TypeError

                if( type_found ):
                    if( debug ):

                        print " adding new type to dtypC_p  from ",dtypObj_temp.type, dtypObj_temp.fftype1,dtypObj_temp.fftype2,dtypObj_temp.fftype3,dtypObj_temp.fftype4

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
                            print " Normalizing dihedral potential "
                            print " finding types for ",pid_i,pid_j
                        NNAB_i = calc_nnab(pid_i,cov_nbindx) - 1
                        NNAB_j = calc_nnab(pid_j,cov_nbindx) - 1

                        dihen_norm = float( NNAB_i + NNAB_j)/2.0

                        if(debug): print " dihen_norm ",dihen_norm

                        dtypObj_temp.normforceconstants(dihen_norm)
                    '''

                    dihObj_o.g_indx = dtypObj_temp.g_indx
                    dtypObj_temp.lmpindx = lmptype_p
                    self.paramC.add_dihtype(dtypObj_temp)                
                    if( debug ):
                        print " %d Dih parameters were found for dih type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                        print " len(dtypC_p) ",len(self.paramC.dihtypes) 
                        print " dtypObj_temp.lmpindx  ",dtypObj_temp.lmpindx
                        print " dtypObj_temp.g_indx  ",dtypObj_temp.g_indx
        #
        # Examine improper dihedrals types
        #
        debug = False
        if( debug):
            for d_all, imptypObj_all  in paramC_o.imptypes.iteritems():
                all_k = imptypObj_all.fftype1 
                all_i = imptypObj_all.fftype2 
                all_j = imptypObj_all.fftype3
                all_l = imptypObj_all.fftype4
                print " all types in parameters ",all_k,all_i,all_j,all_l,imptypObj_all.type
                #sys.exit("  type check debug ")

        imp_cnt = 0 
        for imp_o,impObj_o in self.strucC.impropers.iteritems():
            new_type = True
            lmptype_p = 0
            pid_k = impObj_o.pkey1
            pid_i = impObj_o.pkey2 
            pid_j = impObj_o.pkey3
            pid_l = impObj_o.pkey4 
            fftype_k =  self.strucC.particles[pid_k].properties["fftype"]
            fftype_i =  self.strucC.particles[pid_i].properties["fftype"]
            fftype_j =  self.strucC.particles[pid_j].properties["fftype"]
            fftype_l =  self.strucC.particles[pid_l].properties["fftype"]

            if( debug):
                print " checking ",fftype_k, fftype_i,  fftype_j , fftype_l
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
                    impObj_o.lmpindx = lmptype_p
                    impObj_o.g_indx = imptypObj_p.g_indx

                    if( debug):
                        print " impObj_o.lmpindx  ",impObj_o.lmpindx
                        print " impObj_o.g_indx  ",impObj_o.g_indx
                        print "  previous type ",imptyp_p,p_k,p_i,p_j,p_l,impObj_o.g_indx
                    break
                    
            # If it is not in the parameter set for the struture container
            #  find it in the parameters from the reference parameter file 
            if( new_type ):
                # Find type in btypC_all
                lmptype_p += 1
                cnt_check = 0
                type_found = False 
                # Set type to new type = last type+1
                impObj_o.lmpindx = lmptype_p

                if( debug):
                    print "  new type checking against %d read in parameters "%len(imptypC_all)

                copy_type = False 
                for d_all, imptypObj_all  in paramC_o.imptypes.iteritems():
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
                        #imptypObj_temp.set_g_indx(imptypObj_all.g_indx)
                        type_found = True
                        copy_type = False 
                        if( debug ):
                            print " %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                            print "     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l

                if( not type_found ):
                    if(debug):
                        print " checking  X - FF - FF - FF "
                    copy_type = False 
                    for d_all, imptypObj_all  in paramC_o.imptypes.iteritems():
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
                            #imptypObj_temp.set_g_indx(imptypObj_all.g_indx)
                            type_found = True 
                            copy_type = False 
                            if( debug ):
                                print " %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                                print "     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l

                if( not type_found ):
                    if(debug):
                        print " checking  X - FF - FF - X "
                    copy_type = False 
                    for d_all, imptypObj_all  in paramC_o.imptypes.iteritems():
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
                            #imptypObj_temp.set_g_indx(imptypObj_all.g_indx)
                            type_found = True 
                            copy_type = False 
                            if( debug ):
                                print " %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                                print "     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l

                if( cnt_check < 1 ):
                    raise TypeError(" No Dih parameters were found for dih type %s-%s-%s-%s "%(fftype_k,fftype_i,fftype_j,fftype_l))
                elif( cnt_check > 1 ):
                    print " %d Dih parameters were found for dih type %s-%s-%s-%s please check parameter file  "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                    print imptypObj_temp
                    #imptypObj_temp_list.findtype(fftype_k,fftype_i,fftype_j,fftype_l)
                    if( not use_last  ):
                        raise TypeError

                if( type_found ):
                    if( debug ):

                        print " adding new type to imptypC_p  from ",imptypObj_temp.type, imptypObj_temp.fftype1,imptypObj_temp.fftype2,imptypObj_temp.fftype3,imptypObj_temp.fftype4

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
                            print " Normalizing dihedral potential "
                            print " finding types for ",pid_i,pid_j
                        NNAB_i = calc_nnab(pid_i,cov_nbindx) - 1
                        NNAB_j = calc_nnab(pid_j,cov_nbindx) - 1

                        dihen_norm = float( NNAB_i + NNAB_j)/2.0

                        if(debug): print " dihen_norm ",dihen_norm
                        print " dihen_norm ",dihen_norm

                        imptypObj_temp.normforceconstants(dihen_norm)
                    '''
                    impObj_o.g_indx = imptypObj_temp.g_indx
                    imptypObj_temp.lmpindx = lmptype_p
                    self.paramC.add_imptype(imptypObj_temp,deepcopy = True)
                    if( debug ):
                        print " %d Dih parameters were found for dih type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                        print " len(imptypC_p) ",len(imptypC_p) 
                        print " imptypObj_temp.lmpindx  ",imptypObj_temp.lmpindx
                        print " imptypObj_temp.g_indx  ",imptypObj_temp.g_indx

        debug = False 
        if(debug):
            print " LJ atom types found %d "%(self.n_ljtypes)
            for lj_p, ljObj_p  in self.paramC.ljtypes.iteritems(): 
                print lj_p,ljObj_p.fftype1,ljObj_p.mass,ljObj_p.epsilon,ljObj_p.sigma
            print " Bond types found %d "%(self.n_bondtypes)
            for btyp_p, btypObj_p  in self.paramC.bondtypes.iteritems():
                print btyp_p ,btypObj_p.fftype1 ,btypObj_p.fftype2,btypObj_p.lmpindx,btypObj_p.g_indx
            print " Angle types found %d "%(self.n_angletypes)
            for atyp_p, atypObj_p  in self.paramC.angletypes.iteritems():
                print atyp_p ,atypObj_p.fftype1 ,atypObj_p.fftype2,atypObj_p.fftype3,atypObj_p.lmpindx,atypObj_p.g_indx
            print " Dih types found %d "%(self.n_dihtypes)
            for dtyp_p, dtypObj_p  in self.paramC.dihtypes.iteritems():
                print dtyp_p ,dtypObj_p.fftype1 ,dtypObj_p.fftype2,dtypObj_p.fftype3,dtypObj_p.fftype4,dtypObj_p.lmpindx,dtypObj_p.g_indx
            print " imp Dih types found %d "%(self.n_imptypes)
            for imptyp_p, dtypObj_p  in self.paramC.imptypes.iteritems():
                print imptyp_p ,dtypObj_p.fftype1 ,dtypObj_p.fftype2,dtypObj_p.fftype3,dtypObj_p.fftype4,dtypObj_p.lmpindx,dtypObj_p.g_indx
            sys.exit('find_types')

        debug = False 
        if(debug):
            print "  All particles should have new type labeled as interger stored as a string "
            for pkey_o, particle_o  in self.strucC.particles.iteritems():
                print particle_o.properties["fftype"],particle_o.type
            for d_o,dihObj_o in dihC_o:
                print " lmpindx() g_indx()  ",d_o,dihObj_o.pkey1,dihObj_o.pkey2,dihObj_o.pkey3,dihObj_o.pkey4, dihObj_o.lmpindx ,dihObj_o.g_indx 

        param_out.close()

        return           
        
class CalculationRes(Calculation):
    '''
    Derived type of Calculation for running a calculation on resource.
    In that files will be compressed and moved to scratch directories to run,
    then output files and output data will be compressed and moved to storage.

    '''
    def __init__(self, tag , verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """
        # Base class constructor is called
        Calculation.__init__(self, tag)
        # Computational Resource used for simulation/calculation  
        self.resource = resource.Resource()
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

        print json_data
        
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
                    print "No resource found "
                if( len(res_tag) > 0 ):
                    print "Resource tag found %s "%(res_tag)
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
                        print " Need to set reference calculation type "
                except:
                    print "No references found "
                    


        except IOError:
            logger.warning(" File not found %s in %s "%(json_file,os.getcwd()))
            print " File not found %s in %s "%(json_file,os.getcwd())

        
    def make_dir(self):
        '''
        Check that needed directories exist 
        '''
        logger.debug("Creating directories for resource %s "%(self.tag))
        
        if( self.resource.meta['type'] == "local" ):
            os.chdir(self.dir['home'])
            for dkey,dir_i in self.dir.iteritems():
                if ( not os.path.isdir(dir_i) ):
                    print "Making %s "%(dir_i)
                    os.mkdir(dir_i)
            os.chdir(self.dir['home'])
        elif( self.resource.meta['type'] == "ssh" ):
            # Create local directories 
            os.chdir(self.dir['home'])
            dkey = 'launch'
            try:
                dir_i = self.dir[dkey]
                if ( not os.path.isdir(dir_i) ):
                    print "Making %s "%(dir_i)
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
            print " dir dictionary key %s not set for calculation %s files "%(from_dirkey,self.tag)
        try:
            todir = self.dir[to_dirkey]
        except:
            cpfile = False
            print " dir dictionary key %s not set for calculation %s files "%(to_dirkey,self.tag)
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
        print "> in cp_file ",file_type,file_key,file_name,from_dirkey,to_dirkey
        
        cpfile = True
        try:
            fromdir = self.dir[from_dirkey]
        except:
            cpfile = False
            print " dir dictionary key %s not set for calculation %s files "%(from_dirkey,self.tag)
        try:
            todir = self.dir[to_dirkey]
        except:
            cpfile = False
            print " dir dictionary key %s not set for calculation %s files "%(to_dirkey,self.tag)
        if( cpfile ):
            from_pathfile = os.path.join(fromdir,file_name)
            to_pathfile = os.path.join(todir,file_name)
            print "copying %s to %s "%(from_pathfile,to_pathfile)
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


    def push(self):
        '''
        Push input files to resource 
        '''
        #
        print " Resource type %s "%(self.resource.meta['type'])

        if( self.meta['status'] == 'written' ):
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
                print "runnning push function in %s "%(os.getcwd())
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
                        print "Compressing and copying %s files to scratch directory "%(file_type)
                        self.compress_files(file_type)
                        file_name = self.files[file_type][file_key]
                        # self.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)
                        from_pathfile = os.path.join(fromdir,file_name)
                        to_pathfile = os.path.join(todir,file_name)
                        bash_command = "scp %s %s:%s"%(from_pathfile,ssh_id,todir)
                        os.system(bash_command)
                        # Uncompress file 
                        bash_command = "%s %s "%(self.properties['uncompress'],file_name)
                        bash_command = 'ssh %s \' cd %s ; %s \' '%(ssh_id,todir,bash_command)
                        os.system(bash_command)
                    else:
                        print "No files of type %s present"%(file_type)
            else:
                print " Resource type %s does not need to push files created localy "%(self.resource.meta['type'])
            #
            # Copy reference simulation output over to scratch directory
            #
            file_key = self.properties['comp_key']
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])
            for ref_key,ref_calc in self.references.iteritems():
                print "Copying output of reference calculations %s"%(ref_key)
                file_type = 'output'
                try:
                    file_name = ref_calc.files[file_type][file_key]
                except:
                    print "Fine type %s has no compressed files "%(file_type)
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
                    else:
                        print " Copy from  type %s to type %s not set  "%(ref_calc.resource.meta['type'], self.resource.meta['type'])
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
        
        
        
            
        
        
