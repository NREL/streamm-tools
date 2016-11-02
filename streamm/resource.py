"""
Class data structures for resource data
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"



import copy , os , json , sys

import logging
logger = logging.getLogger(__name__)

class Resource:
    '''
    Data structure for a compute resource  
    '''
    def __init__(self,tag="blank",home_dir='',res_type = 'local'):
        '''
        Constructor for a general Resource object.
        
        Input:
            tag (str) unique id for project
            
        '''
        
        if isinstance(tag, str):
            self.tag = tag
        else:
            raise TypeError("1st arg (tag) in Resource initialization should be string name for the Resource")
        
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

        print "Reading in resource json file %s "%(json_file)
        
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
