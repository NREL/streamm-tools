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
parameters module

Class data structures force-field parameters

Most molecular dynamics codes set each particle to certain type (fftype)
usually each interaction is set using these fftypes

"""


import logging
logger = logging.getLogger(__name__)



import pickle
import copy
import json

try:
    # Import pymatgen Class 
    import pymatgen_core.core.units as units
except:
    raise ImportError("pymatgen import error for pymatgen.core.units module")

# Import streamm dependencies 
from streamm.forcefields.particletype import Particletype
from streamm.forcefields.bondtype import Bondtype 
from streamm.forcefields.angletype import Angletype
from streamm.forcefields.dihtype import Dihtype
from streamm.forcefields.imptype import Imptype


def read_pickle(tag):
    '''    
    Pickle object
    '''
    with open("%s.pkl"%(tag),'rb') as fl:
        return pickle.load( fl )
                

class Parameters(units.ObjectUnits):
    """
    Container for force-field parameters
    

    Kwargs:
        * fftype1 (str): Forcefield key 
        * units_conf (dict): Dictionary of units for each attribute type
                
    """
    

    # 
    @property
    def n_particletypes(self):
        return len(self.particletypes)
    @property
    def n_bondtypes(self):
        return len(self.bondtypes)
    @property
    def n_angletypes(self):
        return len(self.angletypes)
    @property
    def n_dihtypes(self):
        return len(self.dihtypes)
    @property
    def n_imptypes(self):
        return len(self.imptypes)
            

    def __init__(self,tag='blank',unit_conf=units.unit_conf):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        
        self.tag = tag
        self.sufix = 'param'
        # 
        self.particletypes = dict()                           # Creates empty dict struc
        self.bondtypes = dict()                                   # Creates empty dict struc
        self.angletypes = dict()                                # Creates empty dict struc
        self.dihtypes = dict()                                # Creates empty dict struc
        self.imptypes = dict()                                  # Creates empty dict struc
        #     
        #
        # Set defaults 
        #
        self.nbfunc = 1      # Use 1 (Lennard-Jones) or 2 (Buckingham)
        self.combmixrule = 3   # 1 or 3  - geometric; 2 - arithmetic
        self.genpairs = "yes"  # generates 1-4 parameters that are not present in the pair list 
        self.fudgeLJ = 1.0   # the factor by which to multiply Lennard-Jones 1-4 interactions, default 1
        self.fudgeQQ = 1.0   # the factor by which to multiply electrostatic 1-4 interactions, default 1
        

    def __del__(self):
        """
        Destructor, clears object memory
        """
        #
        del self.tag 
        del self.sufix 
        del self.particletypes
        del self.bondtypes
        del self.angletypes
        del self.dihtypes
        del self.imptypes
        # 
        del self.nbfunc
        del self.combmixrule
        del self.genpairs
        del self.fudgeLJ
        del self.fudgeQQ
        
    def __str__(self):
        """
        'Magic' method for printng contents
        """

        strucStr =  "\n"
        strucStr += "    Parameters \n"
        strucStr += "      LJ parameters %d \n"%(self.n_particletypes)
        strucStr += "      Bond parameters %d \n"%(self.n_bondtypes)
        strucStr += "      Angle parameters %d \n"%(self.n_angletypes)
        strucStr += "      Dihedral parameters %d \n"%(self.n_dihtypes)
        strucStr += "      Improper Dihedral parameters %d \n"%(self.n_imptypes)
        return strucStr


    def dump_pickle(self):
        '''    
        Pickle object
        '''
        file_i = open("%s.pkl"%(self.tag),'w')
        pickle.dump(self,file_i)
        file_i.close()
        

        
    def add_particletype(self, particletype_i, deepcopy = True ):
        """
        Add 'Ljtype' object to ljtypes dict in this container and update n_ljtypes accordingly
        """
        if isinstance(particletype_i, Particletype):
            index = self.n_particletypes
            if( deepcopy ):
                self.particletypes[index] = copy.deepcopy(particletype_i) # index 0 -> (N-1)
            else:
                self.particletypes[index] = particletype_i # index 0 -> (N-1)
                
        else:
            raise TypeError("Attempting to add non-paticletype type to container")


    def add_bondtype(self, bondtype_i, deepcopy = True ):
        """
        Add 'Bondtype' object to bondtypes dict in this container and update n_bondtypes accordingly
        """
        if isinstance(bondtype_i, Bondtype):
            index = self.n_bondtypes
            if( deepcopy ):
                self.bondtypes[index] = copy.deepcopy(bondtype_i) # index 0 -> (N-1)
            else:
                self.bondtypes[index] = bondtype_i # index 0 -> (N-1)
        else:
            print "Attempting to add non-Bondtype type to container"
            raise TypeError


    def add_angletype(self, angletype_i, deepcopy = True ):
        """
        Add 'Angletype' object to angletypes dict in this container and update n_angletypes accordingly
        """
        if isinstance(angletype_i, Angletype):
            index = self.n_angletypes
            if( deepcopy ):
                self.angletypes[index] = copy.deepcopy(angletype_i) # index 0 -> (N-1)
            else:
                self.angletypes[index] = angletype_i # index 0 -> (N-1)
        else:
            print "Attempting to add non-Angletype type to container"
            raise TypeError

    def add_dihtype(self, dihtype_i, deepcopy = True ):
        """
        Add 'Dihtype' object to dihtypes dict in this container and update n_dihtypes accordingly
        """
        if isinstance(dihtype_i, Dihtype):
            index = self.n_dihtypes
            if( deepcopy ):
                self.dihtypes[index] = copy.deepcopy(dihtype_i) # index 0 -> (N-1)
            else:
                self.dihtypes[index] = dihtype_i # index 0 -> (N-1)
        else:
            print "Attempting to add non-Dihtype type to container"
            raise TypeError


    def add_imptype(self, imptype_i, deepcopy = True ):
        """
        Add 'Imptype' object to imptypes dict in this container and update n_imptypes accordingly
        """
        if isinstance(imptype_i, Imptype):
            index = self.n_imptypes
            if( deepcopy ):
                self.imptypes[index] = copy.deepcopy(imptype_i) # index 0 -> (N-1)
            else:
                self.imptypes[index] = imptype_i # index 0 -> (N-1)
        else:
            print "Attempting to add non-Imptype type to container"
            raise TypeError

    def __iadd__(self, other ):
        """
        'Magic' method to implement the '+=' operator 
        
        """

        if isinstance(other, Container):
            for ptkey_i, particletypes_i  in other.particletypes.iteritems():
                self.add_particletype(particletypes_i,deepcopy = True )
            for btkey_i,bondtype_i  in other.bondtypes.iteritems():
                self.add_bondtype(bondtype_i,deepcopy = True )
            for atkey_i,angletype_i  in other.angletypes.iteritems():
                self.add_angletype(angletype_i,deepcopy = True )
            for dtkey_i, dihtype_i  in other.dihtypes.iteritems():    
                self.add_dihtype(dihtype_i,deepcopy = True )
            for itkey_i, imptype_i  in other.imptypes.iteritems():    
                self.add_imptype(imptype_i,deepcopy = True )


        else:
            raise TypeError("iadd object should be a Parameter.Container object")
        

        return self


    def update_units(self,new_unit_conf):
        '''
        Update instance values with new units
        
        Args:
            * new_unit_conf (dict): with unit type as the key and the new unit as the value
            
        '''
        # 
        for ptkey_i, particletypes_i  in self.particletypes.iteritems():
            particletypes_i.update_units(new_unit_conf)
        
        for btkey_i,bondtype_i  in self.bondtypes.iteritems():
            bondtype_i.update_units(new_unit_conf)
            
        for atkey_i,angletype_i  in self.angletypes.iteritems():
            angletype_i.update_units(new_unit_conf)
            
        for dtkey_i, dihtype_i  in self.dihtypes.iteritems():    
            dihtype_i.update_units(new_unit_conf)
            
        for itkey_i, imptype_i  in self.imptypes.iteritems():    
            imptype_i.update_units(new_unit_conf)
        
        

    def export_json(self,write_file=True):
        '''    
        Export particles to json
        
        Kwargs:
            * write_file (boolean) to dump json to a file
            
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        #
        json_data = {}
        # Properties
        json_data['nbfunc'] = self.nbfunc
        json_data['combmixrule'] = self.combmixrule
        json_data['genpairs'] = self.genpairs
        json_data['fudgeLJ'] = self.fudgeLJ
        json_data['fudgeQQ'] = self.fudgeQQ
        # particles
        json_data['particletypes']  = {}
        for pk,p in self.particletypes.iteritems():
            json_data['particletypes'][pk] = p.export_json()
        # bonds
        json_data['bondtypes']  = {}
        for bk,b in self.bondtypes.iteritems():
            json_data['bondtypes'][bk] = b.export_json()
        # angles
        json_data['angletypes']  = {}
        for ak,a in self.angletypes.iteritems():
            json_data['angletypes'][ak] = a.export_json()
        # dihedrals
        json_data['dihtypes']  = {}
        for dk,d in self.dihtypes.iteritems():
            json_data['dihtypes'][dk] = d.export_json()
        # impropers
        json_data['imptypes']  = {}
        for ik,i in self.imptypes.iteritems():
            json_data['imptypes'][ik] = i.export_json()
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
                
        self.nbfunc  =  json_data['nbfunc']
        self.combmixrule  =  json_data['combmixrule']
        self.genpairs  =  json_data['genpairs']
        self.fudgeLJ  =  json_data['fudgeLJ']
        self.fudgeQQ  =  json_data['fudgeQQ']

        if( 'particletypes' in json_data.keys() ):
            for pk,json_particletype in sorted(json_data['particletypes'].iteritems()):
                pk = int(pk)
                particletype_i = Particletype()
                particletype_i.import_json(json_particletype)
                self.particletypes[pk] = copy.deepcopy(particletype_i)
                
        if( 'bondtypes' in json_data.keys() ):
            for bk,json_bondtype in json_data['bondtypes'].iteritems():
                bk = int(bk)
                b = Bondtype(json_bondtype['fftype1'],json_bondtype['fftype2'])
                b.import_json(json_bondtype)
                self.bondtypes[bk] = copy.deepcopy(b)
                
        if( 'angletypes' in json_data.keys() ):
            for ak,json_angletype in json_data['angletypes'].iteritems():
                ak = int(ak)
                a = Angletype(json_angletype['fftype1'],json_angletype['fftype2'],json_angletype['fftype3'])
                a.import_json(json_angletype)
                self.angletypes[ak] = copy.deepcopy(a) 
                
        if( 'dihtypes' in json_data.keys() ):
            for dk,json_dihtype in json_data['dihtypes'].iteritems():
                dk = int(dk)
                d = Dihtype(json_dihtype['fftype1'],json_dihtype['fftype2'],json_dihtype['fftype3'],json_dihtype['fftype4'])
                d.import_json(json_dihtype)
                self.dihtypes[dk] = copy.deepcopy(d) 

        if( 'imptypes' in json_data.keys() ):
            for ik,json_imptype in json_data['imptypes'].iteritems():
                ik = int(ik)
                i = Imptype(json_imptype['fftype1'],json_imptype['fftype2'],json_imptype['fftype3'],json_imptype['fftype4'])
                i.import_json(json_imptype)
                self.imptypes[ik] = copy.deepcopy(i) 
            
        #
        return 
        
        