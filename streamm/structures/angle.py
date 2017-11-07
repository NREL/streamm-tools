# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Dr. Travis W. Kemper"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3.3"
__email__ = "organicelectronics@nrel.gov"
__status__ = "Beta"


"""
This module defines the classes relating to angles between atoms
"""
import pymatgen_core.core.units as units
#
class Angle(units.ObjectUnits):
    """Data structure for describing any 3-point associations of Particles

    Args:
        * pkey1   (int):   Index of Particle object in angle
        * pkey2   (int):   Index of Particle object in angle
        * pkey3   (int):   Index of Particle object in angle
        
    Kwargs:
        * unit_conf (dict): Unit types with units used by this object
        
    .. attribute:: cosine (float)

        Cosine of the bond angle
        
    """
    
    def __init__(self, pkey1, pkey2, pkey3,unit_conf=units.unit_conf):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        self.pkey1 = pkey1
        self.pkey2 = pkey2
        self.pkey3 = pkey3
        
        #
        self.cosine = None 

        # Force field
        self.param = None
        # Lammps and gromacs index
        self.param_index = 0 
        self.lammps_index = 0 
        self.gromacs_index = 0             
            
    def __del__(self):
        del self.pkey1
        del self.pkey2
        del self.pkey3
        del self.cosine 
        # Force field 
        del self.param
        del self.param_index
        del self.lammps_index
        del self.gromacs_index

    def __str__(self):
        return " %s - %s - %s"%(self.pkey1,self.pkey2,self.pkey3)

    def export_json(self):
        '''    
        Export object to json
        
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        
        json_data = {}
        json_data['pkey1'] = self.pkey1
        json_data['pkey2'] = self.pkey2
        json_data['pkey3'] = self.pkey3
        json_data['cosine'] = self.cosine
        json_data['param_index'] = self.param_index
        json_data['lammps_index'] = self.lammps_index
        json_data['gromacs_index'] = self.gromacs_index
        #         
        return json_data
        

    def import_json(self,json_data):
        '''    
        Export object to json
        
        Args:
            * json_data (dict) json representation of the object
            
        '''
        self.pkey1  =  json_data['pkey1']
        self.pkey2  =  json_data['pkey2']
        self.pkey3  =  json_data['pkey3']
        self.cosine  =  json_data['cosine']
        self.param_index  =  json_data['param_index']
        self.lammps_index  =  json_data['lammps_index']
        self.gromacs_index  =  json_data['gromacs_index']

