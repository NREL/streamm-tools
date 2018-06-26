# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

"""
This module defines the classes relating to bonds between atoms
"""

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"


import logging
logger = logging.getLogger(__name__)


# Import pymatgen module 
import pymatgen_core.core.units as units 


class Bond(units.ObjectUnits):
    """Data structure for describing any 2-point association of Particles

    Args:
        * pkey1 (int): Index of Particle i.
        * pkey2 (int): Index of Particle i.
        
    Kwargs:
        * units_conf (dict): Dictionary of units for each attribute type
        
    .. attribute:: length (float)

        length of the bond
        unit type (length)
                
    """

    @property
    def length(self):
        return self._property['length']

    @length.setter
    def length(self,value):
        self._property['length'] = value
        
    @property
    def bondorder(self):
        return self._property['bondorder']

    @bondorder.setter
    def bondorder(self,value):
        self._property['bondorder'] = value
        
          
    def __init__(self, pkey1, pkey2,unit_conf=units.unit_conf ):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
                
        self.pkey1 = pkey1
        self.pkey2 = pkey2
        
        self._property['length'] = 0.0 
        self._property_units['length'].append('length')

        self._property['bondorder'] = 1
        
        # Force field
        self.param = None
        # Lammps and gromacs index
        self.param_index = 0 
        self.lammps_index = 0 
        self.gromacs_index = 0 


        #for bkey,bond_i in strucC.bonds.iteritems():
        #    print ">set_bondorder ",bkey,bond_i.pkey1,bond_i.pkey2,bond_i.properties['bondorder'] 
    #
    
    def __del__(self):
        del self.pkey1
        del self.pkey2
        del self.param
        del self.param_index
        del self.lammps_index
        del self.gromacs_index
        
    def __str__(self):
        return " %s - %s"%(self.pkey1,self.pkey2 )


    def export_json(self):
        '''    
        Export object to json
        
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        
        json_data = {}
        json_data['pkey1'] = self.pkey1
        json_data['pkey2'] = self.pkey2
        json_data['length'] = self.length
        json_data['param_index'] = self.param_index
        json_data['lammps_index'] = self.lammps_index
        json_data['gromacs_index'] = self.gromacs_index
        #
        json_data['bondorder'] = self.bondorder
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
        self.length  =  json_data['length']
        self.param_index  =  json_data['param_index']
        self.lammps_index  =  json_data['lammps_index']
        self.gromacs_index  =  json_data['gromacs_index']
        self.bondorder  =  json_data['bondorder']



