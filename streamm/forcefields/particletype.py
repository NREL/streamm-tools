# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

# NoteTK  This should be in structure.particles 

try:
    # Import pymatgen Class 
    import pymatgen_core.core.units as units 
except:
    raise ImportError("pymatgen import error for units object")

class Particletype(object):
    '''
    Particle represented by a Force-field
    
    .. TODO ::
        change fftype1 to fftype_i
        
    '''


    @property
    def unit_conf(self):
        return self._unit_conf
    
    @property
    def epsilon(self):
        return self._property['epsilon'] 


    @epsilon.setter
    def epsilon(self,value):
        self._property['epsilon'] = value

    @property
    def sigma(self):
        return self._property['sigma'] 

    @sigma.setter
    def sigma(self,value):
        self._property['sigma'] = value


    def __init__(self,fftype1='X',unit_conf=units.unit_conf):
         
        self.fftype1 = str(fftype1)

        # Store the units of each attribute type 
        self._unit_conf = unit_conf
        #
        # Default Physical properties
        #
        self._property = {}
        self._property_units = {}
        for unit_type in self._unit_conf.keys():
            self._property_units[unit_type] = []
        #         
        self._property['epsilon'] = 1.0
        self._property['sigma']  = 2.0
        
        self._property_units['energy'].append('epsilon')
        self._property_units['length'].append('sigma')
        
        
        self.lammps_index = None
        self.gromacs_index = None 
    
    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.fftype1
        del self._property
        del self._property_units 
        del self.lammps_index
        del self.gromacs_index
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " {} epsilon:{} sigma:{}".format(self.fftype1,self.epsilon,self.sigma)


    def update_units(self,new_unit_conf):
        '''
        Update instance values with new units
        
        Args:
            new_unit_conf (dict): with unit type as the key and the new unit as the value
            
        '''
        
        self._property,self._unit_conf = units.change_properties_units(self._unit_conf,new_unit_conf,self._property_units,self._property)
        