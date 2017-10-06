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

class Particletype(units.ObjectUnits):
    '''
    Particle represented by a Force-field
    
    Kwargs:
        * fftype1 (str): Forcefield key 
        * units_conf (dict): Dictionary of units for each attribute type
                
    '''
    
    @property
    def mass(self):
        return self._property['mass'] 


    @mass.setter
    def mass(self,value):
        self._property['mass'] = value

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
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        #
        self.fftype1 = str(fftype1)

        #         
        self._property['mass'] = 1.0
        self._property['epsilon'] = 1.0
        self._property['sigma']  = 2.0
        
        self._property_units['mass'].append('mass')
        self._property_units['energy'].append('epsilon')
        self._property_units['length'].append('sigma')
        
        
        self.lammps_index = None
        self.gromacs_index = None 
    
    def __del__(self):
        del self.fftype1
        del self.lammps_index
        del self.gromacs_index
        
    def __str__(self):
        return " {} epsilon:{} sigma:{}".format(self.fftype1,self.epsilon,self.sigma)

