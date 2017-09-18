# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

"""
This module defines the classes relating to bonds between atoms
"""

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

import logging
logger = logging.getLogger(__name__)


# Import pymatgen module 
import pymatgen_core.core.units as units 


class Bond(units.ObjectUnits):
    """Data structure for describing any 2-point association of Particles

    Args:
        pkey1 (int): Index of Particle i.
        pkey2 (int): Index of Particle i.
        
    Kwargs:
        units_conf (dict): Dictionary of units for each attribute type
        
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
        
          
    def __init__(self, pkey1, pkey2,unit_conf=units.unit_conf ):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        
        self.pkey1 = pkey1
        self.pkey2 = pkey2
        
        self._property['length'] = 0.0 
        self._property_units['length'].append('length')

    def __del__(self):
        del self.pkey1
        del self.pkey2
        
    def __str__(self):
        return " %s - %s"%(self.pkey1,self.pkey2 )
