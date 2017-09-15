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
    def __init__(self,fftype1='X',unit_conf=units.unit_conf):
        self.unit_conf = unit_conf
        
        self.fftype1 = str(fftype1) 
        self.epsilon = 1.0
        self.sigma = 2.0 
        self.lammps_index = None
        self.gromacs_index = None 
    
    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.fftype1
        del self.epsilon
        del self.sigma 
        del self.lammps_index
        del self.gromacs_index
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " {} epsilon:{} sigma:{}".format(self.fftype1,self.epsilon,self.sigma)
