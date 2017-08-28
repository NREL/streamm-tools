# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"


class Particletype(object):
    '''
    Particle represented by a Force-field
    '''
    def __init__(self, fftype='X',label='X1'):
        """
        Constructor for a Force-field particle. 
        """
        self.fftype = fftype
        self.label = label
        self.charge = 0.0     
        self.mass  = 0.0
        self.epsilon = 0.0
        self.sigma = 0.0 
        self.lammps_index = -1 
        self.gromacs_index = -1
        # Particle key in structure.containers.Container 
        self.pkey = 0
    
    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.fftype
        del self.label
        del self.charge
        del self.mass
        del self.epsilon
        del self.sigma 
        del self.lammps_index
        del self.gromacs_index
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " {} ({})".format(self.fftype,self.label)
