# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

"""
This module defines the classes relating to general particles 
"""

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

# Standard packages 
from numpy import pi #, dot, transpose, radians
import copy
import math
import random
import sys
import json
import csv
import os
try:
    import cPickle as pickle
except:
    import pickle
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

# Dependency packages 
import numpy as np
import pandas as pd

# pymatgen module
import pymatgen.core.periodic_table as periodictable


class Particle(object):
    """
    Data structure for describing any localized object in QM/MD simulations
    A 'Particle' has a type and dict of specifiers to it's property 
    """
    def __init__(self, type="blank"):
        """
        Constructor for a general particle. 
        """
        self.type = type
        self.tag = "blank"

    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.type
        del self.tag
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " {}".format(self.type)


class ForceField(object):
    '''
    Particle represented by a Force-field
    '''
    def __init__(self, type='X',label='X1'):
        """
        Constructor for a Force-field particle. 
        """
        self.type = type
        self.label = label
        self.charge = 0.0     
        self.mass  = 0.0
        self.lammps_index = -1 
        self.gromacs_index = -1 
    
    
    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.type
        del self.label
        del self.charge
        del self.mass
        del self.lammps_index
        del self.gromacs_index
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " {} ({})".format(self.label,self.type)

