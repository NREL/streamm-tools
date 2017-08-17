#! /usr/bin/env python
"""
This module defines the classes relating to atoms
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
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

# Dependency packages 
import numpy as np
import pandas as pd

# pymatgen module
import pymatgen.core.periodic_table as periodictable

import logging
logger = logging.getLogger(__name__)


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
        # Tags dictionary. To be set by caller
        self.properties=dict()

    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.type
        del self.tag
        del self.properties
    
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s "%(self.type)

# Then you have types particles
class Atom(Particle):
    """
    A derived type of particles for atoms 
    """
    def __init__(self,symbol="X", type="atom"):
        '''
        Constructor for Particle object
        
        Args:
            symbol (str) Atomic symbol
            type   (str) Particle type
            NoteTK: type should probably not be set
        '''
        Particle.__init__(self, type=type)
        # 
        # Get properties of element based on symbol
        # 
        self.properties = periodictable.element_symbol(symbol)

        self.properties["mol"] = 0  
        self.properties["charge"] = 0.0     
        self.properties["fftype"] = self.properties["symbol"]
        self.properties["ffmass"] = self.properties["mass"]
        self.properties["lmpindx"] = -1 
        self.properties["group"] = 0
        self.properties["ring"] = 0
        self.properties["residue"] = 0
        self.properties["resname"] = "RES"
        self.properties["qgroup"] = 0
        self.properties["label"] =  self.properties["symbol"]
        
        self.tag = self.properties["symbol"]
        