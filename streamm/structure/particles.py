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

import logging
logger = logging.getLogger(__name__)

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
        self.mass = 0.0
        self.charge = 0.0

        
        self.mol = 0
        self.ring = 0
        self.residue = 0
        self.resname = "RES"
        self.qgroup = 0
        
    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.type
        del self.tag
        del self.mass
        del self.charge
        # 
        del self.mol
        del self.ring
        del self.residue
        del self.resname
        del self.qgroup
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " {}".format(self.type)


