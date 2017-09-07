# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

"""
This module defines the classes relating to the improper dihedral angles between atoms
"""

class Improper(object):
    """
    Data structure for describing any 4-point associatiaon of Particles
    
    Args:
        pkey1   (int)   Dictionary key of Particle object in improper
        pkey2   (int)   Dictionary key of Particle object in improper
        pkey3   (int)   Dictionary key of Particle object in improper
        pkey4   (int)   Dictionary key of Particle object in improper    
    """

    def __init__(self, pkey1, pkey2, pkey3, pkey4):

        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            raise TypeError("1st arg should be int")

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            raise TypeError("2nd arg should be int type")

        if isinstance(pkey3, int):
            self.pkey3 = pkey3
        else:
            raise TypeError("3rd arg should be int type")

        if isinstance(pkey4, int):
            self.pkey4 = pkey4
        else:
            raise TypeError("4rd arg should be int type")
        #
    def __del__(self):
        del self.pkey1
        del self.pkey2 
        del self.pkey3 
        del self.pkey4

    def __str__(self):
        return " %s - %s - %s - %s"%(self.pkey1,self.pkey2,self.pkey3,self.pkey4 )
        
