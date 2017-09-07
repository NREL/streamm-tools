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
This module defines the classes relating to angles between atoms
"""
#
class Angle(object):
    """Data structure for describing any 3-point associatiaon of Particles

    Args:
        pkey1   (int):   Index of Particle object in angle
        pkey2   (int):   Index of Particle object in angle
        pkey3   (int):   Index of Particle object in angle    
    """
    def __init__(self, pkey1, pkey2, pkey3):
        self.pkey1 = pkey1
        self.pkey2 = pkey2
        self.pkey3 = pkey3
        
        #
        self.cosine = None 

    def __del__(self):
        del self.pkey1
        del self.pkey2
        del self.pkey3
        del self.cosine 


    def __str__(self):
        return " %s - %s - %s"%(self.pkey1,self.pkey2,self.pkey3)
