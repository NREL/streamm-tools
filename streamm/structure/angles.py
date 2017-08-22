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
    """
    Data structure for describing any 3-point associatiaon of Particle/s
    """
    def __init__(self, pkey1, pkey2, pkey3):
        """
        Constructor for a general angle. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pkey1   (int)   Index of Particle object in angle
            pkey2   (int)   Index of Particle object in angle
            pkey3   (int)   Index of Particle object in angle
        """
        
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

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2
        del self.pkey3


    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s - %s"%(self.pkey1,self.pkey2,self.pkey3)
