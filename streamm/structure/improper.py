#! /usr/bin/env python
"""
This module defines the classes relating to the improper dihedral angles between atoms
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

class Improper(object):
    """
    Data structure for describing any 4-point associatiaon of Particle/s
    """

    def __init__(self, pkey1, pkey2, pkey3, pkey4):
        """
        Constructor for a general improper. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pkey1   (int)   Dictionary key of Particle object in improper
            pkey2   (int)   Dictionary key of Particle object in improper
            pkey3   (int)   Dictionary key of Particle object in improper
            pkey4   (int)   Dictionary key of Particle object in improper
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

        if isinstance(pkey4, int):
            self.pkey4 = pkey4
        else:
            raise TypeError("4rd arg should be int type")
        #
        self.lmpindx = 0
        self.g_indx = 0
        #
        self.properties = dict()

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2 
        del self.pkey3 
        del self.pkey4
        del self.lmpindx
        del self.g_indx 
        del self.properties

    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s - %s - %s "%(self.pkey1,self.pkey2,self.pkey3,self.pkey4 )
        
