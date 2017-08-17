#! /usr/bin/env python
"""
This module defines the classes relating to bonds between atoms
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

class Bond(object):
    """
    Data structure for describing any 2-point association of Particle/s
    """
    def __init__(self, pkey1, pkey2):
        """
        Constructor for a general bond. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pkey1   (int)   Dictionary key of Particle object in bond
            pkey2   (int)   Dictionary key of Particle object in bond
        """
        
        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            raise TypeError("1st arg should be int")

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            raise TypeError("2nd arg should be int type")

        self.lmpindx = 0 
        self.g_indx = 0

        self.properties = dict()
        # self.border = 1 # 1-single,2-double,3-triple

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2
        del self.lmpindx
        del self.g_indx
        del self.properties

    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s  "%(self.pkey1,self.pkey2 )
