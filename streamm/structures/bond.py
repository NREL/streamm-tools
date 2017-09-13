#! /usr/bin/env python
"""
This module defines the classes relating to bonds between atoms
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

class Bond(object):
    """Data structure for describing any 2-point association of Particles

    Args:
        pkey1 (int): Index of Particle i.
        pkey2 (int): Index of Particle i.
        

    .. attribute:: length (float)

        length of the bond 
                
    """
    def __init__(self, pkey1, pkey2):
        self.pkey1 = pkey1
        self.pkey2 = pkey2
        
        self.length = None 

    def __del__(self):
        del self.pkey1
        del self.pkey2
        del self.length 

    def __str__(self):
        return " %s - %s"%(self.pkey1,self.pkey2 )

        