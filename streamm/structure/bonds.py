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
            pkey1   (int)   Index of Particle object in bond
            pkey2   (int)   Index of Particle object in bond
        """
        
        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            raise TypeError("1st arg should be int")

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            raise TypeError("2nd arg should be int type")

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2

    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s"%(self.pkey1,self.pkey2 )


class ForceField(object):
    '''
    Particle represented by a Force-field
    '''
    def __init__(self, type='X',label='X1'):
        """
        Constructor for a Force-field particle. 
        """
        self.lammps_index = -1 
        self.gromacs_index = -1 
    
    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.lammps_index
        del self.gromacs_index
        