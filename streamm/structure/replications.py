#! /usr/bin/env python
"""
This module defines the classes relating to replicating a structure container 
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

        
class Replication(object):
    '''
    Object to record the replication of structure 
    '''
    def __init__(self,name_i,name_j,name_ij,method,n):

        self.name_i = name_i
        self.name_j = name_j
        self.name_ij = name_ij
        self.method = method
        self.n = n
    
    def __del__(self):
        """
        'Magic' method for deleting contents of container
        """
        del self.name_ij
        del self.name_i
        del self.name_j
        del self.method
        del self.n
        
    def __str__(self):
        """
        'Magic' method for printng contents of container 
        """
        return " %s +  %s x %d ( %s ) -> %s "%(self.name_i,self.name_j,self.n,self.method,self.name_ij)
