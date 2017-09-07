# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

'''
This module defines the classes relating to neighbor lists 
'''

import logging
logger = logging.getLogger(__name__)
              
class NBlist(object):
    """
    Class for neighbor list.
    
    A neighbor list has two Instance variable:
    * list (list): 
    * index (list):
    
    Where the ``list`` is a list of all neighboring indexes in sequence.
    The ``index`` is list of the beginning of the neighbors sequence based on index.
    
    Example:
    You have particles with indexes 0, 1, 2, 3 and 4. And the following bonds:
    
    * 0 - 1
    * 1 - 2
    * 2 - 3
    * 2 - 4
    
    The ``list`` would be ``[1,0,2,1,3,4,2,2]`` and the ``index`` would be ``[0,1,3,6,7,8]``.
    So if we want the neighbor indexes of particle ``2``,
    we get the beginning position in ``list`` from ``index[2]`` to be 3 and the ending positin in ``list``
    for  ``index[2+1]-1`` to be 5. This gives us the neighbor list of particle ``1`` to be ``[1,3,4]``.
    """
    def __init__(self):
        self.list = []
        self.index = []
        self.cnt = -1 

    def __del__(self):
        del self.list
        del self.index
        del self.cnt
        
    def __str__(self):
        return " NBlist of {} particle with {} connections".format(len(self.index),len(self.list))
    

    def getnbs(self,key_i):
        '''
        Return list of in of neighbors of key_i
        
        Args:
            key_i (int): Index of item in neighbor list
        Returns:
            nb_list (list): list of indexes of neighbors of key_i
        '''
        try:
            return self.list[self.index[key_i]:self.index[key_i+1]]
        except:
            raise KeyError(" Neighbor list not set ")
            return []
            

    def calc_nnab(self,key_i):
        '''
        Calculate the number of neighbors for particle key_i
        
        Args:
            key_i (int): Index of item in neighbor list
        Returns:
            (int): number of neighbors of key_i        
        '''
        try:
            nnab_i = self.index[key_i+1] - self.index[key_i]
        except:
            raise KeyError(" Neighbor list not set ")
            return []
            

