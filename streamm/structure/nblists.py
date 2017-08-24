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
    Class for neighbor list
    """
    def __init__(self):
        """
        Constructor
        """
        self.list = []
        self.index = []
        self.cnt = -1 

    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        del self.list
        del self.index
        del self.cnt
        
    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " NBlist of {} particle with {} connections".format(len(self.index),len(self.list))
    

    def getnbs(self,key_i):
        '''
        Return list of keys of neighbors of key_i        
        '''
        try:
            nbs_i = self.list[self.index[key_i]:self.index[key_i+1]]
        except:
            logger.warning(" Neighbor list not set ")
            nbs_i = []
            
        return nbs_i

    def calc_nnab(self,key_i):
        '''
        Calculate the number of neighbors for particle key_i
        '''
        try:
            nnab_i = self.index[key_i+1] - self.index[key_i]
        except:
            logger.warning(" Neighbor list not set ")
            nnab_i = 0
        #
        return nnab_i

