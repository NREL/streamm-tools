"""
Project is a set of calcultions

Projects track:
  file IO
  simulations ran
  computational resources used for each simulation

Data for each action is stored in a json file

"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

import copy , os , json , sys
import time, datetime
from string import replace

import logging
logger = logging.getLogger(__name__)

# streamm
from resource import Resource 
import resource,buildingblock, calculation, structure


class Project:
    '''
    Data structure for a project 
    '''
    def __init__(self,tag):
        """
        Constructor for a general project object.
        
        Args:
            tag (str): String identifier for object
        """
        
        if isinstance(tag, str):
            self.tag = tag   # String name of simulation code (eg GaussianJuly21)
        else:
            raise TypeError("1st arg (tag) in %s Project initialization should be string"%(__name__))
        
        self.calculations = dict()
        
    def __del__(self):
        """
        Delete Calculation object
        """
        del self.tag
        del self.calculations
