#! /usr/bin/env python
"""
This module defines the classes relating to classical particles in box

The default units are

distance - Angstroms 
mass - AMU

Units are changed when a structure is associated with a simulation 

"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

# Standard packages 
from numpy import pi #, dot, transpose, radians
import copy
import math
import random
import sys
import json
import csv
import os

try:
    import cPickle as pickle
except:
    import pickle
    
from datetime import datetime

# Dependency packages 
import numpy as np
import pandas as pd

# pymatgen module
import pymatgen.core.periodic_table as periodictable

import logging
logger = logging.getLogger(__name__)



'''
Extra functions  
'''
                      
def prop_list(propkey,keys,dic):
    '''
    Calculate bond lengths all bond according to list of keys
    '''
    if( len(keys) == 0 ):
        # If list is not specified use all bonds
        keys = dic.keys()
    prop_list = []
    for key_i in keys:
        dic_entry = dic[key_i]
        prop_list.append(dic_entry.properties[propkey])

    return prop_list


def hterm_Csp3(hb_length,r_i,r_ij_array):
    """
    Hydrogen terminate segment 

     (j)     (l)
        \   /
         (i)   
        /   \
     (k)     (m)



    """

    add_jk = np.zeros(3)
    for r_ij in r_ij_array:
        add_jk +=  r_ij
    add_jk  = -1.0*add_jk
    add_scale = hb_length*add_jk/np.linalg.norm(add_jk)

    r_l = r_i + add_scale

    return r_l 

def hterm_Csp2(hb_length,r_i,r_ij_array):
    """
    Hydrogen terminate conjugated atom

    (j)  
        \  
         (i)  - (l)
        /  
    (k)
    
    Hydrogens will be added to sp2 carbons in a plainer configuration

    """
    debug = False

    add_jk = -1.0*( r_ij_array[0] + r_ij_array[1] )
    add_scale = hb_length*add_jk/np.linalg.norm(add_jk)
    r_l = r_i + add_scale

    return r_l


