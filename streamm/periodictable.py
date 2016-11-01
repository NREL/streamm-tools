#! /usr/bin/env python
"""
Class of elements with associated atomic properties 
"""


__author__ = "Travis W. Kemper"
__version__ = "0.0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Alpha"

import os,json 

import logging
logger = logging.getLogger(__name__)


# Loads element data from json file
with open(os.path.join(os.path.dirname(__file__), "periodic_table.json"), "rt"
          ) as f:
    elements = json.load(f)

verbose = False  
n_elements = len(elements)
# Create an empty element for unknow values
def unknow_el():
      el_i = dict()
      el_i["symbol"] = str("X")
      el_i["number"] = int(-1)
      el_i["mass"] = float(0.0)
      el_i["cov_radii"] = float(0.0)
      el_i["vdw_radii"] = float(0.0)
      return el_i

def element_mass(mass_i):
    """
    Find element based on atomic mass
    """
    mass_i_int = int(mass_i)
    for el_symb in elements.keys():
        el = elements[el_symb]
        el_mass_int = int( el['mass'] )
        if( mass_i_int == el_mass_int ):
            return el
    logger.warning("No element found for atomi with mass %s "%(mass_i))
    el_empty = unknow_el()
    el_empty['mass']  = mass_i
    return el_empty

def element_number(atomic_number_i):
    """
    Find element based on atomic number 
    """
    for el_symb in elements.keys():
        el = elements[el_symb]
        el_n_i = int(el['number'] )
        if( el_n_i == atomic_number_i ):
            return el
        
    logger.warning("No element found for atomic with number %s "%(atomic_number_i))
    el_empty = unknow_el()
    el_empty['number']  = atomic_number_i
    return el_empty

def element_symbol(symbol_i):
    """
    Find element based on atomic symbol  
    """

    for el_symb in elements.keys():
        el = elements[el_symb]
        if( symbol_i == el['symbol'] ):
            return el

    logger.warning("No element found for atomi symbol %s "%(symbol_i))
    el_empty = unknow_el()
    el_empty['symbol']  = symbol_i
    return el_empty



