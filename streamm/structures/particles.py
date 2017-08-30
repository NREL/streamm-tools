# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

"""
This module defines the classes relating to general particles 
"""

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

import logging
logger = logging.getLogger(__name__)

try:
    # Import pymatgen Class 
    import pymatgen.core.periodic_table as periodictable
except:
    logger.warning("pymatgen import error for Lattice object")
    exit()
    
ATOMIC_MASS_PRECISION = 0

'''    
# Notes to add mendeleev dependency 
try:
    import mendeleev.element
    USE_MENDELEEV = True
except:
    USE_MENDELEEV = False

        if( USE_MENDELEEV ):
            self.element = mendeleev.element(symbol)
            self.symbol = self.element.symbol
            self.atomic_number = self.element.atomic_number
            self.atomic_weight = self.element.atomic_weight
            self.covalent_radius = self.element.covalent_radius
            self.vdw_radius = self.element.vdw_radius
'''        

class Particle(object):
    """
    
    Data structure for describing any localized object in QM/MD simulations
    A 'Particle' has a type and dict of specifiers to it's property
    
    """

    def set_element(self,symbol=None,number=None,mass=None):
        '''
        Set the element of particle 
        '''
        if( symbol != None ):
            self.symbol = str(symbol)
            logger.info("Finding element by atomic symbol {}".format(self.symbol))
            self.element = periodictable.Element(self.symbol)
            
        elif( number != None ):
            number = int(number)
            logger.info("Finding element by atomic number {}".format(number))
            self.element = periodictable.Element.from_Z(number)
            self.symbol = self.element.symbol
            
        elif( mass != None ):
            mass = float(mass)
            logger.info("Finding element by atomic mass {} ".format(mass))
            for symbol, data in periodictable._pt_data.items():
                if round(data["Atomic mass"],ATOMIC_MASS_PRECISION) == round(mass,ATOMIC_MASS_PRECISION):
                    self.symbol = symbol
                    self.element = periodictable.Element(symbol)
                    break
                
            if( self.symbol  == None ):
                logger.warning("Atomic mass of {} was not found in periodic table with precision of {}".format(mass,ATOMIC_MASS_PRECISION))
                return 
        else:
            logger.warning("No arguments supplied to function element will not be set")
            return
        #
        # Set properties based on element properties 
        #
        self.mass = self.element.atomic_mass
        self.bonded_radius = self.element.atomic_radius_calculated
        self.nonbonded_radius = self.element.van_der_waals_radius

        # Set values to be the same as mendeleev for easy
        # upgrade in next revision 
        self.element.atomic_weight = self.element.atomic_mass
        self.element.covalent_radius = self.element.atomic_radius_calculated
        self.element.vdw_radius = self.element.van_der_waals_radius
        
        return
    
    def __init__(self,type='atom',label=None,symbol = None ):
        """
        Constructor for a general particle. 
        """
        #
        self.type = type 
        #
        # Default Physical properties
        #
        self.mass             = 1.0  # Choose reasonable values for initialization 
        self.charge           = 0.0  # Choose reasonable values for initialization     
        #
        self.bonded_radius    = 1.0  # Choose reasonable values for initialization 
        self.nonbonded_radius = 2.0  # Choose reasonable values for initialization 
        # 
        self.mol      = 0
        self.ring     = 0
        self.residue  = 0
        self.resname  = "RES"
        self.qgroup   = 0
        self.index    = None
        #   
        # Force field 
        self.ff = None
        #
        # Atomic properties 
        self.symbol = symbol
        self.element = None
        if( symbol != None ):
            self.set_element(symbol=symbol)
        # If no label is given use symbol 
        if( label == None and symbol != None ):
            self.label  = symbol
        else:     
            self.label  = label
           
        
    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.type
        del self.label
        #
        del self.mass
        del self.charge
        #
        del self.bonded_radius
        del self.nonbonded_radius 
        # 
        del self.mol
        del self.ring
        del self.residue
        del self.resname
        del self.qgroup
        # Atomic properties 
        del self.symbol 
        del self.element 
        # Force field 
        del self.ff 
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return "{}:{}".format(self.type,self.label)
    
