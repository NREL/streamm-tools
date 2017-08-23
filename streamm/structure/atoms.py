# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

"""
This module defines the classes relating to atoms
"""

import logging
logger = logging.getLogger(__name__)

try:
    # Import streamm Classes 
    from streamm.structure.particles import Particle
except:
    logger.warning("streamm is not installed test will use relative path")
    from particles import Particle

try:
    # Import pymatgen Class 
    import pymatgen.core.periodic_table as pymatgen_pt
except:
    logger.warning("pymatgen import error for Lattice object")
    exit()
    

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

# Then you have types particles
class Atom(Particle):
    """
    A derived type of particles for atoms in a simulation
    
    """
    def __init__(self,symbol="X",type="atom"):
        '''
        Constructor for Atom object
        
        Args:
            symbol (str) Atomic symbol
            label  (str) Atomic label
            
        '''
        Particle.__init__(self, type=type)
        # 
        # Get properties of element based on symbol
        #
        self.element = pymatgen_pt.Element(symbol)
        self.symbol = self.element.symbol
        # Set values to be the same as mendeleev for easy
        # upgrade in next revision 
        self.element.atomic_weight = self.element.atomic_mass
        self.element.covalent_radius = self.element.atomic_radius_calculated
        self.element.vdw_radius = self.element.van_der_waals_radius

        self.mol = 0
        self.ring = 0
        self.residue = 0
        self.resname = "RES"
        self.qgroup = 0
        

    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.element
        del self.mol
        del self.ring
        del self.residue
        del self.resname
        del self.qgroup
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " {}".format(self.element.symbol)

        
