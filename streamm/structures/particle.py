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
    raise ImportError("pymatgen import error for pymatgen.core.periodic_table module")


# Import streamm dependencies 
import streamm.util.units as units
# import streamm.structures.lattice as lattice
# import streamm.structures.lattice as lattice

DIMENSIONS = 3

class Particle(object):
    """Data structure for describing any localized object in QM/MD simulation, such as an atom.
    

    Kwargs:
        * type   (str): Particle type
        * label  (str): Identifier in output files
        * symbol (str): Atomic symbol
            
    Particles have the fundmental attributes of
    
    .. attribute:: mass (float)

        Mass of the particle 
    
    .. attribute:: charge (float)
    
        Charge of the particle in e 

    .. attribute:: bonded_radius (float)
    
        Radius in (angstroms) used in determining if two particles are bonded 

    .. attribute:: nonbonded_radius (float)

        Radius in (angstroms) used in determining if two particles are interacting
    
    which are used through out the streamm code.  Additional attributes include:
    
    .. attribute:: mol (int)
        
        Molecule index
    
    .. attribute:: ring (int)
        
        Ring index
    
    .. attribute:: residue (int)

        Residue index 
    
    .. attribute:: resname (str)

        Residue name 
    
    .. attribute:: qgroup (int)

        Charge group index 
    
    .. attribute:: index (int)

        Index of particle in particle dictionary 
    
    .. attribute:: ffkey (str)

        Particletype key
        
    
    The ``element`` attribute can be set to an element object from the
    ``periodic_table`` object in ``pymatgen``.
    
    The ``ff`` attribute can be set to a ``Particletype`` object. 
    
    .. TODO:
        * create function p_i.get_properties()
        * set position in particle object 

    """
    

    def set_element(self,symbol=None,number=None,mass=None,mass_precision=0):
        '''
        Set the element of property of the particle
        
        Kwargs:
            * symbol (str): Atomic symbol.
            * number (int): Atomic number
            * mass   (float): Atomic mass (AMU)
            * mass_precision (int): precision to match the mass with value from
            periodic table 
            
        This will set the ``symbol`` property of the particle to the Atomic symbol,
        the ``mass`` property of the particle to the ``atomic_mass``,
        the ``bonded_radius`` property of the particle to the ``atomic_radius_calculated`` and 
        the ``nonbonded_radius`` property of the particle to the ``van_der_waals_radius``
        
        '''
        if( symbol != None ):
            self.symbol = str(symbol)
            logger.debug("Finding element by atomic symbol {}".format(self.symbol))
            self.element = periodictable.Element(self.symbol)
            
        elif( number != None ):
            number = int(number)
            logger.debug("Finding element by atomic number {}".format(number))
            self.element = periodictable.Element.from_Z(number)
            self.symbol = self.element.symbol
            
        elif( mass != None ):
            mass = float(mass)
            logger.debug("Finding element by atomic mass {} ".format(mass))
            for symbol, data in periodictable._pt_data.items():
                if round(data["Atomic mass"],mass_precision) == round(mass,mass_precision):
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
        if( self.element.atomic_radius_calculated != None ):
            self.bonded_radius = units.Length(self.element.atomic_radius_calculated,'ang')                        
        if( self.element.van_der_waals_radius != None ):
            self.nonbonded_radius = units.Length(self.element.van_der_waals_radius,'ang')
        logger.info("Particle[{}] has been set to C with mass:{} bonded_radius:{} nonbonded_radius:{}".format(self.index,self.mass,self.bonded_radius,self.nonbonded_radius))
        # Set values to be the same as mendeleev for easy
        # upgrade in next revision 
        self.element.atomic_weight = self.element.atomic_mass
        self.element.covalent_radius = self.element.atomic_radius_calculated
        self.element.vdw_radius = self.element.van_der_waals_radius
        
        return
    
    def __init__(self,type='atom',label=None,symbol = None,unit_conf=units.unit_conf ):
        #
        logger.debug("Particle created type:{} label:{} symbol:{}".format(type,label,symbol))
        self.type = type
        self.label = label
        self.symbol = symbol
        # Store the units of each attribute type 
        self.unit_conf = unit_conf  
        #
        # Default Physical properties
        #
        self.mass             = 1.0  # Choose reasonable values for initialization 
        self.charge           = 0.0  # Choose reasonable values for initialization     
        # 
        self.bonded_radius    = 1.0   # Choose reasonable values for initialization 
        self.nonbonded_radius = 2.0   # Choose reasonable values for initialization 
        # 
        self.mol      = 0
        self.ring     = 0
        self.residue  = 0
        self.resname  = "RES"
        self.qgroup   = 0
        self.index    = None
        #   
        # Force field
        self.ffkey = None 
        self.ff = None
        # Reactive site type
        self.rsite = '' 
        #
        # Atomic properties 
        self.symbol = symbol
        self.element = None
        #
        if( symbol != None ):
            self.set_element(symbol=symbol)
            logger.debug("No label is given using symbol as label")
            if( label == None  ):
                self.label  = symbol
            
        elif( label != None ):
            logger.debug("No symbol is given using label as symbol")
            self.symbol = label
        
    def __del__(self):
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
        #
        del self.rsite
        
    def __str__(self):
        return "{}[{}] {} ({})".format(self.type,self.index,self.label,self.symbol,)
    
