# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

"""
This module defines the classes relating to general particles 
"""

__author__ = "Travis W. Kemper, Ph.D."
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"


import logging
logger = logging.getLogger(__name__)

import pymatgen_core.core.periodic_table as periodictable
import pymatgen_core.core.units as units

class Particle(units.ObjectUnits):
    """Data structure for describing any localized object in QM/MD simulation, such as an atom.
    

    Kwargs:
        * type   (str): Particle type
        * label  (str): Identifier in output files
        * symbol (str): Atomic symbol
        * unit_conf (dict): Unit types with units used by this object
            
    Particles have fundamental attributes which are used throughout the streamm code:
    
    .. attribute:: mass (float)

        Mass of the particle 
    
    .. attribute:: charge (float)
    
        Charge of the particle in e 

    .. attribute:: bonded_radius (float)
    
        Radius in used in determining if two particles are bonded 

    .. attribute:: nonbonded_radius (float)

        Radius in used in determining if two particles are interacting
    
    Additional attributes include:
    
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
    

    """
 
    @property
    def mass(self):
        return self._property['mass'] 
    
    @mass.setter
    def mass(self,value):
        self._property['mass']  = value
        
    @property
    def charge(self):
        return self._property['charge']
    
    @charge.setter
    def charge(self,value):
        self._property['charge'] = value


    @property
    def bonded_radius(self):
        return self._property['bonded_radius']
    
    @bonded_radius.setter
    def bonded_radius(self,value):
        self._property['bonded_radius'] = value

    @property
    def nonbonded_radius(self):
        return self._property['nonbonded_radius']
    
    @nonbonded_radius.setter
    def nonbonded_radius(self,value):
        self._property['nonbonded_radius'] = value

    def set_element(self,symbol=None,number=None,mass=None,mass_precision=0):
        '''
        Set the element of property of the particle
        
        Kwargs:
            * symbol (str): Atomic symbol.
            * number (int): Atomic number
            * mass   (float): Atomic mass (AMU)
            * mass_precision (int): precision to match the mass with value from periodic table 
            
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
        self.mass = float(self.element.atomic_mass) # amu 
        if( self.element.atomic_radius_calculated != None ):
           self.bonded_radius = self.element.atomic_radius_calculated
        else:
           self.nonbonded_radius = 1.5
        if( self.element.van_der_waals_radius != None ):
           self.nonbonded_radius = self.element.van_der_waals_radius
        else:
           self.nonbonded_radius = 2.5
            
        logger.debug("Particle has been set to {} with mass:{} bonded_radius:{} nonbonded_radius:{}".format(self.symbol,self.mass,self.bonded_radius,self.nonbonded_radius))
        # Set values to be the same as mendeleev for easy
        # upgrade in next revision 
        self.element.atomic_weight = self.element.atomic_mass
        self.element.covalent_radius = self.element.atomic_radius_calculated
        self.element.vdw_radius = self.element.van_der_waals_radius
        
        return
    

    def __init__(self,type='atom',label=None,symbol = None,unit_conf=units.unit_conf):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        #
        logger.debug("Particle created type:{} label:{} symbol:{}".format(type,label,symbol))
        # 
        self.type = type
        self.label = label
        self.symbol = symbol
        # 
        self._property['mass']             = 1.0  # Choose reasonable values for initialization
        self._property['charge']           = 0.0  # Choose reasonable values for initialization     
        # 
        self._property['bonded_radius']    = 1.0   # Choose reasonable values for initialization 
        self._property['nonbonded_radius'] = 2.0   # Choose reasonable values for initialization
        #
        self._property_units['mass'].append('mass')
        self._property_units['charge'].append('charge')
        self._property_units['length'].append('bonded_radius')
        self._property_units['length'].append('nonbonded_radius')
        # 
        self.mol      = 0
        self.ring     = 0
        self.residue  = 0
        self.resname  = "RES"
        self.qgroup   = 0
        #   
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
            

        # Force field
        self.paramkey = None 
        self.param = None
        # Lammps and gromacs index
        self.param_index = 0 
        self.lammps_index = 0 
        self.gromacs_index = 0             
            
    def __del__(self):
        del self.type
        del self.label
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
        del self.paramkey
        del self.param
        del self.param_index
        del self.lammps_index
        del self.gromacs_index
        #
        del self.rsite
        
    
    def __str__(self):
        return "{} {} ({})".format(self.type,self.label,self.symbol)
    
    def show_attributes(self):
        '''
        Like __str__ but with all the values of the instance's attributes 
        '''
        property_msg = " type:{} ".format(self.type)
        property_msg += "\n label:{}".format(self.label)
        property_msg += "\n symbol:{}".format(self.symbol)
        property_msg += "\n mass:{} ({})".format(self.mass,self._unit_conf['mass'])
        property_msg += "\n charge:{} ({})".format(self.charge,self._unit_conf['charge'])
        property_msg += "\n bonded_radius:{} ({})".format(self.bonded_radius,self._unit_conf['length'])
        property_msg += "\n nonbonded_radius:{} ({})".format(self.nonbonded_radius,self._unit_conf['length'])
        
        return property_msg
    
    def export_json(self):
        '''    
        Export object to json
        
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        
        json_data = {}
        json_data['type'] = self.type
        json_data['label'] = self.label
        json_data['symbol'] = self.symbol
        json_data['mass'] = self.mass
        json_data['charge'] = self.charge
        json_data['bonded_radius'] = self.bonded_radius
        json_data['nonbonded_radius'] = self.nonbonded_radius
        #
        json_data['mol'] = self.mol
        json_data['ring'] = self.ring
        json_data['residue'] = self.residue
        json_data['resname'] = self.resname
        json_data['qgroup'] = self.qgroup
        json_data['paramkey'] = self.paramkey
        json_data['param_index'] = self.param_index
        json_data['lammps_index'] = self.lammps_index
        json_data['gromacs_index'] = self.gromacs_index
        json_data['rsite'] = self.rsite
        
        return json_data
        

    def import_json(self,json_data):
        '''    
        Export object to json
        
        Args:
            * json_data (dict) json representation of the object
            
        '''
        
        self.type  =          json_data['type'] 
        self.label  =          json_data['label'] 
        self.symbol  =          json_data['symbol']
        # Set element 
        self.set_element(symbol = self.symbol )
        
        self.mass  =          json_data['mass'] 
        self.charge  =          json_data['charge'] 
        self.bonded_radius  =          json_data['bonded_radius'] 
        self.nonbonded_radius  =          json_data['nonbonded_radius'] 
        self.mol  =          json_data['mol'] 
        self.ring  =          json_data['ring'] 
        self.residue  =          json_data['residue'] 
        self.resname  =          json_data['resname'] 
        self.qgroup  =          json_data['qgroup'] 
        self.paramkey  =          json_data['paramkey'] 
        self.param_index  =          json_data['param_index'] 
        self.lammps_index  =          json_data['lammps_index'] 
        self.gromacs_index  =          json_data['gromacs_index'] 
        self.rsite  =          json_data['rsite'] 
        
        
        
