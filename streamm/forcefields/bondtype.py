# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"


try:
    # Import pymatgen Class 
    import pymatgen_core.core.units as units 
except:
    raise ImportError("pymatgen import error for units object")


class Bondtype(object):
    """
    Set of bond parameters

    Args:
         fftype1  (str):   Atom type 
         fftype2  (str):   Atom type 
         type    (str):   Bond type
         
    """

    @property
    def unit_conf(self):
        return self._unit_conf
    
    @property
    def r0(self):
        return self._property['r0'] 
    
    @r0.setter
    def r0(self,value):
        self._property['r0']  = value
        
    @property
    def kb(self):
        return self._property['kb'] 
    
    @kb.setter
    def kb(self,value):
        self._property['kb']  = value
        
    def __init__(self, fftype1='blank', fftype2='blank', type="harmonic",unit_conf=units.unit_conf):
        # Store the units of each attribute type 
        self._unit_conf = unit_conf
        #
        # Default Physical properties
        #
        self._property = {}
        self._property_units = {}
        for unit_type in self._unit_conf.keys():
            self._property_units[unit_type] = []        

        self.fftype1 = fftype1
        self.fftype2 = fftype2
        self.type = type

        # Set default values for parameters
        self._property['kb'] = 1.0
        self._property['r0'] = 1.0
        # 
        self._property_units['harm_bond_coeff'].append('kb')
        self._property_units['length'].append('r0')
        
        # Lammps and gromacs index
        self.lammps_index = 0 
        self.gromacs_index = 0 


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.fftype1 
        del self.fftype2
        del self.type 
        del self.r0
        del self.kb 
        del self.lammps_index 
        del self.gromacs_index 
        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " bond  %s - %s type %s "%(self.fftype1,self.fftype2,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic r_0 = %f K = %f lammps index %d  gromacs index %d  " %(self.r0 ,self.kb,self.lammps_index ,self.gromacs_index )
                
        return strucStr

    def setharmonic(self, r0, kb):
        """
        set Harmonic parameters

        E = kb( r - r0 )^2 

        Args:
            r0 (float) distance 
            kb (float) force constant
            
        """

        
        if isinstance(r0, float):
            self.r0 = r0
        else:
            raise TypeError("1st arg should be float")

        if isinstance(kb, float):
            self.kb = kb
        else:
            raise TypeError("2nd arg should be float")
    

    def update_units(self,new_unit_conf):
        '''
        Update instance values with new units
        
        Args:
            new_unit_conf (dict): with unit type as the key and the new unit as the value
            
        '''
        
        self._property,self._unit_conf = units.change_properties_units(self._unit_conf,new_unit_conf,self._property_units,self._property)
            