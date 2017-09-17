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

class Imptype(object):
    """
    Set of improper dihedral angle parameters

    Args:
         fftype1  (str)   Atom type 
         fftype2  (str)   Atom type 
         fftype3  (str)   Atom type 
         fftype4  (str)   Atom type 
         type    (str)   Bond type
         
    """

    @property
    def unit_conf(self):
        return self._unit_conf
        
    @property 
    def e0(self): 
        return self._property['e0']
    @e0.setter 
    def e0(self,value): 
        self._property['e0'] = value
    @property 
    def ke(self): 
        return self._property['ke']
    @ke.setter 
    def ke(self,value): 
        self._property['ke'] = value
    @property 
    def pn(self): 
        return self._property['pn']
    @pn.setter 
    def pn(self,value): 
        self._property['pn'] = value
    
    
    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank", fftype4="blank" , type="improper" ,unit_conf=units.unit_conf):

        # Store the units of each attribute type 
        self._unit_conf = unit_conf
        #
        # Default Physical properties
        #
        self._property = {}
        self._property_units = {}
        for unit_type in self._unit_conf.keys():
            self._property_units[unit_type] = []
        # 
        self.fftype1 = fftype1
        self.fftype2 = fftype2
        self.fftype3 = fftype3
        self.fftype4 = fftype4
        self.type = type

        # Set default values for parameters

        self._property['e0']  = 1
        self._property['ke']  = 1.0
        self._property['pn']  = 0.0
        
        self._property_units['energy'].append('ke')
        self._property_units['angle'].append('e0')
        #         
        # Lammps and gromacs index
        self.lammps_index = 0 
        self.gromacs_index = 0 

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.fftype1
        del self.fftype2 
        del self.fftype3
        del self.fftype4
        del self.type
        del self._property
        del self._property_units
        del self.lammps_index
        del self.gromacs_index 

        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " improper  %s - %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.fftype4,self.type)
        
        if( self.type ==  "improper" ):
            strucStr += "\n  imp e0 = %f ke = %f lammps index %d  gromcas index %d " %(self.e0,self.ke,self.lammps_index ,self.gromacs_index )

        return strucStr        


    def setimp(self, e0, ke):
        """
        set Harmonic parameters

        E = kb ( e_{lijk} - e0 )^2 

        Args:
            e0     (float) 
            kb     (float) force constant    kcal/mol
        """

        if isinstance(e0, float):
            self.e0 = e0
        else:
            raise TypeError("1st arg should be float")

        if isinstance(ke, float):
            self.ke = ke
        else:
            raise TypeError("2nd arg should be float")

