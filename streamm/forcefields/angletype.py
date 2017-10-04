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

   
class Angletype(units.ObjectUnits):
    """
    Set of Angle parameters


    Kwargs:
        * fftype1  (str):   Atom type 
        * fftype2  (str):   Atom type 
        * fftype3  (str):   Atom type 
        * type    (str):  Bond type
        * units_conf (dict): Dictionary of units for each attribute type

    """

    @property
    def kb(self):
        return self._property['kb'] 
    
    @kb.setter
    def kb(self,value):
        self._property['kb']  = value
        
    @property
    def theta0(self):
        return self._property['theta0'] 
    
    @theta0.setter
    def theta0(self,value):
        self._property['theta0']  = value
        
    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank" , type="harmonic" ,unit_conf=units.unit_conf):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        #
        
        self.fftype1 = fftype1
        self.fftype2 = fftype2
        self.fftype3 = fftype3
        self.type = type

        # Set default values for parameters
        self._property['kb']  = 0.0 
        self._property['theta0']  = 0.0

        self._property_units['energy'].append('kb')
        self._property_units['angle'].append('theta0')
        
        # Lammps and gromacs index
        self.lammps_index = 0 
        self.gromacs_index = 0 


    def __del__(self):
        
        del self.fftype1
        del self.fftype2 
        del self.fftype3
        del self.type
        del self.lammps_index
        del self.gromacs_index 


    def __str__(self):
        strucStr =  " angle  %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic theta_0 = %f K = %f lammps index %d  gromacs index %d  " %(self.theta0 ,self.kb,self.lammps_index ,self.gromacs_index )

        return strucStr


    def setharmonic(self, theta0, kb):
        """
        set Harmonic angle parameters

        Args:
            * theta0 (float): angle           
            * kb     (float): force constant  
            
        .. math ::
            E = kb( theta - theta_0 )^2 

        """

        if isinstance(theta0, float):
            self.theta0 = theta0
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(kb, float):
            self.kb = kb
        else:
            print "2nd arg should be float"
            raise TypeError

