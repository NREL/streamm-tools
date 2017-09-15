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

    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank", fftype4="blank" , type="improper" ,unit_conf=units.unit_conf):
        self.unit_conf = unit_conf
     
        self.fftype1 = fftype1
        self.fftype2 = fftype2
        self.fftype3 = fftype3
        self.fftype4 = fftype4
        self.type = type

        # Set default values for parameters
        self.e0 = 0.0
        self.ke = 1.0
        self.pn = 0.0    # For periodicimproper

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
        del self.ke
        del self.e0
        del self.pn
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
            print "1st arg should be float"
            raise TypeError

        if isinstance(ke, float):
            self.ke = ke
        else:
            print "2nd arg should be float"
            raise TypeError

