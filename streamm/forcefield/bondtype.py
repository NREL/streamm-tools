# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"





class Bondtype(object):
    """
    Set of bond parameters 
    """

    def __init__(self, fftype1='blank', fftype2='blank', type="harmonic" ):
        """
        Constructor for a bond parameter.
        
        Args:
             fftype1  (str)   Atom type 
             fftype2  (str)   Atom type 
             type    (str)   Bond type 
        """

        if isinstance(fftype1, str):
            self.fftype1 = fftype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(fftype2, str):
            self.fftype2 = fftype2
        else:
            print "2nd arg should be str"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "3rd arg should be str"
            raise TypeError

        # Set default values for parameters
        self.r0 = 0.0
        self.kb = 0.0 

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
            strucStr += "\n  harmonic r_0 = %f K = %f lammps index %d  gromcas index %d  " %(self.r0 ,self.kb,self.lammps_index ,self.gromacs_index )
                
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
            print "1st arg should be float"
            raise TypeError

        if isinstance(kb, float):
            self.kb = kb
        else:
            print "2nd arg should be float"
            raise TypeError
    