# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"




   
class Angletype(object):
    """
    Set of Angle parameters 
    """

    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank" , type="harmonic" ):
        """
        Constructor for a angle parameter.
        
        Args:
             fftype1  (str)   Atom type 
             fftype2  (str)   Atom type 
             fftype3  (str)   Atom type 
             type    (str)   Bond type 
        """

        self.fftype1 = fftype1
        self.fftype2 = fftype2
        self.fftype3 = fftype3
        self.type = type

        # Set default values for parameters
        self.theta0 = 0.0
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
        del self.fftype3
        del self.type
        del self.theta0
        del self.kb 
        del self.lammps_index
        del self.gromacs_index 


    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " angle  %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic theta_0 = %f K = %f lammps index %d  gromcas index %d  " %(self.theta0 ,self.kb,self.lammps_index ,self.gromacs_index )

        return strucStr


    def setharmonic(self, theta0, kb):
        """
        set Harmonic angle parameters

        E = kb( theta - theta_0 )^2 

        Args:
            theta0 (float) angle             deg  
            kb     (float) force constant    kcal/mol
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

