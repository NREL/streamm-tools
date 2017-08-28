# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"




class Dihtype(object):
    """
    Set of Dihedral angle parameters 
    """

    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank", fftype4="blank" , type="multiharmonic" ):
        """
        Constructor for a angle parameter.
        
        Args:
             fftype1  (str)   Atom type 
             fftype2  (str)   Atom type 
             fftype3  (str)   Atom type 
             fftype4  (str)   Atom type 
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

        if isinstance(fftype3, str):
            self.fftype3 = fftype3
        else:
            print "3rd arg should be str"
            raise TypeError

        if isinstance(fftype4, str):
            self.fftype4 = fftype4
        else:
            print "4th arg should be str"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "5th arg should be str"
            raise TypeError

        # Set default values for parameters
        self.paths = 1 
        self.d = 1.0    # gamma 
        self.mult = 0.0 # n 
        self.kb = 0.0
        self.theat_s = 0.0

        # Fourier opls coefficients 
        self.k1 = 0.0 
        self.k2 = 0.0 
        self.k3 = 0.0 
        self.k4 = 0.0 

        # Ryckaert-Bellemans function coefficients 
        self.C0 = 0.0 
        self.C1 = 0.0 
        self.C2 = 0.0 
        self.C3 = 0.0 
        self.C4 = 0.0 
        self.C5 = 0.0


        self.e0 = 0.0
        self.ke = 1.0 

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
        del self.d
        del self.mult
        del self.kb
        del self.theat_s
        del self.k1
        del self.k2
        del self.k3
        del self.k4
        del self.C0
        del self.C1
        del self.C2
        del self.C3
        del self.C4
        del self.C5
        del self.e0
        del self.ke
        del self.lammps_index
        del self.gromacs_index 

    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " dihedral  %s - %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.fftype4,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic d = %f mult = %f K = %f theat_s = %f lammps index %d  gromcas index %d " %(self.d,self.mult ,self.kb,self.theat_s,self.lammps_index ,self.gromacs_index )
        if( self.type ==  "multiharmonic" ):
            strucStr += "\n  harmonic d = %f mult = %f K = %f theat_s = %f lammps index %d  gromcas index %d " %(self.d,self.mult ,self.kb,self.theat_s,self.lammps_index ,self.gromacs_index )
        if( self.type ==  "opls" ):
            strucStr += "\n  k1 = %f k2 = %f k3 = %f k4 = %f lammps index %d  gromcas index %d " %(self.k1,self.k2,self.k3,self.k4,self.lammps_index ,self.gromacs_index )
        if( self.type ==  "rb" ):
            strucStr += "\n  C0 = %f  C1 = %f C2 = %f C3 = %f C4 = %f  C5 = %f lammps index %d  gromcas index %d " %(self.C0,self.C1,self.C2,self.C3,self.C4,self.C5 ,self.lammps_index ,self.gromacs_index)

        return strucStr


    def setharmonic(self, mult, kb,theat_s):
        """
        set MultiHarmonic parameters
        dihedral_style charmm

        E = kb[ 1 - cos( mult theta - theat_s ) ]  gromacs
        E = kb[ 1 - cos( n theta - d ) ]           lammps 

        Args:
            mult     (float) 
            kb     (float) force constant    kcal/mol
            theat_s     (float) angle degrees 
        """

        if isinstance(kb, float):
            self.kb = kb
        else:
            print "1nd arg kb should be float"
            raise TypeError

        if isinstance(mult, float):
            self.mult = mult
        else:
            print "2rd arg mult should be float"
            raise TypeError

        if isinstance(theat_s, float):
            self.theat_s = theat_s
        else:
            print "3th arg theat_s should be float"
            raise TypeError

    def setopls(self,k1,k2,k3,k4):
        """
        set opls parameters

        E = 1/2 k1[1+cos(theta)]+1/2 k2[1-cos(2 theta)]+1/2 k3[1+cos(3 theta)]+1/2 k4[1-cos(4 theta)]

        Args:
            k1     (float) force constant    kcal/mol
            k2     (float) force constant    kcal/mol
            k3     (float) force constant    kcal/mol
            k4     (float) force constant    kcal/mol
        """

        if isinstance(k1, float):
            self.k1 = k1
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(k2, float):
            self.k2 = k2
        else:
            print "2nd arg should be float"
            raise TypeError

        if isinstance(k3, float):
            self.k3 = k3
        else:
            print "3rd arg should be float"
            raise TypeError

        if isinstance(k4, float):
            self.k4 = k4
        else:
            print "4th arg should be float"
            raise TypeError

        # Translate to Ryckaert-Bellemans function
        self.C0 = k2 + 0.5*(k1+k3)
        self.C1 = 0.5*(-1.0*k1+3.0*k3)
        self.C2 = -1.0*k2 + 4.0*k3
        self.C3 = -2.0*k3
        self.C4 = -4.0*k4
        self.C5 = 0.0


    def setrb(self,C0,C1,C2,C3,C4,C5):
        """
        set Ryckaert-Bellemans parameters

        V_{rb}(theta) = \sum_n=0^5 C_n [ cos(theata - 180 )^n ]

        Args:
            C0     (float) force constant    kcal/mol
            C1     (float) force constant    kcal/mol
            C2     (float) force constant    kcal/mol
            C3     (float) force constant    kcal/mol
            C4     (float) force constant    kcal/mol
            C5     (float) force constant    kcal/mol
        """

        if isinstance(C0, float):
            self.C0 = C0
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(C1, float):
            self.C1 = C1
        else:
            print "2nd arg should be float"
            raise TypeError

        if isinstance(C2, float):
            self.C2 = C2
        else:
            print "3rd arg should be float"
            raise TypeError

        if isinstance(C3, float):
            self.C3 = C3
        else:
            print "4th arg should be float"
            raise TypeError

        if isinstance(C4, float):
            self.C4 = C4
        else:
            print "5th arg should be float"
            raise TypeError

        if isinstance(C5, float):
            self.C5 = C5
        else:
            print "6th arg should be float"
            raise TypeError
        # Translate to opls 
        self.k1 = -1.0*( 2.0*C1 + 3.0*C3/2.0)
        self.k2 = -1.0*( C2 + C4)
        self.k3 = -0.5*C3
        self.k4 = -0.25*C4

        