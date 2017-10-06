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


class Dihtype(units.ObjectUnits):
    """
    Set of Dihedral angle parameters

    Args:
         * fftype1  (str)   Atom type 
         * fftype2  (str)   Atom type 
         * fftype3  (str)   Atom type 
         * fftype4  (str)   Atom type 
         * type     (str)   Bond type
         
    """

    @property 
    def paths(self): 
        return self._property['paths']
    @paths.setter 
    def paths(self,value): 
        self._property['paths'] = value
    @property 
    def kb(self): 
        return self._property['kb']
    @kb.setter 
    def kb(self,value): 
        self._property['kb'] = value
    @property 
    def d(self): 
        return self._property['d']
    @d.setter 
    def d(self,value): 
        self._property['d'] = value
    @property 
    def mult(self): 
        return self._property['mult']
    @mult.setter 
    def mult(self,value): 
        self._property['mult'] = value
    @property 
    def theta_s(self): 
        return self._property['theta_s']
    @theta_s.setter 
    def theta_s(self,value): 
        self._property['theta_s'] = value
    @property 
    def k1(self): 
        return self._property['k1']
    @k1.setter 
    def k1(self,value): 
        self._property['k1'] = value
    @property 
    def k2(self): 
        return self._property['k2']
    @k2.setter 
    def k2(self,value): 
        self._property['k2'] = value
    @property 
    def k3(self): 
        return self._property['k3']
    @k3.setter 
    def k3(self,value): 
        self._property['k3'] = value
    @property 
    def k4(self): 
        return self._property['k4']
    @k4.setter 
    def k4(self,value): 
        self._property['k4'] = value
    @property 
    def C0(self): 
        return self._property['C0']
    @C0.setter 
    def C0(self,value): 
        self._property['C0'] = value
    @property 
    def C1(self): 
        return self._property['C1']
    @C1.setter 
    def C1(self,value): 
        self._property['C1'] = value
    @property 
    def C2(self): 
        return self._property['C2']
    @C2.setter 
    def C2(self,value): 
        self._property['C2'] = value
    @property 
    def C3(self): 
        return self._property['C3']
    @C3.setter 
    def C3(self,value): 
        self._property['C3'] = value
    @property 
    def C4(self): 
        return self._property['C4']
    @C4.setter 
    def C4(self,value): 
        self._property['C4'] = value
    @property 
    def C5(self): 
        return self._property['C5']
    @C5.setter 
    def C5(self,value): 
        self._property['C5'] = value
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
        
    

    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank", fftype4="blank" , type="multiharmonic" ,unit_conf=units.unit_conf):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)

        #        
        self.fftype1 = fftype1
        self.fftype2 = fftype2
        self.fftype3 = fftype3
        self.fftype4 = fftype4
        self.type = type
        #
        # Set default values for parameters
        self._property['paths']  = 1
        self._property['d']  = 1.0
        self._property['mult']  = 0.0 
        self._property['kb']  = 0.0
        self._property['theta_s']  = 0.0 
        #
        self._property_units['energy'].append('kb')
        self._property_units['angle'].append('theta_s')
        #         
        # Fourier opls coefficients 
        self._property['k1']  = 0.0
        self._property['k2']  = 0.0
        self._property['k3']  = 0.0
        self._property['k4']  = 0.0

        self._property_units['energy'].append('k1')
        self._property_units['energy'].append('k2')
        self._property_units['energy'].append('k3')
        self._property_units['energy'].append('k4')

        # Ryckaert-Bellemans function coefficients 
        self._property['C0']  = 0.0
        self._property['C1']  = 0.0
        self._property['C2']  = 0.0
        self._property['C3']  = 0.0
        self._property['C4']  = 0.0
        self._property['C5']  = 0.0
        
        self._property_units['energy'].append('C0')
        self._property_units['energy'].append('C1')
        self._property_units['energy'].append('C2')
        self._property_units['energy'].append('C3')
        self._property_units['energy'].append('C4')
        self._property_units['energy'].append('C5')
        
        self._property['e0']  = 0.0
        self._property['ke']  = 0.0
        #
        self._property_units['angle'].append('e0')
        self._property_units['energy'].append('ke')

        # Lammps and gromacs index
        self.lammps_index = 0 
        self.gromacs_index = 0

    def __del__(self):
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
        strucStr =  " dihedral  %s - %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.fftype4,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic d = %f mult = %f K = %f theta_s = %f lammps index %d  gromcas index %d " %(self.d,self.mult ,self.kb,self.theta_s,self.lammps_index ,self.gromacs_index )
        if( self.type ==  "multiharmonic" ):
            strucStr += "\n  harmonic d = %f mult = %f K = %f theta_s = %f lammps index %d  gromcas index %d " %(self.d,self.mult ,self.kb,self.theta_s,self.lammps_index ,self.gromacs_index )
        if( self.type ==  "opls" ):
            strucStr += "\n  k1 = %f k2 = %f k3 = %f k4 = %f lammps index %d  gromcas index %d " %(self.k1,self.k2,self.k3,self.k4,self.lammps_index ,self.gromacs_index )
        if( self.type ==  "rb" ):
            strucStr += "\n  C0 = %f  C1 = %f C2 = %f C3 = %f C4 = %f  C5 = %f lammps index %d  gromcas index %d " %(self.C0,self.C1,self.C2,self.C3,self.C4,self.C5 ,self.lammps_index ,self.gromacs_index)

        return strucStr


    def setharmonic(self, mult, kb,theta_s):
        """
        set MultiHarmonic parameters
        dihedral_style charmm

        Args:
            * mult     (float) 
            * kb     (float) force constant    kcal/mol
            * theta_s     (float) angle degrees
            
        gromacs
        
        .. math::
            E = kb[ 1 - cos( mult*theta - theta_s ) ]  

        lammps
        
        .. math::
            E = kb[ 1 - cos( n*theta - d ) ]            
            
        """

        if isinstance(kb, float):
            self.kb = kb
        else:
            raise TypeError("1nd arg kb should be float")

        if isinstance(mult, float):
            self.mult = mult
        else:
            raise TypeError("2rd arg mult should be float")

        if isinstance(theta_s, float):
            self.theta_s = theta_s
        else:
            raise TypeError("3th arg theta_s should be float")

    def setopls(self,k1,k2,k3,k4):
        """
        set opls parameters
        
        Args:
            * k1     (float) force constant    kcal/mol
            * k2     (float) force constant    kcal/mol
            * k3     (float) force constant    kcal/mol
            * k4     (float) force constant    kcal/mol

        .. math::
            E = 1/2k1[1+cos(theta)]

        .. math::
            +1/2k2[1-cos(2theta)]

        .. math::
            +1/2k3[1+cos(3theta)]

        .. math::
            +1/2k4[1-cos(4theta)]

        """

        if isinstance(k1, float):
            self.k1 = k1
        else:
            raise TypeError("1st arg should be float")

        if isinstance(k2, float):
            self.k2 = k2
        else:
            raise TypeError("2nd arg should be float")

        if isinstance(k3, float):
            self.k3 = k3
        else:
            raise TypeError("3rd arg should be float")

        if isinstance(k4, float):
            self.k4 = k4
        else:
            raise TypeError("4th arg should be float")

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

        Args:
            C0     (float) force constant    kcal/mol
            C1     (float) force constant    kcal/mol
            C2     (float) force constant    kcal/mol
            C3     (float) force constant    kcal/mol
            C4     (float) force constant    kcal/mol
            C5     (float) force constant    kcal/mol
            
        .. math::
            V_{rb}(theta) = \sum_{n=0}^5 C_n [ cos(theta - 180 )^n ]

        """

        if isinstance(C0, float):
            self.C0 = C0
        else:
            raise TypeError("1st arg should be float")

        if isinstance(C1, float):
            self.C1 = C1
        else:
            raise TypeError("2nd arg should be float")

        if isinstance(C2, float):
            self.C2 = C2
        else:
            raise TypeError("3rd arg should be float")

        if isinstance(C3, float):
            self.C3 = C3
        else:
            raise TypeError("4th arg should be float")

        if isinstance(C4, float):
            self.C4 = C4
        else:
            raise TypeError("5th arg should be float")

        if isinstance(C5, float):
            self.C5 = C5
        else:
            raise TypeError("6th arg should be float")
        
        # Translate to opls 
        self.k1 = -1.0*( 2.0*C1 + 3.0*C3/2.0)
        self.k2 = -1.0*( C2 + C4)
        self.k3 = -0.5*C3
        self.k4 = -0.25*C4

        

    def update_units(self,new_unit_conf):
        '''
        Update instance values with new units
        
        Args:
            * new_unit_conf (dict): with unit type as the key and the new unit as the value
            
        '''
        
        self._property,self._unit_conf = units.change_properties_units(self._unit_conf,new_unit_conf,self._property_units,self._property)
                