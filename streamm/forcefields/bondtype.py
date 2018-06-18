# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"



try:
    # Import pymatgen Class 
    import pymatgen_core.core.units as units 
except:
    raise ImportError("pymatgen import error for units object")


class Bondtype(units.ObjectUnits):
    """
    Set of bond parameters

    Kwargs:
        * fftype1  (str):   Atom type 
        * fftype2  (str):   Atom type 
        * type    (str):   Bond type
        * units_conf (dict): Dictionary of units for each attribute type
                
    """

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
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        #
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

        Args:
            * r0 (float) distance 
            * kb (float) force constant
    
        .. math::
            E = kb( r - r_0 )^2 

            
        """

        
        if isinstance(r0, float):
            self.r0 = r0
        else:
            raise TypeError("1st arg should be float")

        if isinstance(kb, float):
            self.kb = kb
        else:
            raise TypeError("2nd arg should be float")

    def export_json(self):
        '''    
        Export object to json
        
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        
        json_data = {}
        json_data['fftype1'] = self.fftype1
        json_data['fftype2'] = self.fftype2
        json_data['type'] = self.type
        json_data['r0'] = self.r0
        json_data['kb'] = self.kb
        json_data['lammps_index'] = self.lammps_index
        json_data['gromacs_index'] = self.gromacs_index
        # unit_conf
        json_data['unit_conf'] = self.unit_conf
        
        return json_data

    def import_json(self,json_data):
        '''    
        Export object to json
        
        Args:
            * json_data (dict) json representation of the object
            
        '''
        
        self.fftype1  =  json_data['fftype1']
        self.fftype2  =  json_data['fftype2']
        self.type  =  json_data['type']
        self.r0  =  json_data['r0']
        self.kb  =  json_data['kb']
        self.lammps_index  =  json_data['lammps_index']
        self.gromacs_index  =  json_data['gromacs_index']
        

        # Read in Unit config 
        if( 'unit_conf' in json_data.keys() ):               
            self._unit_conf = json_data['unit_conf']
        else:
            logger.warning('unit_conf not in json ')
                
        
 