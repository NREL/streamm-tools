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

class Imptype(units.ObjectUnits):
    """
    Set of improper dihedral angle parameters

    Args:
         * fftype1  (str)   Atom type 
         * fftype2  (str)   Atom type 
         * fftype3  (str)   Atom type 
         * fftype4  (str)   Atom type 
         * type    (str)   Bond type
         
    """

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
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)

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
        del self.fftype1
        del self.fftype2 
        del self.fftype3
        del self.fftype4
        del self.type
        del self.gromacs_index 
        del self.lammps_index 

        
    def __str__(self):
        strucStr =  " improper  %s - %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.fftype4,self.type)
        
        if( self.type ==  "improper" ):
            strucStr += "\n  imp e0 = %f ke = %f lammps index %d  gromcas index %d " %(self.e0,self.ke,self.lammps_index ,self.gromacs_index )

        return strucStr        


    def setimp(self, e0, ke):
        """
        set Harmonic parameters
        
        .. math::
            E = kb ( e_{lijk} - e_0 )^2 

        Args:
            * e0     (float) 
            * kb     (float) force constant    (energy)
        """

        if isinstance(e0, float):
            self.e0 = e0
        else:
            raise TypeError("1st arg should be float")

        if isinstance(ke, float):
            self.ke = ke
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
        json_data['fftype3'] = self.fftype3
        json_data['fftype4'] = self.fftype4
        json_data['type'] = self.type
        json_data['e0'] = self.e0
        json_data['pn'] = self.pn
        json_data['ke'] = self.ke
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
        # 
        self.fftype1  =  json_data['fftype1']
        self.fftype2  =  json_data['fftype2']
        self.fftype3  =  json_data['fftype3']
        self.fftype4  =  json_data['fftype4']
        self.type  =  json_data['type']
        self.e0  =  json_data['e0']
        self.pn  =  json_data['pn']
        self.ke  =  json_data['ke']
        self.lammps_index  =  json_data['lammps_index']
        self.gromacs_index  =  json_data['gromacs_index']



        # Read in Unit config 
        if( 'unit_conf' in json_data.keys() ):               
            self._unit_conf = json_data['unit_conf']
        else:
            logger.warning('unit_conf not in json ')