# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals


__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"


"""
This module defines the classes relating to the improper dihedral angles between atoms
"""

# Import pymatgen module 
import pymatgen_core.core.units as units 

class Improper(units.ObjectUnits):
    """
    Data structure for describing any 4-point association of Particles
    
    Args:
        * pkey1   (int)   Dictionary key of Particle object in improper
        * pkey2   (int)   Dictionary key of Particle object in improper
        * pkey3   (int)   Dictionary key of Particle object in improper
        * pkey4   (int)   Dictionary key of Particle object in improper

    Kwargs:
        * units_conf (dict): Dictionary of units for each attribute type
                
    """

    def __init__(self, pkey1, pkey2, pkey3, pkey4,unit_conf=units.unit_conf ):
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        

        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            raise TypeError("1st arg should be int")

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            raise TypeError("2nd arg should be int type")

        if isinstance(pkey3, int):
            self.pkey3 = pkey3
        else:
            raise TypeError("3rd arg should be int type")

        if isinstance(pkey4, int):
            self.pkey4 = pkey4
        else:
            raise TypeError("4rd arg should be int type")
        #

        # Force field
        self.param = None
        # Lammps and gromacs index
        self.param_index = 0 
        self.lammps_index = 0 
        self.gromacs_index = 0 
        
    def __del__(self):
        del self.pkey1
        del self.pkey2 
        del self.pkey3 
        del self.pkey4
        del self.param
        del self.param_index
        del self.lammps_index
        del self.gromacs_index

    def __str__(self):
        return " %s - %s - %s - %s"%(self.pkey1,self.pkey2,self.pkey3,self.pkey4 )
        

    def export_json(self):
        '''    
        Export object to json
        
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        
        json_data = {}
        json_data['pkey1'] = self.pkey1
        json_data['pkey2'] = self.pkey2
        json_data['pkey3'] = self.pkey3
        json_data['pkey4'] = self.pkey4
        json_data['cosine'] = self.cosine
        json_data['param_index'] = self.param_index
        json_data['lammps_index'] = self.lammps_index
        json_data['gromacs_index'] = self.gromacs_index
        #         
        return json_data
        

    def import_json(self,json_data):
        '''    
        Export object to json
        
        Args:
            * json_data (dict) json representation of the object
            
        '''
        self.pkey1  =  json_data['pkey1']
        self.pkey2  =  json_data['pkey2']
        self.pkey3  =  json_data['pkey3']
        self.pkey4  =  json_data['pkey4']
        self.cosine  =  json_data['cosine']
        self.param_index  =  json_data['param_index']
        self.lammps_index  =  json_data['lammps_index']
        self.gromacs_index  =  json_data['gromacs_index']
