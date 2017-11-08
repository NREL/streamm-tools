# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import math
import itertools
import warnings

from six.moves import map, zip

import random 
import json
        
import numpy as np
from numpy.linalg import inv
from numpy import pi, dot, transpose, radians

from monty.json import MSONable
from monty.dev import deprecated

from pymatgen_core.util.num import abs_cap
import pymatgen_core.core.units as units 

import logging
logger = logging.getLogger(__name__)


"""
This module defines the classes relating to 3D lattices.

.. Changes ::
    * 
    * 
    * 
    * 
"""


__author__ = "Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import pymatgen_core.core.units as units



class Lattice(units.ObjectUnits):
    """
    A lattice object.  Essentially a matrix with conversion matrices. In
    general, it is assumed that length units are in Angstroms and angles are in
    degrees unless otherwise stated.
    """

    # Properties lazily generated for efficiency.

    @property
    def unit_conf(self):
        return self._unit_conf
    
    @property
    def matrix(self):
        return self._property['matrix'] 

    @property
    def constants(self):
        return self._property['constants']
        
    @property
    def angles(self):
        return self._property['angles'] 
    
    
    @property
    def lengths(self):
        return self._property['lengths'] 
    
    
    @matrix.setter
    def matrix(self,matrix):

        m = np.array(matrix, dtype=np.float64).reshape((self.n_dim, self.n_dim))
        lengths = np.sqrt(np.sum(m ** 2, axis=1))
        angles = np.zeros(self.n_dim)
        for i in range(self.n_dim):
            j = (i + 1) % self.n_dim
            k = (i + 2) % self.n_dim
            angles[i] = abs_cap(dot(m[j], m[k]) / (lengths[j] * lengths[k]))
                    
        self._property['lengths'] = lengths
        self._property['angles'] = np.arccos(angles) 
        
        
        if( self.unit_conf['angle'] == 'degree' ):
            self._property['angles'] = self._property['angles']* 180/np.pi
            
        self._property['constants'] = np.zeros(6)
        self._property['constants'][0:3] = self._property['lengths']
        self._property['constants'][3:6] = self._property['angles']
        
        self._property['matrix'] = m
        self.is_orthogonal = all([abs(a - 90) < 1e-5 for a in self.angles])
    

    def get_rad_angles(self):
        '''Return angles in radians
        
        '''

        if( self._unit_conf['angle'] == 'degree'):
            alpha_r = np.deg2rad(self._property['angles'][0] )
            beta_r  = np.deg2rad(self._property['angles'][1]  )
            gamma_r = np.deg2rad(self._property['angles'][2] )
        elif( self._unit_conf['angle'] == 'radian'):
            alpha_r = self._property['angles'][0]
            beta_r  = self._property['angles'][1]
            gamma_r = self._property['angles'][2]
        else:
            raise KeyError('angle unit {} not supported '.format( self._unit_conf['angle'] ))
                    
        return alpha_r,beta_r,gamma_r
        
    def constants2matrix(self,a,b,c,alpha_r,beta_r,gamma_r):
        '''Calculate the matrix based on the lattice constants
        
        Args:
            * a  (float): lattice constant a
            * b   (float): lattice constant b
            * c   (float): lattice constant c
            * alpha   (float): lattice angle alpha in radians
            * beta    (float): lattice angle beta in radians
            * gamma   (float): lattice angle gamma in radians
            
        '''
        # 
        val = (np.cos(alpha_r) * np.cos(beta_r) - np.cos(gamma_r))\
            / (np.sin(alpha_r) * np.sin(beta_r))
        # Sometimes rounding errors result in values slightly > 1.
        val = abs_cap(val)
        gamma_star = np.arccos(val)
        self._property['matrix'][0] = [a * np.sin(beta_r), 0.0, a * np.cos(beta_r)]
        self._property['matrix'][1] = [-b * np.sin(alpha_r) * np.cos(gamma_star),
                    b * np.sin(alpha_r) * np.sin(gamma_star),
                    b * np.cos(alpha_r)]
        self._property['matrix'][2] = [0.0, 0.0, float(c)]
            
    
    @angles.setter
    def angles(self,angles):

        self._property['angles'][0] = angles[0] 
        self._property['angles'][1] = angles[1] 
        self._property['angles'][2] = angles[2] 

        
        a = self._property['lengths'][0] 
        b = self._property['lengths'][1] 
        c = self._property['lengths'][2] 
        
        alpha_r,beta_r,gamma_r = self.get_rad_angles()
        self.constants2matrix(a,b,c,alpha_r,beta_r,gamma_r)


        
    @lengths.setter
    def lengths(self,lengths):

        
        self._property['lengths'][0] = lengths[0] 
        self._property['lengths'][1] = lengths[1] 
        self._property['lengths'][2] = lengths[2] 

        a = lengths[0] 
        b = lengths[1] 
        c = lengths[2]
        
        alpha_r,beta_r,gamma_r = self.get_rad_angles()       
        self.constants2matrix(a,b,c,alpha_r,beta_r,gamma_r)


    @property
    def constants(self):

        self._property['constants'][0:3] = self._property['lengths']
        self._property['constants'][3:6] = self._property['angles']
        
        return self._property['constants']
    
    @constants.setter
    def constants(self,constants):
        """Set lattice with lattice constants.  
        
        Args:
            * constants[0]     (float): lattice constant a
            * constants[1]     (float): lattice constant b
            * constants[2]     (float): lattice constant c
            * constants[3]     (float): lattice angle alpha 
            * constants[4]     (float): lattice angle beta
            * constants[5]     (float): lattice angle gamma
            
        """
        
        self._property['constants'] = constants

        self._property['lengths'][0] = constants[0] 
        self._property['lengths'][1] = constants[1] 
        self._property['lengths'][2] = constants[2] 

        self._property['angles'][0] = constants[3] 
        self._property['angles'][1] = constants[4] 
        self._property['angles'][2] = constants[5] 
        
        a = constants[0] 
        b = constants[1] 
        c = constants[2]
        
        alpha_r,beta_r,gamma_r = self.get_rad_angles()        
        self.constants2matrix(a,b,c,alpha_r,beta_r,gamma_r)
                
        
    def __init__(self, matrix=[100.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,100.0],unit_conf = units.unit_conf):
        """
        Create a lattice from any sequence of 9 numbers. Note that the sequence
        is assumed to be read one row at a time. Each row represents one
        lattice vector.

        Args:
            * matrix: Sequence of numbers in any form. Examples of acceptable
                input.
                i) An actual numpy array.
                ii) [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                iii) [1, 0, 0 , 0, 1, 0, 0, 0, 1]
                iv) (1, 0, 0, 0, 1, 0, 0, 0, 1)
                Each row should correspond to a lattice vector.
                E.g., [[10, 0, 0], [20, 10, 0], [0, 0, 30]] specifies a lattice
                with lattice vectors [10, 0, 0], [20, 10, 0] and [0, 0, 30].
        """
        # init object's units dictionaries 
        units.ObjectUnits.__init__(self,unit_conf=unit_conf)
        # 
        self.n_dim = int(3) # Number of spatial dimensions
        # 
        self.matrix = matrix
        # Set the unit types for each property with units 
        self._property_units['length'].append('matrix')
        self._property_units['length'].append('lengths')
        self._property_units['angle'].append('angles')
         
        
        self._inv_matrix = None
        self._metric_tensor = None
        self._diags = None
        self._lll_matrix_mappings = {}
        self._lll_inverse = None
        
        self.pbcs = [ False for d in range(self.n_dim) ] # Initialize Periodic boundary conditions as off
 
    def __del__(self):
        del self.n_dim
        del self.pbcs


        
        
    def __format__(self, fmt_spec=''):
        """
        Support format printing. Supported formats are:

        1. "l" for a list format that can be easily copied and pasted, e.g.,
           ".3fl" prints something like
           "[[10.000, 0.000, 0.000], [0.000, 10.000, 0.000], [0.000, 0.000, 10.000]]"
        2. "p" for lattice parameters ".1fp" prints something like
           "{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}"
        3. Default will simply print a 3x3 matrix form. E.g.,
           10.000 0.000 0.000
           0.000 10.000 0.000
           0.000 0.000 10.000
        """
        m = self.matrix.tolist()
        if fmt_spec.endswith("l"):
            fmt = "[[{}, {}, {}], [{}, {}, {}], [{}, {}, {}]]"
            fmt_spec = fmt_spec[:-1]
        elif fmt_spec.endswith("p"):
            fmt = "{{{}, {}, {}, {}, {}, {}}}"
            fmt_spec = fmt_spec[:-1]
            m = self.lengths_and_angles
        else:
            fmt = "{} {} {}\n{} {} {}\n{} {} {}"
        return fmt.format(*[format(c, fmt_spec) for row in m
                            for c in row])


    def __str__(self):
        return "\n".join([" ".join(["%.6f" % i for i in row])
                          for row in self.matrix])


    def copy(self):
        """Deep copy of self."""
        return self.__class__(self.matrix.copy())
        

    def fractoreal(self,frac_o):
        '''
        Translate fractional coordinates to real
        
        Args:
            * frac_o (numpy array): fraction coordinates
            
        Returns:
            * pos_o (numpy array): real cartesian coordinates
            
        '''
        pos_o = np.zeros(self.n_dim)
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                pos_o[m] += self.matrix[n][m]*frac_o[n]
                        
        return pos_o
            
        
    def d_pos(self,pos_i,pos_j):
        """
        Difference between two positions
        
        Args:
            * pos_i (list): cartesian coordinates [x,y,z] for position i
            * pos_j (list): cartesian coordinates [x,y,z] for position j

        Returns:
            * dr_ij (list): cartesian vector [x,y,z] between position i and j
            
        If pbcs are set to true this calculates the fractional coordinates according to the
        inverse of the matrix multiplied by the position 
        
        .. math ::
        
            F = P M^{-1}
        
        Then subtractions the fractional coordinates and rounds them to be within the lattice box.
        Then translates the fractional coordinates back to real coordinates to get the separation vector
        
        """
        if( any ( pbcs_i for pbcs_i in self.pbcs ) ):
            frac_i = np.zeros(self.n_dim)
            frac_j = np.zeros(self.n_dim)
            frac_ij = np.zeros(self.n_dim)
            d_pos = np.zeros(self.n_dim)
    
            AA = self.matrix[0][0]
            AB = self.matrix[0][1]
            AC = self.matrix[0][2]
            BA = self.matrix[1][0]
            BB = self.matrix[1][1]
            BC = self.matrix[1][2]
            CA = self.matrix[2][0]
            CB = self.matrix[2][1]
            CC = self.matrix[2][2]
            
            X = pos_i[0]
            Y = pos_i[1]
            Z = pos_i[2]
            frac_i[0] = -((-BC*CB*X+BB*CC*X+BC*CA*Y-BA*CC*Y-BB*CA*Z+BA*CB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
            frac_i[1] = -((AC*CB*X-AB*CC*X-AC*CA*Y+AA*CC*Y+AB*CA*Z-AA*CB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
            frac_i[2] = -((-AC*BB*X+AB*BC*X+AC*BA*Y-AA*BC*Y-AB*BA*Z+AA*BB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
            X = pos_j[0]
            Y = pos_j[1]
            Z = pos_j[2]
            frac_j[0] = -((-BC*CB*X+BB*CC*X+BC*CA*Y-BA*CC*Y-BB*CA*Z+BA*CB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
            frac_j[1] = -((AC*CB*X-AB*CC*X-AC*CA*Y+AA*CC*Y+AB*CA*Z-AA*CB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
            frac_j[2] = -((-AC*BB*X+AB*BC*X+AC*BA*Y-AA*BC*Y-AB*BA*Z+AA*BB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
            
            frac_ij[0] = frac_j[0] - frac_i[0]
            frac_ij[1] = frac_j[1] - frac_i[1]
            frac_ij[2] = frac_j[2] - frac_i[2]
            
            if( self.pbcs[0] ): frac_ij[0] = frac_ij[0] - round(frac_ij[0])
            if( self.pbcs[1] ): frac_ij[1] = frac_ij[1] - round(frac_ij[1])
            if( self.pbcs[2] ): frac_ij[2] = frac_ij[2] - round(frac_ij[2])
                    
            d_pos[0] = AA*frac_ij[0] + BA*frac_ij[1] + CA*frac_ij[2]
            d_pos[1] = AB*frac_ij[0] + BB*frac_ij[1] + CB*frac_ij[2]
            d_pos[2] = AC*frac_ij[0] + BC*frac_ij[1] + CC*frac_ij[2]
            
        else:
            d_pos  = np.array(pos_j) -  np.array(pos_i)
            
        return  d_pos


    def delta_pos(self,pos_i,pos_j):
        """
        Difference between two positions

        Args:
            * pos_i (list): cartesian coordinates [x,y,z] for position i
            * pos_j (list): cartesian coordinates [x,y,z] for position j        

        Returns:
            * dr_ij (list): cartesian vector [x,y,z] between position i and j
            * mag_dr_ij (float): magnitude of vector dr_ij
            
        """
        
        dr_ij  = self.d_pos(pos_i,pos_j)
        mag_dr_ij = np.sqrt(dr_ij.dot(dr_ij))

        return dr_ij,mag_dr_ij
            

    def norm_d_pos(self,pos_i,pos_j):
        '''
        Normalized difference between two positions in cubic lattice
        
        Args:
            * pos_i (list): cartesian coordinates [x,y,z] for position i
            * pos_j (list): cartesian coordinates [x,y,z] for position j
            
        Returns:
            * dr_ij (list): normalized cartesian vector [x,y,z] between position i and j
            
        '''
        dr_ij = self.d_pos(pos_i, pos_j)
        return (dr_ij)/np.linalg.norm(dr_ij)
        

            
    def delta_npos(self,npos_i,npos_j):
        """
        Difference between two position lists in  cubic lattice
        
        Args:
            * pos_i (list of lists): list of cartesian coordinates [x,y,z] for positions i
            * pos_j (list of lists): list of cartesian coordinates [x,y,z] for positions j        

        Returns:
            * npos_ij (m x n numpy array): array of cartesian vector [x,y,z] between position i and j where the npos_ij[i][j] value is the difference between pos_i[i] and pos_j[j] 
            * mag_dr_ij (m x n numpy array): array of the magnitudes of vector dr_ij where the mag_dr_ij[i][j] value is the distance between pos_i[i] and pos_j[j]
            
        """
        n_i = len(npos_i)
        n_j = len(npos_j)
        key_list_i = range(n_i)
        key_list_j = range(n_j)
        #
        n_ij = n_i*n_j
        #
        logger.debug(" Taking difference %d x %d "%(n_i,n_j))
        
        npos_ij = np.zeros(n_ij*self.n_dim,dtype='float64')
        #npos_ij = np.zeros(shape=(n_ij,self.n_dim),dtype='float64')
        nd_ij=  np.zeros(shape=(n_i,n_j),dtype='float64')

        # if any pbc's are on 
        if( any ( pbcs_i for pbcs_i in self.pbcs ) ):
            logging.debug(" Using cubic pbcs")
            for m in key_list_i:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.d_pos(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    #npos_ij[m][n] = dr_ij
                    nd_ij[m][n] = np.sqrt(dot_dr_ij)
        else:
            logging.debug(" Using no pbcs")
            for m in key_list_i:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.d_pos(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    #npos_ij[m][n] = dr_ij
                    nd_ij[m][n] = np.sqrt(dot_dr_ij)
          
        return npos_ij,nd_ij


    def proximitycheck(self,npos_i,npos_j,pos_cut,p=None):
        """
        Difference between two position lists in  cubic lattice

        Args:
            * pos_i (list of lists): list of cartesian coordinates [x,y,z] for positions i
            * pos_j (list of lists): list of cartesian coordinates [x,y,z] for positions j
            
        Returns:
            * True: if no overlap between particles was found
            * False: if overlap between particles was found
            
        
        """
        n_i = len(npos_i)
        n_j = len(npos_j)
        key_list_i = range(n_i)
        key_list_j = range(n_j)
        pos_cut_sq = pos_cut*pos_cut
        #
        # MPI setup
        #
        if( p == None ):
            rank = 0
            size = 0
            key_list_i_p = key_list_i
        else:
            rank = p.getRank()
            size = p.getCommSize()
            key_list_i_p =  p.splitListOnProcs(key_list_i)
            
        overlap_p = 0
        
        # if any pbc's are on 
        for m in key_list_i_p:
            pos_i = npos_i[m]
            for n in key_list_j:
                pos_j = npos_j[n]
                dr_ij = self.d_pos(pos_i,pos_j)
                dot_dr_ij = dr_ij.dot(dr_ij)
                if( dot_dr_ij < pos_cut_sq ):
                    overlap_p += 1
                    
        logger.debug(overlap_p)
        if( p == None ):
            if( overlap_p == 0 ):
                return True
            else:
                return False
        else:
            overlap_sum = p.allReduceSum(overlap_p)
            if( overlap_sum == 0 ):
                return True
            else:
                return False
                                
                
    def random_pos(self):
        '''
        Generate random position in lattice.
        
        Assume random.seed(seed) has been initialized
        
        '''
        
        n_dec = int(6)
        int_mult = 10**float(n_dec)
        pos_o = np.zeros(self.n_dim)
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                if( self.matrix[m][n] > 0.0 ):
                    pos_o[m] += float(random.randrange(0,int(int_mult*self.matrix[m][n]) ))/float(int_mult)
                    
        return pos_o


    def expand_matrix(self,exlat_frac):
        '''
        Increase size of lattice by certain fraction

        Args:
            * exlat_frac (float): fraction to increase lattice by
            
        '''
        matrix_i = self.matrix
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                matrix_i[m][n] += matrix_i[m][n]*exlat_frac
                
        self.matrix = matrix_i
        
        return 
        


    def set_cubic(self,len_o):
        '''
        Set lattice to cubic with lattice constant
        
        Args:
            * len_o (float): Length of cubic lattice constant a
        
        '''
        matrix = np.zeros((self.n_dim, self.n_dim))
        for d in range(self.n_dim):
            matrix[d][d] = len_o
        self.matrix = matrix        
                    

        
        
    def export_json(self,tag,write_file=True):
        '''    
        Export object to json
        
        Args:
            * tag (str) ID of file to be written
        Kwargs:
            * write_file (boolean) to dump json to a file
            
        Returns:
            * json_data (dict) json representation of the object
            
        '''
        #
        json_data = {}
        json_data['unit_conf'] = self.unit_conf
        json_data['matrix'] = [ ]
        for x in self.matrix:
            for y in x:
                json_data['matrix'].append(y)
        #
        json_data['pbcs'] = self.pbcs
        #
        if( write_file ):
            with open("%s_lat.json"%(tag),'wb') as fl:
                json.dump(json_data,fl,indent = 2)
        #
        return json_data

    def import_json(self,tag,json_data={},read_file=True):
        '''    
        Export object to json
        
        Args:
            * tag (str) ID of file to be read
            
        Kwargs:
            * read_file (boolean) to read json from a file
            * json_data (dict) json representation of the object
            
        '''
        # 
        if( read_file ):
            file_name = "%s_lat.json"%(tag)
            print("Reading {}".format(file_name))
            with open(file_name,'rb') as fl:
                json_data = json.load(fl)
        #
        logger.debug("Set object properties based on json")
        
        # Read in Unit config 
        if( 'unit_conf' in json_data.keys() ):               
            self._unit_conf = json_data['unit_conf']
        else:
            logger.warning('unit_conf not in json ')
        # Read in pbcs
        if( 'matrix' in json_data.keys() ):               
            self.matrix  = json_data['matrix']
        else:
            logger.warning('matrix not in json ')
        
        # Read in pbcs
        if( 'pbcs' in json_data.keys() ):               
            self.pbcs  = json_data['pbcs']
        else:
            logger.warning('pbcs not in json ')
        #
        #
        #
        return 