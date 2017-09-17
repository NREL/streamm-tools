# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import math
import itertools
import warnings

from six.moves import map, zip

import random 
 
import numpy as np
from numpy.linalg import inv
from numpy import pi, dot, transpose, radians

from monty.json import MSONable
from monty.dev import deprecated

from pymatgen_core.util.num import abs_cap


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



class Lattice(MSONable):
    """
    A lattice object.  Essentially a matrix with conversion matrices. In
    general, it is assumed that length units are in Angstroms and angles are in
    degrees unless otherwise stated.
    """

    # Properties lazily generated for efficiency.


    def __init__(self, matrix=[100.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,100.0]):
        """
        Create a lattice from any sequence of 9 numbers. Note that the sequence
        is assumed to be read one row at a time. Each row represents one
        lattice vector.

        Args:
            matrix: Sequence of numbers in any form. Examples of acceptable
                input.
                i) An actual numpy array.
                ii) [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                iii) [1, 0, 0 , 0, 1, 0, 0, 0, 1]
                iv) (1, 0, 0, 0, 1, 0, 0, 0, 1)
                Each row should correspond to a lattice vector.
                E.g., [[10, 0, 0], [20, 10, 0], [0, 0, 30]] specifies a lattice
                with lattice vectors [10, 0, 0], [20, 10, 0] and [0, 0, 30].
        """
        self.n_dim = int(3) # Number of spatial dimensions
        self.set_matrix(matrix)
        self.pbcs = [ False for d in range(self.n_dim) ] # Initialize periodic boundries as off


    def __del__(self):
        del self.n_dim
        del self.pbcs
        
    def set_matrix(self,matrix):
        '''
        Set the values of the matrix of lattice vectors
        and reset lengths and angles accordingly
        
        Args:
            matrix (list) of length n_dim x n_dim
        
        
        '''        
        m = np.array(matrix, dtype=np.float64).reshape((self.n_dim, self.n_dim))
        lengths = np.sqrt(np.sum(m ** 2, axis=1))
        angles = np.zeros(self.n_dim)
        for i in range(self.n_dim):
            j = (i + 1) % self.n_dim
            k = (i + 2) % self.n_dim
            angles[i] = abs_cap(dot(m[j], m[k]) / (lengths[j] * lengths[k]))

        self._angles = np.arccos(angles) * 180. / pi
        self._lengths = lengths
        self._matrix = m
        self._inv_matrix = None
        self._metric_tensor = None
        self._diags = None
        self._lll_matrix_mappings = {}
        self._lll_inverse = None
        self.is_orthogonal = all([abs(a - 90) < 1e-5 for a in self._angles])

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

    def copy(self):
        """Deep copy of self."""
        return self.__class__(self.matrix.copy())



    def set_consts(self,a, b, c, alpha, beta, gamma):
        """Set lattice with lattice constants.  
        
        Args:
            a     (float): lattice constant a
            b     (float): lattice constant b
            c     (float): lattice constant c
            alpha (float): lattice angle alpha
            beta  (float): lattice angle beta
            gamma (float): lattice angle gamma
            
        """
        alpha_rad = np.deg2rad(alpha)
        beta_rad  = np.deg2rad(beta )
        gamma_rad = np.deg2rad(gamma)

        ax = a
        bx = np.cos( gamma_rad )* b
        by = np.sin( gamma_rad )* b
        cx = c*np.cos(beta_rad)
        cy = c*(np.cos(alpha_rad)- np.cos(beta)*np.cos(gamma_rad))/np.sin(gamma_rad)
        cz = np.sqrt( c*c - cx*cx - cy*cy )

        self._lengths[0] = a
        self._lengths[1] = b
        self._lengths[2] = c

        self._angles[0] = alpha
        self._angles[1] = beta
        self._angles[2] = gamma
        
        
        # Set lattice vectors
        self._matrix =  np.array([ np.zeros(self.n_dim) for dim in range(self.n_dim) ])
        self._matrix[0][0] = ax
        self._matrix[1][0] = bx
        self._matrix[1][1] = by
        self._matrix[2][0] = cx
        self._matrix[2][1] = cy
        self._matrix[2][2] = cz
        
        return         
        
    def set_box(self,box):
        '''
        Set lattice constants based on box list.
        
        Args:
            box (list): List of lattice constants (a, b, c, alpha, beta, gamma)
        '''
        if( len(box) != 6 ):
            logger.warning("box variable does not length 6 for values a, b, c, alpha, beta, gamma. The Lattice will not be set ")
            return
        
        a      = box[0]
        b      = box[1]
        c      = box[2]
        alpha  = box[3]
        beta   = box[4]
        gamma  = box[5]
        
        self.set_consts(a, b, c, alpha, beta, gamma)
        
        return 
        
    def deltasq_pos(self,pos_i,pos_j):
        """
        Difference between two positions
        
        Args:
            pos_i (list): cartesian coordinates [x,y,z] for position i
            pos_j (list): cartesian coordinates [x,y,z] for position j

        Returns:
            dr_ij (list): cartesian vector [x,y,z] between position i and j
            
        """
        
        dr_ij  = np.array(pos_j) -  np.array(pos_i)

        return  dr_ij

    def deltasq_pos_c(self,pos_i,pos_j):
        """
        Difference between two positions in cubic lattice

        Args:
            pos_i (list): cartesian coordinates [x,y,z] for position i
            pos_j (list): cartesian coordinates [x,y,z] for position j
        Returns:
            dr_ij (list): cartesian vector [x,y,z] between position i and j
        """
        dr_ij = self.deltasq_pos(pos_i,pos_j)
        
        dr_ij[0] = dr_ij[0] - self._matrix[0][0] * round( dr_ij[0]/  self._matrix[0][0] )
        dr_ij[1] = dr_ij[1] - self._matrix[1][1] * round( dr_ij[1]/  self._matrix[1][1] )
        dr_ij[2] = dr_ij[2] - self._matrix[2][2] * round( dr_ij[2]/  self._matrix[2][2] )
        
        return dr_ij

    def norm_delta_pos_c(self,pos_i,pos_j):
        '''
        Normalized difference between two positions in cubic lattice
        

        Args:
            pos_i (list): cartesian coordinates [x,y,z] for position i
            pos_j (list): cartesian coordinates [x,y,z] for position j
            
        Returns:
            dr_ij (list): normalized cartesian vector [x,y,z] between position i and j
            
        '''
        dr_ij = self.deltasq_pos_c(pos_i, pos_j)
        return (dr_ij)/np.linalg.norm(dr_ij)
        
        

    def delta_pos_c(self,pos_i,pos_j):
        """
        Difference between two positions in cubic lattice

        Args:
            pos_i (list): cartesian coordinates [x,y,z] for position i
            pos_j (list): cartesian coordinates [x,y,z] for position j

        Returns:
            dr_ij (list): cartesian vector [x,y,z] between position i and j
            mag_dr_ij (float): magnitude of vector dr_ij    
        """
        dr_ij = self.deltasq_pos_c(pos_i,pos_j)
        mag_dr_ij = np.sqrt(dr_ij.dot(dr_ij))

        return dr_ij,mag_dr_ij


    def delta_pos(self,pos_i,pos_j):
        """
        Difference between two positions

        Args:
            pos_i (list): cartesian coordinates [x,y,z] for position i
            pos_j (list): cartesian coordinates [x,y,z] for position j        

        Returns:
            dr_ij (list): cartesian vector [x,y,z] between position i and j
            mag_dr_ij (float): magnitude of vector dr_ij        """
        
        dr_ij  = self.deltasq_pos(pos_i,pos_j)
        mag_dr_ij = np.sqrt(dr_ij.dot(dr_ij))

        return dr_ij,mag_dr_ij
            
    def delta_npos(self,npos_i,npos_j):
        """
        Difference between two position lists in  cubic lattice
        
        Args:
            pos_i (list of lists): list of cartesian coordinates [x,y,z] for positions i
            pos_j (list of lists): list of cartesian coordinates [x,y,z] for positions j        

        Returns:
            npos_ij (m x n numpy array): array of cartesian vector [x,y,z] between position i and j
            where the npos_ij[i][j] value is the difference between pos_i[i] and pos_j[j] 
            mag_dr_ij (m x n numpy array): array of the magnitudes of vector dr_ij      
            where the mag_dr_ij[i][j] value is the distance between pos_i[i] and pos_j[j] 
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
                    dr_ij = self.deltasq_pos_c(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    #npos_ij[m][n] = dr_ij
                    nd_ij[m][n] = np.sqrt(dot_dr_ij)
        else:
            logging.debug(" Using no pbcs")
            for m in key_list_i:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.deltasq_pos(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    #npos_ij[m][n] = dr_ij
                    nd_ij[m][n] = np.sqrt(dot_dr_ij)
          
        return npos_ij,nd_ij


    def proximitycheck(self,npos_i,npos_j,pos_cut,p=None):
        """
        Difference between two position lists in  cubic lattice

        Args:
            pos_i (list of lists): list of cartesian coordinates [x,y,z] for positions i
            pos_j (list of lists): list of cartesian coordinates [x,y,z] for positions j
            
        Returns:
            True: if no overlap between particles was found
            False: if overlap between particles was found
            
        .. Todo::
            move to utilities 
        
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
        if( any ( pbcs_i for pbcs_i in self.pbcs ) ):
            logging.debug(" Using cubic pbcs")
            for m in key_list_i_p:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.deltasq_pos_c(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    if( dot_dr_ij < pos_cut_sq ):
                        overlap_p += 1 
        else:
            logging.debug(" Using no pbcs")
            for m in key_list_i_p:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.deltasq_pos(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    logger.debug(pos_i)
                    logger.debug(pos_j)
                    logger.debug(pos_cut_sq)
                    logger.debug(" {} - {} : {} ".format(m,n,dot_dr_ij))
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
                if( self._matrix[m][n] > 0.0 ):
                    pos_o[m] += float(random.randrange(0,int(int_mult*self._matrix[m][n]) ))/float(int_mult)
                    
        return pos_o


    def expand_matrix(self,exlat_frac):
        '''
        Increase size of lattice by certain fraction

        Args:
            exlat_frac (float): fraction to increase lattice by
            
        '''
        matrix_i = self._matrix
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                matrix_i[m][n] += matrix_i[m][n]*exlat_frac
                
        self.set_matrix(matrix_i)
        
        return 
        

    def fractoreal(self,frac_o):
        '''
        Translate fractional coordinates to real
        
        Args:
            frac_o (numpy array): fraction coordinates
            
        Returns:
            pos_o (numpy array): real cartesian coordinates
            
        '''
        pos_o = np.zeros(self.n_dim)
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                pos_o[m] += self._matrix[n][m]*frac_o[n]
                        
        return pos_o
    

    def set_cubic(self,len_o):
        '''
        Set lattice to cubic with lattice constant
        
        Args:
            len_o (float): Length of cubic lattice constant a
        
        '''
        matrix = np.zeros((self.n_dim, self.n_dim))
        for d in range(self.n_dim):
            matrix[d][d] = len_o
        self.set_matrix(matrix)        
                    

    def update_units(self,new_unit_conf):
        '''
        Update instance values with new units
        
        Args:
            new_unit_conf (dict): with unit type as the key and the new unit as the value
            
        '''
        
        self._property,self._unit_conf = units.change_properties_units(self._unit_conf,new_unit_conf,self._property_units,self._property)
        