# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

"""
This module defines the classes relating to a periodic lattice

"""

import logging
logger = logging.getLogger(__name__)

try:
    from pymatgen.core.lattice import Lattice as pymatgen_lat
except:
    logger.warning("pymatgen import error for Lattice object")
    exit()
    
import numpy as np 


PRECISION = 8


class Lattice(pymatgen_lat):
    
    def __init__(self,matrix=[100.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,100.0]):
        """        
        Create a lattice 
        """
        pymatgen_lat.__init__(self, matrix=matrix)
        #
        self.pbcs = [ False for d in range(self.n_dim) ] # Initialize periodic boundries as off


    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.n_dim
        del self.pbcs 
        #
        
    def set_percision(self):
        '''
        Set values to correct precision 
        '''
        logger.info("Setting all values to {} decimal places ")
        for i in range(self.n_dim):
            self._lengths[i] = round(self._lengths[i],PRECISION)
            self._angles[i] = round(self._angles[i],PRECISION)
            for j in range(self.n_dim):
                self._matrix[i][j]  = round(self._matrix[i][j] ,PRECISION)

        return 

    def set_consts(self,a, b, c, alpha, beta, gamma):
        """
        Convert lattice constants  to lattice  vectors.  
        
        Args:
            math:: box = (list) [a,b,c,\alpha (degree),\beta (degree),\gamma (degree)]
        
        math::
        
        matrix_{0,0} = a
        matrix_{1,0} = cos(\gamma) x b
        matrix_{1,1} = sin(\gamma) x b
        matrix_{2,0} = 
        matrix_{2,1} = 
        matrix_{2,2} = 
        matrix_{0,0} = 
        
        Need to double check from materials book
        
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
        
        # Round values to correct precision     
        self.set_percision()
        
        return         
        
    def set_box(self,box):
        '''
        Set lattice constants based on box list
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
        """
        
        dr_ij  = np.array(pos_j) -  np.array(pos_i)

        return  dr_ij

    def deltasq_pos_c(self,pos_i,pos_j):
        """
        Difference between two positions in cubic lattice  
        """
        dr_ij = self.deltasq_pos(pos_i,pos_j)
        
        dr_ij[0] = dr_ij[0] - self._matrix[0][0] * round( dr_ij[0]/  self._matrix[0][0] )
        dr_ij[1] = dr_ij[1] - self._matrix[1][1] * round( dr_ij[1]/  self._matrix[1][1] )
        dr_ij[2] = dr_ij[2] - self._matrix[2][2] * round( dr_ij[2]/  self._matrix[2][2] )
        
        return dr_ij

    def norm_delta_pos_c(self,pos_i,pos_j):
        '''
        Normalized difference between two positions in cubic lattice 
        '''
        dr_ij = self.deltasq_pos_c(pos_i, pos_j)
        return (dr_ij)/np.linalg.norm(dr_ij)
        
        

    def delta_pos_c(self,pos_i,pos_j):
        """
        Difference between two positions in cubic lattice  
        """
        dr_ij = self.deltasq_pos_c(pos_i,pos_j)
        mag_dr_ij = np.sqrt(dr_ij.dot(dr_ij))

        return dr_ij,mag_dr_ij


    def delta_pos(self,pos_i,pos_j):
        """
        Difference between two positions 
        """
        
        dr_ij  = self.deltasq_pos(pos_i,pos_j)
        mag_dr_ij = np.sqrt(dr_ij.dot(dr_ij))

        return dr_ij,mag_dr_ij


    def proximitycheck(self,npos_i,npos_j,pos_cut,p=None):
        """
        Difference between two position lists in  cubic lattice  
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
                    pos_j = npos_i[n]
                    dr_ij = self.deltasq_pos(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    if( dot_dr_ij < pos_cut_sq ):
                        overlap_p += 1 
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
                                
                            
    def delta_npos(self,npos_i,npos_j):
        """
        Difference between two position lists in  cubic lattice  
        """
        n_i = len(npos_i)
        n_j = len(npos_j)
        key_list_i = range(n_i)
        key_list_j = range(n_j)
        #
        n_ij = n_i*n_j
        #print n_i,n_j
        #
        start_dpos = datetime.now()
        logging.debug(" Taking difference %d x %d "%(n_i,n_j))
        
        npos_ij = np.empty(n_ij*self.n_dim,dtype='float64')
        #nd_ij = np.empty(n_ij,dtype='float64')
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

        finish_dpos = datetime.now()
        delt_t = finish_dpos - start_dpos
        logging.debug(" Finished taking difference %s sec "%(delt_t.seconds))

          
        return npos_ij,nd_ij

    def random_pos(self):
        '''
        Generate random position in lattice

        random.seed(seed) needs to be initialized
        
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

        Arguments:
        exlat_frac (float) fraction to increase lattice by
        '''
        matrix_i = self._matrix
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                matrix_i[m][n] += matrix_i[m][n]*exlat_frac
                
        self.set_matrix(matrix_i)
        
        return 
        # 
    def fractoreal(self,frac_o):
        '''
        Translate fractional coordinates to real 

        Arguments:
            frac_o (np.array) fraction coordinates
        '''
        pos_o = np.zeros(self.n_dim)
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                pos_o[m] += self._matrix[n][m]*frac_o[n]
                        
        return pos_o
                  