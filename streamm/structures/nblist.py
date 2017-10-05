# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

'''
This module defines the classes relating to neighbor lists 
'''

import logging
logger = logging.getLogger(__name__)
     
import csv
              
class NBlist(object):
    """Class for neighbor list.

    .. attribute:: list (list) 

        list of all neighboring indexes in sequence
        
    .. attribute:: index (list) 
        
        list of the location in ``list`` of the beginning of the neighbors sequence
        
    Example:
    You have ethane with particle labels of ``C1``, ``H1``, ``H2``, ``H3``, ``C2``, ``H4``, ``H5`` and ``H6``\
    with particle indexes of  0, 1, 2, 3, 4, 5, 6 and 7 respectively.
    
    ::
        
           H1         H4
             \       /
        H2 -- C1 -- C2 -- H5
             /       \ 
           H3         H6
    
    
    
    You would have the following bonds:
    
    * C1[0] - C2[4]
    * C1[0] - H1[1]
    * C1[0] - H2[2]
    * C1[0] - H3[3]
    * C2[4] - H4[5]
    * C2[4] - H5[6]
    * C2[4] - H6[7]

    The ``list`` would be ``[1, 2, 3, 4, 0, 0, 0, 0, 5, 6, 7, 4, 4, 4]`` and the ``index`` would be ``[0, 4, 5, 6, 7, 11, 12, 13, 14]``.
    So if we want the neighbor indexes of ``C2`` (particle ``4``),
    we get the beginning position in ``list`` from ``index[4]`` to be 7 and the ending position in ``list``
    for  ``index[7+1]-1`` to be 10. So from the ``index`` we know ``C2`` has the neighbors ``list[7:10]`` ([0, 5, 6, 7]),
    which are ``C2``, ``H4``, ``H5`` and ``H6``. The last value of 14 in the ``index`` list
    is to give the last position in the ``list`` for the last particle (``index[7+1]-1``). 
    """
    def __init__(self):
        self.list = []
        self.index = []
        self.cnt = -1 

    def __del__(self):
        del self.list
        del self.index
        del self.cnt
        
    def __str__(self):
        return " NBlist of {} particle with {} connections".format(len(self.index)-1,len(self.list))
    

    def calc_nnab(self,key_i):
        '''
        Calculate the number of neighbors for particle key_i
        
        Args:
            key_i (int): Index of item in neighbor list
        Returns:
            (int): number of neighbors of key_i        
        '''
        try:
            return self.index[key_i+1] - self.index[key_i]
        except:
            raise KeyError(" Neighbor list not set ")
            return []
            
    def getnbs(self,key_i):
        '''
        Return list of in of neighbors of key_i
        
        Args:
            key_i (int): Index of item in neighbor list
        Returns:
            nb_list (list): list of indexes of neighbors of key_i
        '''
        try:
            return self.list[self.index[key_i]:self.index[key_i+1]]
        except:
            raise KeyError(" Neighbor list not set ")
            return []
            

    def radii_nblist(self,lat,positions,radii,radii_buffer=1.25,write_dr=True,del_drmatrix=False):
        '''Create neighbor list of positions based on distance and radius of each particle 
        
        Args:
            * lat (Lattice) object 
            * positions (list) of  particle position numpy arrays
            * radii (list) of  radius of particle
            
        Kwargs:
            * radii_buffer (float) to multiply radii cut off
            * write_dr (boolean) 
            * del_drmatrix (boolean) 
            
        '''
        
        self.index = []
        self.list = []
        self.cnt = -1

        # Record rd
        if( write_dr ):
            
            dr_file = 'dr.csv'
            fout = open(dr_file,'wb')
            pair_writer = csv.writer(fout,delimiter=str(","))
            header = ['key_i','key_j','dr']
            #if( rank == 0 ):
            pair_writer.writerow(header)            
        
        # Create 2D list of lists of inter particle distances
        npos_i = positions
        npos_j = positions
        self.dr_matrix, self.dist_matrix  = lat.delta_npos(npos_i,npos_j)
        # Loop over all particles
        for key_i  in range(len(npos_i)):
            self.index.append(self.cnt + 1)
            radi_i = radii[key_i]
            for key_j  in range(len(npos_j)):
                radi_j = radii[key_j]
                if( key_i != key_j):
                    dr_cut = radi_i + radi_j
                    dr_cut = dr_cut*radii_buffer
                    dr = self.dist_matrix[key_i,key_j] 
                    if( dr <= dr_cut ):
                        self.cnt += 1
                        self.list.append(key_j)
                        if( write_dr  ):
                            row_i = [key_i,key_j,dr]
                            pair_writer.writerow(row_i)


        # Record rd
        if( write_dr ):
            fout.close()
                        
        # Add extra index positions for key+1 call made by final key 
        self.index.append(self.cnt + 1)
        if( del_drmatrix ):
            # Clear list from memory
            del self.dr_matrix
            del self.dist_matrix
        else:
            logger.debug(" Saving dr and dist matrix for groupset ")
        

