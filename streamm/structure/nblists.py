# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

'''
This module defines the classes relating to neighbor lists 
'''
              
class NBlist(object):
    """
    Class for neighbor list
    """
    def __init__(self):
        """
        Constructor
        """
        self.list = []
        self.index = []
        self.cnt = -1 

    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        del self.list
        del self.index
        del self.cnt
        
    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " NBlist of {} particle with {} connections".format(len(self.index),len(self.list))
    

    def getnbs(self,key_i):
        '''
        Return list of keys of neighbors of key_i        
        '''
        try:
            nbs_i = self.list[self.index[key_i]:self.index[key_i+1]]
        except:
            logger.warning(" Neighbor list not set ")
            nbs_i = []
            
        return nbs_i

    def calc_nnab(self,key_i):
        '''
        Calculate the number of neighbors for particle key_i
        '''
        try:
            nnab_i = self.index[key_i+1] - self.index[key_i]
        except:
            logger.warning(" Neighbor list not set ")
            nnab_i = 0
        #
        return nnab_i

    def build_nblist(self,particles,bonds):
        """
        Create neighbor list of bonded particles

        Arguments:
            particles (dict) of particles with int keys
            bonds (dict) of bonds with int keys of particles 

        """
        self.list = []
        self.index = []
        self.cnt = -1 
        
        # Create 2D list of lists for each particle 
        nd2D = [ [] for pkey_i  in particles.keys() ]
        # Fill each particle list with it's neighbors based on the bonds
        for bkey_i, bond_i  in bonds.iteritems():            
            nd2D[bond_i.pkey2].append( bond_i.pkey1 )
            nd2D[bond_i.pkey1].append( bond_i.pkey2 )
        # Loop over all particles and add it's neighbors to 1D list  (NBlist.list)          
        #   while tracking the index in the 1D in the index list  (NBlist.index)
        for pkey_i  in particles.keys():
            self.index.append(self.cnt + 1)
            for pkey_j in nd2D[pkey_i]:
                if( pkey_i != pkey_j):
                    self.cnt += 1
                    self.list.append(pkey_j)
                    
        # Add extra index positions for key+1 call made by final key 
        self.index.append(self.cnt + 1)
        # Clear 2D list from memory 
        del nd2D 

    def guess_nblist(self,lat,particles,positions,radii_key,radii_buffer=1.25):
        """
        Create neighbor list of particles based on distance and radius of each particle 
        
        Arguments:
            lat (Lattice) object 
            particles (dict) of particles with int keys
            positions (list) of  particle position numpy arrays
            radii_key (str) of key for particle radius in the particles dict
            radii_buffer (float) to multiply radii cut off 

        """
        self.list = []
        self.index = []
        self.cnt = -1 
        # Create 2D list of lists of inter particle distances
        npos_i = positions
        npos_j = positions
        dr_matrix, dist_matrix  = lat.delta_npos(npos_i,npos_j)
        # Loop over all particles
        for pkey_i,particle_i  in particles.iteritems():
            radii_i = particle_i.properties[radii_key]
            self.index.append(self.cnt + 1)
            for pkey_j,particle_j in particles.iteritems():
                if( pkey_i != pkey_j):
                    radii_j = particle_j.properties[radii_key]
                    dr_cut = radii_i + radii_j
                    dr_cut = dr_cut*radii_buffer
                    print "Particles  i_%d - j_%d dr %f cut %f "%(pkey_i,pkey_j,dist_matrix[pkey_i,pkey_j],dr_cut)
                    if( dist_matrix[pkey_i,pkey_j] <= dr_cut ):
                        self.cnt += 1
                        self.list.append(pkey_j)
                    
        # Add extra index positions for key+1 call made by final key 
        self.index.append(self.cnt + 1)
        # Clear list from memory 
        del dr_matrix
        del dist_matrix


    def radii_nblist(self,lat,positions,radii,radii_buffer=1.25,write_dr=True,del_drmatrix=False):
        """
        Create neighbor list of particles based on distance and radius of each particle 
        
        Arguments:
            lat (Lattice) object 
            positions (list) of  particle position numpy arrays
            radii (list) of  radius of particle
            radii_buffer (float) to multiply radii cut off 

        """
        self.index = []
        self.list = []
        self.cnt = -1

        # Record rd
        if( write_dr ):
            dr_file = 'dr.csv'
            fout = open(dr_file,'wb')
            pair_writer = csv.writer(fout,delimiter=',')
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
        
