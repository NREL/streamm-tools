# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"
                

''' Functions to replicate structure containers '''                
           

import numpy as np
import copy
from datetime import datetime
import random
import math 


        
import streamm.mpi.mpiBase as mpiBase 

import logging
logger = logging.getLogger(__name__)


class Replication(object):
    '''
    Object to record the replication of structure 
    '''
    def __init__(self,name_i,name_j,name_ij,method,n):

        self.name_i = name_i
        self.name_j = name_j
        self.name_ij = name_ij
        self.method = method
        self.n = n
    
    def __del__(self):
        """
        'Magic' method for deleting contents of container
        """
        del self.name_ij
        del self.name_i
        del self.name_j
        del self.method
        del self.n
        
    def __str__(self):
        """
        'Magic' method for printng contents of container 
        """
        return " %s +  %s x %d ( %s ) -> %s "%(self.name_i,self.name_j,self.n,self.method,self.name_ij)

       
def add_struc(strucC_i,other,n_i,seed,p=None,tag="blank"):
    """
    Add structure other to strucC_i n times via random placement
    """
    import random
    #
    # MPI setup
    #
    rank = 0
    size = 0
    if( p != None ):
        rank = p.getRank()
        size = p.getCommSize()
    #
    # Initialize random
    
    if( rank == 0 ):            
        logging.info('Starting add_struc %s to %s with seed %d '%(other.tag,strucC_i.tag,seed))
        logging.info(' into a simulation cell of size {} '.format(strucC_i.lat))
        
    random.seed(seed)
    #
    pos_cut = 10.0 # minimum distance between particles of added structure and current structure
    ang_acc = 1000  # number of digets in random angle
    max_sys = 3
    max_mol_place = 10 
    exlat_frac = 0.10 # Fraction to expand the lattice by if max_sys is exceeded 
    #
    # Initilize
    exlat_cnt = 0                   # Number of lattice expantions 
    reset_cnt = 0                   # Number of times the new structure  has been reset
    struc_add_cnt = 0                # Total number of structures added to the new structure
    strucC_new = copy.deepcopy(strucC_i) # Create copy to reset to if addition of other fails
    strucC_new.tag = strucC_i.tag 
    strucC_new.lat = strucC_i.lat 
    # Set lattice to be recoppied over the strucC_new during initial add
    other.lat = strucC_i.lat         
    other.calc_mass()
    n_mol_o = strucC_i.n_molecules()
    if( n_mol_o > 0 ):
        # Add one here since we count from zero when adding new structures  
        n_mol_o +=  1
    # 
    # Create a list of atomic indices for each processor 
    # particle_keys = strucC_i.particles.keys()  
    # particle_keys_p  = p.splitListOnProcs(particle_keys)
    #
    #
    # Start adding molecules to the system
    #   
    add_strucC = True
    while ( add_strucC ):
        #
        # Initialize 
        #
        add_strucC = True
        poxpass  = False 
        placement_cnt = 0
        # 
        # For each structure add to 
        # while ( overlap_sum != 0  ):
        while ( not poxpass ):
            placement_cnt += 1
            if( rank == 0 ):
                logging.debug('-Placement %d %d '%(placement_cnt,rank))

            rot_angle_i_o = 0.0 
            rot_angle_j_o = 0.0 
            r_random_o  = np.zeros(strucC_new.lat.n_dim) 

            if ( rank == 0 ):
                #
                #  Get random rotation angles from single processor  
                #
                rot_angle_i_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                rot_angle_j_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                #
                #  Get random translation from single processor 
                r_random_o = strucC_new.lat.random_pos()
            if( size > 0 ):
                # Broadcast random rotation angles and translations to all processors 
                p.barrier() # Barrier for MPI_COMM_WORLD
                rot_angle_i = p.bcast(rot_angle_i_o)
                rot_angle_j = p.bcast(rot_angle_j_o)
                r_random = p.bcast(r_random_o)
                p.barrier() # Barrier for MPI_COMM_WORLD
            else:
                rot_angle_i = rot_angle_i_o
                rot_angle_j = rot_angle_j_o
                r_random = r_random_o

            if( rank == 0 ):
                logging.debug('Angle i %f angle j %f on %d '%(rot_angle_i,rot_angle_j,rank))
            #
            # Get coordinates of randomly rotated and shifted
            #
            other.calc_center_mass()
            other.shift_pos(-1.0*other.center_mass)  # Place center of mass at origin
            other.rotate_xy(rot_angle_i)
            other.rotate_xz(rot_angle_j)
            other.shift_pos(r_random)  # Place center of mass at random position in lattice 
            #
            # Calculate inter particle seperation between particles in the structure and structure to be added
            #
            poxpass = True
            if( strucC_new.n_particles > 0 ):
                npos_i = strucC_new.positions
                npos_j = other.positions
                poxpass = strucC_new.lat.proximitycheck(npos_i,npos_j,pos_cut,p=p)
                
            # If no overlap detected add molecule to the system
            if( poxpass ):
                if ( rank == 0 ):
                    # logging.debug('No overlap found adding structure %d'%(struc_add_cnt))
                    print 'No overlap found adding structure %s'%(struc_add_cnt)
                # Create copy of structure to add to preserve original state 
                structoadd= copy.deepcopy(other)
                # Update mol number
                for pkey_i, particle_i  in structoadd.particles.iteritems():
                    particle_i.mol = n_mol_o + struc_add_cnt
                struc_add_cnt += 1
                strucC_new += structoadd
                
            if( placement_cnt >= max_mol_place ):

                # If attempts to place molecule into the system exceed max set by max_mol_place
                #   reset system and star over 
                if(  rank == 0  ):
                    msg = 'Max placments %d exceeded resetting to original system '%(max_mol_place)
                    logging.info(msg)
                    print msg
                reset_cnt += 1                 # Number of times the new structure  has been reset
                struc_add_cnt = 0                    # Number of structures add to the new structure 
                placement_cnt = 0
                matrix_i = copy.deepcopy(strucC_new.lat.matrix) # Save state of lattice
                # Reset new structure to original state 
                strucC_new = copy.deepcopy(strucC_i) # Create copy to add new structures to 
                strucC_new.lat.matrix = matrix_i  # Set lattice to saved state 

                    
            if( reset_cnt >= max_sys  ):
                if(  rank == 0  ):
                    logging.debug('Max resets %d exceeded expanding box size by %f '%(max_sys,exlat_frac))

                exlat_cnt += 1
                # If the new structure has been reset over max_sys times expand the box size by lc_expand
                matrix_i = copy.deepcopy(strucC_new.lat.matrix) # Save state of lattice
                # Reset new structure to original state 
                strucC_new = copy.deepcopy(strucC_i) # Create copy to add new structures to 
                strucC_new.lat.matrix  = matrix_i  # Set lattice to saved state
                strucC_new.lat.expand_matrix(exlat_frac)

                reset_cnt = 0
                struc_add_cnt = 0                    # Number of structures add to the new structure 
                placement_cnt = 0



        if(  rank == 0  ):

                log_line =  "  %d  added after "%(struc_add_cnt)
                log_line += "%d  expantions "%(exlat_cnt)
                log_line += "%d resets "%(reset_cnt)
                log_line += "%d placements "%(placement_cnt)
                logging.info(log_line)


        if( struc_add_cnt ==  n_i  ):
            if(  rank == 0  ):
                logging.info(" All structures added with %d particles "%(strucC_new.n_particles))
            # If all the molecule have been added exit while loop and print system 
            add_strucC = False


    # Record replication
    strucC_new.tag = tag
    replication_i = Replication(strucC_i.tag,other.tag,strucC_new.tag,"random",n_i)
    strucC_new.replications.append(replication_i)        

    if( rank == 0 ):
        logging.info("Finished %s "%(datetime.now()))
       
    return strucC_new


def add_struc_grid(strucC_i,other,n_i,p=None,tag="blank",calc_overlap = True ):
    """
    Add structure other to strucC_i n times on a grid
    
    Args:
        n_i  (int) number of times to replicate 
        
    """
    #
    # Initialize 
    rank = 0
    size = 0
    pos_cut = 10.0 
    ang_acc = 1000  # number of digets in random angle
    exlat_cnt = 0                   # Number of lattice expantions 
    exlat_indx = 0 # Initial lattice vector index to expand 
    exlat_frac = 0.10 # Fraction to expand the lattice by if max_sys is exceeded 
    #
    # Initilize
    strucC_new = copy.deepcopy(strucC_i) # Create copy to reset to if addition of other fails 
    strucC_new.lat = strucC_i.lat 
    # Set lattice to be recoppied over the strucC_new during initial add
    other.lat = strucC_new.lat
    n_mol_o = strucC_i.n_molecules()
    if( n_mol_o > 0 ):
        # Add one here since we count from zero when adding new structures  
        n_mol_o +=  1
        
    # 
    # Create a list of atomic indices for each processor 
    # particle_keys = strucC_i.particles.keys()  
    # particle_keys_p  = p.splitListOnProcs(particle_keys)
    #
    # Guess size of volume need for each added structure 
    other.calc_mass()
    strucC_new.calc_volume()
    n_vol = strucC_new.volume/float(n_i)  # Number of other strucC per volume
    l_n = n_vol**(1.0/3.0)
    n_lat = [] #np.zeros(strucC_new.lat.n_dim)
    for d in range(strucC_new.lat.n_dim):
        n_lat.append(int(math.ceil( strucC_new.lat.lengths[d]/l_n)))
    if(  rank == 0 ):
        log_line = " Adding %d into volume %f \n"%(n_i,strucC_new.volume)
        log_line += "  %d along v_1 \n"%(n_lat[0])
        log_line += "  %d along v_2 \n"%(n_lat[1])
        log_line += "  %d along v_3 \n"%(n_lat[2])
        log_line += "  %d lattice positions \n"%(n_lat[0]*n_lat[1]*n_lat[2])
        logger.debug(log_line)
    #
    # Start adding molecules to the system
    add_strucC = True
    while ( add_strucC ):
        add_strucC = True
        struc_add_cnt = 0                    # Number of structures add to the new structure
        for indx_i in range(n_lat[0]):
            if( not add_strucC): break 
            for indx_j in range(n_lat[1]):
                if( not add_strucC): break 
                for indx_k in range(n_lat[2]):
                    # Calculate fractional coord
                    f_i = float(indx_i)/float(n_lat[0])
                    f_j = float(indx_j)/float(n_lat[1])
                    f_k = float(indx_k)/float(n_lat[2])
                    frac_o  = np.array([f_i,f_j,f_k])
                    pos_o = strucC_new.lat.fractoreal(frac_o)
                    #
                    if ( rank == 0 ):
                        #
                        #  Get random rotation angles from single processor  
                        #
                        rot_angle_i_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                        rot_angle_j_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                    if( size > 0 ):
                        # Broadcast random rotation angles and translations to all processors 
                        p.barrier() # Barrier for MPI_COMM_WORLD
                        rot_angle_i = p.bcast(rot_angle_i_o)
                        rot_angle_j = p.bcast(rot_angle_j_o)
                        r_random = p.bcast(pos_o)
                        p.barrier() # Barrier for MPI_COMM_WORLD
                    else:
                        rot_angle_i = rot_angle_i_o
                        rot_angle_j = rot_angle_j_o
                        r_random = pos_o
                    #
                    # Get coordinates of randomly rotated and shifted
                    other.calc_center_mass()
                    other.shift_pos(-1.0*other.center_mass)  # Place center of mass at origin
                    other.rotate_xy(rot_angle_i)
                    other.rotate_xz(rot_angle_j)
                    other.shift_pos(r_random)  # Place center of mass at random position in lattice 
                    #
                    # Calculate inter particle seperation between particles in the structure and structure to be added
                    overlap = 0                        
                    if( strucC_new.n_particles > 0 and calc_overlap ):
                        npos_i = other.positions
                        npos_j = strucC_new.positions
                        npos_ij,nd_ij = strucC_new.lat.delta_npos(npos_i,npos_j)
                        # If the particles are further than the cut off
                        for key_i in other.particles.keys():
                            for key_j in strucC_new.particles.keys():
                                if( nd_ij[key_i][key_j] <= pos_cut ):
                                    overlap = 1
                    if( size > 0 ):
                        #
                        # Reduce sum the overlap variable from all the processors
                        #   if it is zero everywhere there was no overlap detected 
                        overlap_sum = p.allReduceSum(overlap)
                        p.barrier() # Barrier for MPI_COMM_WORLD
                    else:
                        overlap_sum = overlap
                    #
                    # If no overlap detected add molecule to the system
                    if( overlap_sum == 0 ):
                        # Create copy of structure to add to preserve original state 
                        structoadd= copy.deepcopy(other)
                        # Update mol number
                        for pkey_i, particle_i  in structoadd.particles.iteritems():
                            particle_i.mol = n_mol_o + struc_add_cnt
                        struc_add_cnt += 1
                        strucC_new += structoadd
                        logger.info( "Molecule %d/%d added "%(struc_add_cnt,n_i))
                    if( struc_add_cnt == n_i ):
                        add_strucC = False

                        # Record replication
                        strucC_new.tag = tag 
                        replication_i = Replication(strucC_i.tag,other.tag,strucC_new.tag,"grid",n_i)
                        strucC_new.replications.append(replication_i)        
                        
                        return strucC_new
        #
        # If attempts to place structure into the new structure failed
        exlat_cnt += 1
        #
        # If the new structure has been reset over max_sys times expand the box size by lc_expand
        matrix_i = copy.deepcopy(strucC_new.lat.matrix) # Save state of lattice
        # Reset new structure to original state 
        strucC_new = copy.deepcopy(strucC_i) # Create copy to add new structures to 
        strucC_new.lat.matrix = matrix_i  # Set lattice to saved state
        strucC_new.lat.expand_matrix(exlat_frac)
        #
        # Add one row to a single dimension of the lattice grid 
        n_lat[exlat_indx] += 1
        exlat_indx += 1
        if( exlat_indx > 2 ):
                exlat_indx = 0
        #
        log_line =  "  %d  added after "%(struc_add_cnt)
        log_line += "%d  exlat_indx "%(exlat_indx)
        log_line += "  n_lat 0 %f "%(n_lat[0])
        log_line += "  n_lat 1 %f "%(n_lat[1])
        log_line += "  n_lat 2 %f "%(n_lat[2])
        log_line += "%d  expantions "%(exlat_cnt)
        logger.info(log_line)
                    

    return None
     