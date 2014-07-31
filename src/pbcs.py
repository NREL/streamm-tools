#! /usr/bin/env python
"""
Vector opperations for particles in a box 
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/07/2014
# travis.kemper@nrel.gov

import numpy as np


def sq_drij_c(r_i,r_j,latticevec):
    """
    Reutrn square magnitude of distance between to vectors  using cubic periodic boundry conditions 
    """
    
    # Find magnitude of dr_ij
    r_ij  = r_j - r_i
    
    r_x = r_ij[0] - latticevec[0][0] * round( r_ij[0]/  latticevec[0][0] )
    r_y = r_ij[1] - latticevec[1][1] * round( r_ij[1]/  latticevec[1][1] )
    r_z = r_ij[2] - latticevec[2][2] * round( r_ij[2]/  latticevec[2][2] )
    
    dr_pbc = np.array( [r_x,r_y,r_z] )
    
    sq_dr = np.dot( dr_pbc,dr_pbc)
    
    return sq_dr
    
    

def norm_r_ij(r_i,r_j):
    """
    Normailze difference between two vectors 
    """
    delta_ij = r_j - r_i
    return (delta_ij)/np.linalg.norm(delta_ij)

def getAngle(r_i,r_j):
    """
    Calcuate angle
      k - i - j 
      r_i = r_ik
      r_j = r_ij
      cos( \theta ) = ( a dot b ) / ( |a| |b| )
    """
    #

    r_i_norm = r_i/np.linalg.norm(r_i)
    r_j_norm = r_j/np.linalg.norm(r_j) 
    dot_ij = np.dot(r_i_norm,r_j_norm)
    
    if( dot_ij >= 1.0 ):
       ang_deg = 0.0
    elif(  dot_ij <= -1.0 ):
       ang_deg = 180.0
    else:    
        cos_ang = np.arccos(dot_ij )
        ang_deg = np.rad2deg( cos_ang )
    
    
    return ang_deg
