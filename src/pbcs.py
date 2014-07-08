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
    
    r_ij = r_j - r_i
    
    r_x = r_ij[0] - latticevec[0][0] * round( r_ij[0]/  latticevec[0][0] )
    r_y = r_ij[1] - latticevec[1][1] * round( r_ij[1]/  latticevec[1][1] )
    r_z = r_ij[2] - latticevec[2][2] * round( r_ij[2]/  latticevec[2][2] )
    
    dr_pbc = np.array( [r_x,r_y,r_z] )
    
    sq_dr = np.dot( dr_pbc,dr_pbc)
    
    return sq_dr
    
    
