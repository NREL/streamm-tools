"""
Read and write xmol files 
"""

# Dr. Travis Kemper
# NREL
# Initial Date 6/27/2014
# travis.kemper@nrel.gov


from particles import Particle
from particles import ParticleContainer


def read(xmol_file):
    """
    Read in xmol file
    
    Args:
        xmol_file (str) xmol file name

    Returns
        struc_array (list) of ParticleContainer's for each frame
        
    """
    # Initialize structure array for containing multiple frames
    struc_array = []
    
    # Initialize number of particles to zero for frame checking
    NP = 0
    # Set charge and mass to zero since they are not in xmol file
    q_i = 0.0
    m_i = 0.0 
    # Initialize line count and structure count 
    line_cnt = 0
    struc_cnt = 0 
    # Read in file line by line 
    with open(xmol_file) as f:
        for line in f:
            line_cnt += 1
            col = line.split()
            if( line_cnt == 1  ):
                # Read first line to record number of particles NP
                NP = int(col[0])
                struc_i =  ParticleContainer()
                struc_cnt += 1 
            elif( line_cnt > 2 and len(col) >= 4 ):
                # Read lines and add particles to structure 
                atomic_symb = col[0]
                r_i = [ col[1],col[2],col[3] ]
                part_i = Particle( r_i,atomic_symb,q_i,m_i )
                struc_i.put(part_i)
                del part_i
                
            if( line_cnt > 1  and line_cnt > struc_cnt*(NP + 2) ):
                # xmol file contains multiple frames
                #   store in structure array 
                struc_array.append(struc_i)
                #   reinitialize struc_i ParticleContainer
                del struc_i
                struc_i =  ParticleContainer()
                struc_cnt += 1 

    # Append last structure
    struc_array.append(struc_i)
    del struc_i
    
    return(struc_array)
            
    
def write(struc_array,xmol_file):
    """
    Write a structure array to an xmol file

    Args:
        struc_array  (list) containing structures to write
        xmol_file    (str) xmol file name
    Reutrns
        null
    """
    # Initialize frame cnt
    frame_cnt = 0 
    # Open xmol file 
    F = open(xmol_file,"w")
    # Loop over structures
    for struc_i in struc_array:
        frame_cnt += 1 
        NP = len( struc_i )
        F.write(" %d \n" % NP )
        F.write(" Frame   %d \n"%frame_cnt)
        for pid, ptclObj  in struc_i:
            r_i = ptclObj.position
            atomic_symb = ptclObj.type
            F.write( "%5s %16.8f %16.8f %16.8f \n"  % (atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]) ) )
    F.close()

    
        
