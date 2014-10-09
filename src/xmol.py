"""
Read and write xmol files 
"""

# Dr. Travis Kemper
# NREL
# Initial Date 6/27/2014
# travis.kemper@nrel.gov


from particles import Particle
from particles import ParticleContainer


def read(xmol_file,format):
    """
    Read in xmol file
    
    Args:
        xmol_file (str) xmol file name

    Returns
        struc_array (list) of ParticleContainer's for each frame
        
    """

    verbose = False
    
    # Initialize structure array for containing multiple frames
    ptclC_array = []
    
    # Initialize number of particles to zero for frame checking
    NP = 0
    # Set charge and mass to zero since they are not in xmol file
    q_i = 0.0
    m_i = 0.0 
    # Initialize line count and structure count 
    line_cnt = 0
    ptclC_cnt = 0 
    # Read in file line by line 
    with open(xmol_file) as f:
        for line in f:
            line_cnt += 1
            col = line.split()
            if( line_cnt == 1  ):
                # Read first line to record number of particles NP
                NP = int(col[0])
                ptclC_i =  ParticleContainer()
                ptclC_cnt += 1 
            elif( line_cnt > 2 and len(col) >= 4 ):
                # Read lines and add particles to structure
                if( format == "atomic_symb" ):
                    atomic_symb = col[0]
                elif( format == "lammps" ):
                    part_type = str(col[0])
                    atomic_symb = "C"
                else:
                    atomic_symb = "C"
                    
                r_i = [ col[1],col[2],col[3] ]
                part_i = Particle( r_i,atomic_symb,q_i,m_i )
                ptclC_i.put(part_i)
                del part_i
                
            if( line_cnt > 1  and line_cnt > ptclC_cnt*(NP + 2) ):
                # xmol file contains multiple frames
                #   store in structure array 
                ptclC_array.append(ptclC_i)
                #   reinitialize struc_i ParticleContainer
                del ptclC_i
                ptclC_i =  ParticleContainer()
                ptclC_cnt += 1 

    # Append last structure
    ptclC_array.append(ptclC_i)
    del ptclC_i



    for pid,pt_i in  struc_o.ptclC:
        print pid,pt_i.type, pt_i.postion
    sys.exit("debug xmol read in 1")
    
    return(ptclC_array)
            

def write(ptclC, xmol_file,comment,append):

    """
    Write a structure  to an xmol file

    Args:
        xmol_file    (str) xmol file name
        comment  (str) for comment line 
        append  (boolean) to append or create a new file 
    Reutrns
        null
    """
    # Open xmol file 
    if(append):
        F = open(xmol_file,"a")
    else:
        F = open(xmol_file,"w")

    # Loop over structures
    NP = len( ptclC )
    F.write(" %d \n" % NP )
    F.write(" %s \n"%comment)
    for pid, ptclObj  in ptclC:
        r_i = ptclObj.position
        atomic_symb = ptclObj.tagsDict['symbol']
        F.write( " %5s %16.8f %16.8f %16.8f \n"  % (atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]) ) )   
    F.close()


