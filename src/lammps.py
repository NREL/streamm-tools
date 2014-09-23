#! /usr/bin/env python
"""
subroutines for lammps file manipulation
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov

from particles  import Particle
from bonds      import Bond,     BondContainer
from angles     import Angle,    AngleContainer
from dihedrals  import Dihedral, DihedralContainer

def read_lmpdata( strucC , data_file):
    """
    Read Lammps data file

    Arguments:
        # strucC     (StructureContainer)
        # parmC      (ParameterContainer)
        data_file  (str) data file

    ReturnS:
        strucC     (StructureContainer)
        parmC      (ParameterContainer)
    
    """
    
    
    import file_io
    import sys, numpy
    
    debug = 0
    verbose = 0

    if( file_io.file_exists(data_file )):
        if( debug ): print "Reading in ",data_file
    else:
        print " Specified file ",data_file," does not exisit "
        sys.exit("Invalid file ")
        

    F = open(data_file , 'r' )
    lines = F.readlines()
    F.close()
    
    #
    # Read in data header with number of parameters 
    #
    for line in lines:
        col = line.split()
        if ( len(col) >=2 ):
            # Read in number of each topolgical component  
            if( col[1] == "atoms" ):
                n_atoms = int( col[0] )
            elif( col[1] == "bonds" ):
                n_bonds = int( col[0] )
            elif( col[1] == "angles" ):
                n_angles = int( col[0] )
            elif( col[1] == "dihedrals" ):
                n_dihedrals = int( col[0] )
            elif( col[1] == "impropers" ):
                n_impropers = int( col[0] )

        if ( len(col) >= 3 ):
            # Read in number of each parameter type 
            if( col[1] == "atom" and   col[2] == "types" ):
                n_atypes = int( col[0] )
            elif( col[1] == "bond" and   col[2] == "types"  ):
                n_btypes = int( col[0] )
            elif( col[1] == "angle" and   col[2] == "types"  ):
                n_angtypes = int( col[0] )
            elif( col[1] == "dihedral" and   col[2] == "types"  ):
                n_dtypes = int( col[0] )
            elif( col[1] == "improper" and   col[2] == "types"  ):
                n_imptypes = int( col[0] )
                
        # Read in box size    
        if ( len(col) >= 4 ):
            if( col[2]  == "xlo"  and col[3]  == "xhi"  ):
                strucC.latvec[0][0] = float( col[1] ) - float( col[0] ) 
            if( col[2]  == "ylo"  and col[3]  == "yhi"  ):
                strucC.latvec[1][1] = float( col[1] ) - float( col[0] ) 
            if( col[2]  == "zlo"  and col[3]  == "zhi"  ):
                strucC.latvec[2][2] = float( col[1] ) - float( col[0] )

    # Prind debug
    if( verbose ):
        print " atoms ",n_atoms
        print " n_bonds ",n_bonds
        print " n_angles ",n_angles
        print " n_dihedrals ",n_dihedrals
        print ""
        print "n_atom_types",n_atypes
        print "n_bond_types",n_btypes
        print "n_angle_types",n_angtypes
        print "n_dihedral_types",n_dtypes

    # Check to see if a previous read has occured
    pt_overwrite = False
    if( len(strucC.ptclC) > 0 ):
        pt_overwrite = True
    # Check of conistent number of atoms
    if( pt_overwrite ):
        if(  len(strucC.ptclC) + 1 != n_atoms):
            sys.exit(" Inconsistent number of atoms " )

    bonds_overwrite = False
    if( len(strucC.bondC) > 0 ):
        bonds_overwrite = True
    angles_overwrite = False
    if( len(strucC.angleC) > 0 ):
        angles_overwrite = True
    dih_overwrite = False
    if( len(strucC.dihC) > 0 ):
        dih_overwrite = True
            
    #
    # Intialize
    #   - read in boolean to off
    #
    read_Masses = 0
    read_Pair = 0
    read_Bond_coeff = 0
    read_Angle_coeff = 0
    read_Dihedral_coeff = 0
    read_Improper_coeff = 0
    
    read_Atoms = 0
    read_Bonds = 0
    read_Angles = 0
    read_Dihedrals = 0
    read_Impropers = 0
    
    #   - lists as indecise can be out of order in data file
    ATYPE_REF = n_atypes*[""]
    ATYPE_MASS = numpy.zeros(n_atypes)
    ATYPE_EP = numpy.zeros(n_atypes)
    ATYPE_SIG = numpy.zeros(n_atypes)
    
    BTYPE_REF = n_btypes*[2*[""]]
    BONDTYPE_R0 = numpy.zeros(n_btypes)
    BONDTYPE_K = numpy.zeros(n_btypes)
    
    ANGTYPE_REF = n_angtypes*[3*[""]]
    ANGLETYPE_R0 = numpy.zeros(n_angtypes)
    ANGLETYPE_K = numpy.zeros(n_angtypes)
    
    DTYPE_REF = n_dtypes*[4*[""]]
    DIHTYPE_C = numpy.zeros((n_dtypes,4))
    DIHTYPE_F = numpy.zeros(n_dtypes)
    DIHTYPE_K = numpy.zeros(n_dtypes)
    DIHTYPE_PN = numpy.zeros(n_dtypes)
    DIHTYPE_PHASE = numpy.zeros(n_dtypes)
    
    MOLNUMB = n_atoms*[0]
    ATYPE_IND  = n_atoms*[0]
    CHARGES  = numpy.zeros(n_atoms)
    R = n_atoms*[numpy.zeros(3)]
    ATYPE  = n_atoms*[""]
    
    BONDS = n_bonds*[[0,0]]
    BTYPE_IND = n_bonds*[0]
    
    ANGLES = n_angles*[[0,0,0]]
    ANGTYPE_IND = n_angles*[0]
    
    DIH = n_dihedrals*[[0,0,0,0]]
    DTYPE_IND = n_dihedrals*[0]
    
    
    #
    # Read in data parameters 
    #
    
    for line in lines:
        col = line.split()
        if ( len(col) >= 1  ):
            if( col[0] == "Masses" ):
                read_Masses = 1
                cnt_Masses = 0
                
            if( col[0] == "Atoms" ):
                read_Atoms = 1
                cnt_Atoms = 0
                
            if( col[0] == "Bonds" ):
                read_Bonds = 1
                cnt_Bonds = 0
                
            if( col[0] == "Angles" ):
                read_Angles = 1
                cnt_Angles  = 0
                
            if( col[0] == "Dihedrals" ):
                read_Dihedrals = 1
                cnt_Dihedrals = 0
                
            if( col[0] == "Impropers" ):
                read_Impropers = 1
                cnt_Impropers = 0
                
        if ( len(col) >= 2 ):
            if( col[0] == "Pair" and col[1] == "Coeffs" ):
                read_Pair = 1
                cnt_Pair = 0
            if( col[0] == "Bond" and col[1] == "Coeffs" ):
                read_Bond_coeff = 1
                cnt_Bond_coeff = 0
            if( col[0] == "Angle" and col[1] == "Coeffs" ):
                read_Angle_coeff  = 1
                cnt_Angle_coeff = 0
            if( col[0] == "Dihedral" and col[1] == "Coeffs" ):
                read_Dihedral_coeff  = 1
                cnt_Dihedral_coeff  = 0
                
            if( col[0] == "Improper" and col[1] == "Coeffs" ):
                read_Improper_coeff  = 1
                cnt_Improper_coeff  = 0
                
        if( read_Masses and  len(col) >= 2 ):
            cnt_Masses += 1
            ind = int(col[0]) - 1
            ATYPE_MASS[ind] = float(col[1])
            if(cnt_Masses ==  n_atypes ):
                read_Masses = 0
            if( len(col) >= 4 ):
                ATYPE_REF[ind] = col[3]
            else:
                ATYPE_REF[ind] = "??"

        if( read_Pair and  len(col) >= 3 ):
            cnt_Pair += 1
            ind = int(col[0]) - 1
            ATYPE_EP[ind] = float(col[1])
            ATYPE_SIG[ind] = float(col[2])            
            if(cnt_Pair >=  n_atypes ):
                read_Pair = 0

            
        
        if( read_Bond_coeff and  len(col) >= 3 ):
            cnt_Bond_coeff += 1
            #AT_i = int(col[0])
            #AT_j = int(col[1])
            b_ind = int( col[0]) - 1
            if( b_ind > n_btypes ):
                print sys.exit(" Error in data file index of bond parameter exceeds number of bond parameters specified with bond types ")
                
            BTYPE_REF[b_ind][0] = "??"
            BTYPE_REF[b_ind][1] = "??"
            BONDTYPE_K[b_ind] = float(col[1])
            BONDTYPE_R0[b_ind] = float(col[2])
            
            if( cnt_Bond_coeff >=  n_btypes ):
                read_Bond_coeff = 0
                
        
        if( read_Angle_coeff and  len(col) >= 3 ):
            cnt_Angle_coeff += 1
            #AT_i = int(col[0])
            #AT_j = int(col[1])
            a_ind = int( col[0]) - 1
            if( a_ind > n_angtypes ):
                print sys.exit(" Error in data file index of angle parameter exceeds number of angle parameters specified with angle types ")
                
            ANGTYPE_REF[a_ind][0] = "??"
            ANGTYPE_REF[a_ind][1] = "??"
            ANGTYPE_REF[a_ind][2] = "??"
            ANGLETYPE_K[a_ind] = float(col[1])
            ANGLETYPE_R0[a_ind] = float(col[2])
            
            if( cnt_Angle_coeff >=  n_angtypes ):
                read_Angle_coeff = 0
                
        
        
        if( read_Dihedral_coeff and  len(col) >= 3 ):
            cnt_Dihedral_coeff += 1
            #AT_i = int(col[0])
            #AT_j = int(col[1])
            d_ind = int( col[0]) - 1
            
            if( debug): print " reading dih type ",d_ind," cnt ",cnt_Dihedral_coeff," of ",n_dtypes
            
            if( d_ind > n_dtypes ):
                print sys.exit(" Error in data file index of dihedral parameter exceeds number of dihedral parameters specified with dihedral types ")
                
            DTYPE_REF[d_ind][0] = "??"
            DTYPE_REF[d_ind][1] = "??"
            DTYPE_REF[d_ind][2] = "??"
            DTYPE_REF[d_ind][3] = "??"
            
            # Assume OPLS dihedral type
            DIHTYPE_F[d_ind] = 3
            DIHTYPE_C[d_ind][0] = float(col[1])
            DIHTYPE_C[d_ind][1] = float(col[2])
            DIHTYPE_C[d_ind][2] = float(col[3])
            DIHTYPE_C[d_ind][3] = float(col[4])
            
            if( cnt_Dihedral_coeff >=  n_dtypes ):
                read_Dihedral_coeff = 0
                
        if(read_Atoms and len(col) >= 7 ):
            cnt_Atoms += 1
            ind = int( col[0]) - 1
            if( ind > n_atoms ):
                print sys.exit(" Error in data file index of atoms exceeds number of atoms specified with atoms ")
                
            chain_i = int(col[1]) 
            lmptype_i = int(col[2]) - 1
            q_i = float(col[3])
            m_i = ATYPE_MASS[lmptype_i]
            r_i =  [ float(col[4]),float(col[5]),float(col[6])] 
            type_i = "??"

            tagsD = {"chain":chain_i}
            if( pt_overwrite ):
                pt_i = strucC.ptclC[cnt_Atoms]
                pt_i.setTagsDict(tagsD)
                pt_i.position = r_i
                pt_i.charge = q_i
                pt_i.mass = m_i
            else:
                pt_i = Particle( r_i,type_i,q_i,m_i)
                pt_i.setTagsDict(tagsD)
                strucC.ptclC.put(pt_i)

            if( cnt_Atoms >=  n_atoms ):
                read_Atoms = 0

        if(read_Bonds and len(col) >= 4 ):
            cnt_Bonds += 1
            ind = int( col[0]) - 1

            if( ind > n_bonds ):
                print sys.exit(" Error in data file index of bonds exceeds number of bonds specified with bonds ")
                
            BTYPE_IND[ind] = int(col[1] ) - 1
            BONDS[ind] = [int(col[2]) - 1 , int(col[3]) - 1 ]
            i_o = int(col[2])
            j_o = int(col[3])

            if( bonds_overwrite ):
                bondObj = strucC.bondC[cnt_Bonds]
                bondObj.pgid1 = i_o
                bondObj.pgid2 = j_o
            else:
                b_i = Bond( i_o, j_o )            
                strucC.bondC.put(b_i)
                
            
            if( cnt_Bonds >=  n_bonds ):
                read_Bonds = 0

            
                
            
        if(read_Angles and len(col) >= 5 ):
            cnt_Angles += 1
            ind = int( col[0]) - 1
            ANGTYPE_IND[ind] = int(col[1] ) - 1
            ANGLES[ind] = [int(col[2]) - 1, int(col[3]) - 1, int(col[4]) - 1 ]
            
            k_o = int(col[2])
            i_o = int(col[3])
            j_o = int(col[4])
            
            if( cnt_Angles >=  n_angles ):
                read_Angles = 0

            if( angles_overwrite ):
                angleObj = strucC.angleC[cnt_Angles]
                angleObj.pgid1 = k_o
                angleObj.pgid2 = i_o
                angleObj.pgid3 = j_o
            else:
                angle_i = Angle( k_o,i_o, j_o )            
                strucC.angleC.put(angle_i)
                
            
        if(read_Dihedrals and len(col) >= 6 ):
            cnt_Dihedrals += 1
            ind = int( col[0]) - 1
            
            DTYPE_IND[ind] = int(col[1] ) - 1
            DIH[ind] = [int(col[2]) - 1,int(col[3]) - 1, int(col[4]) - 1,int(col[5]) - 1]

            k_o = int(col[2])
            i_o = int(col[3])
            j_o = int(col[4])
            l_o = int(col[5])
            
            if( dih_overwrite ):
                dObj = strucC.dihC[cnt_Dihedrals]
                dObj.pgid1 = k_o
                dObj.pgid2 = i_o
                dObj.pgid3 = j_o
                dObj.pgid4 = l_o
            else:
                dih_i = Dihedral( k_o,i_o, j_o,l_o )            
                strucC.dihC.put(dih_i)


            if( read_Dihedrals >=  n_dihedrals ):
                read_Dihedrals = 0
                
        #    
        #if(read_Impropers and len(col) >= 2 ):
        #    cnt_Bonds += 1
        #    ind = int( col[0]) - 1
        #    BTYPE_IND[ind] = int(col[1] ) - 1
        #    BONDS[ind][0] = int(col[2])
        #    if( cnt_Bonds >=  n_atoms ):
        #        read_Bonds = 0
        #        
            
            
        
    if( debug):
        for ind in range(len(ATYPE_MASS)):
            print ind+1,ATYPE_MASS[ind]
            
        for ind in range(len(ATYPE_EP)):
            print ind+1,ATYPE_EP[ind],ATYPE_SIG[ind]
        
        for ind in range(n_btypes):
            print ind+1,BONDTYPE_R0[ind],BONDTYPE_K[ind]
                
        for ind in range(n_angtypes):
            print ind+1,ANGLETYPE_R0[ind],ANGLETYPE_K[ind]
        
        for ind in range(n_dtypes):
            print ind+1,DIHTYPE_C[ind]
        
    debug =0
    
    if( debug):
        for ind in range(len(BONDS)):
            print ind+1,BONDS[ind]
            
    if(debug):
        sys.exit("debug 1 ")
    #
    #      
    return (strucC)

