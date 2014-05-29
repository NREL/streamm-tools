#! /usr/bin/env python
"""
I/O of lammps files 
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov


def lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH):
    """
    Find bonded types 
    """
    import sys
    
    debug  = 0 
    #
    # for lammmps files
    #
    ATYPE_IND = []  
    ATYPE_REF = []
    ATYPE_MASS = []
    BTYPE_IND = []
    BTYPE_REF = []
    ANGTYPE_IND = []
    ANGTYPE_REF = []
    DTYPE_IND = []
    DTYPE_REF = []
    IMPYPE_IND = []
    IMPTYPE_REF = []
    #
    # Count atom types
    #
    n_atoms = len(ATYPE)
    TYPE_CNT = -1 # so array is 0-N
    debug = 0
    for i in range(n_atoms):
        new_type = 1
        for ind in range( len(ATYPE_REF) ):
            if( ATYPE[i] == ATYPE_REF[ind] ):
                new_type = 0
                
                ATYPE_IND.append( ind  )
                #if(debug): print ' maches previous ',i,ATYPE[i],ATYPE_REF[ind],ind
        if( new_type ):
            ATYPE_REF.append( ATYPE[i] )
            ATYPE_MASS.append( float(AMASS[i]) )

            TYPE_CNT = TYPE_CNT + 1
            ATYPE_IND.append( TYPE_CNT )
            if(debug): print ' new type found ',i,ELN[i],ATYPE[i],TYPE_CNT,ATYPE_MASS[TYPE_CNT],len(ATYPE_REF)

    debug = 0
    if(debug): sys.exit('lmp_types')

    #
    # Count bonds types
    #
    if( len(BONDS)  > 0):
        TYPE_CNT = 0
        AT_i = ATYPE[ BONDS[0][0] ]
        AT_j = ATYPE[ BONDS[0][1] ]
        BTYPE_IND.append( TYPE_CNT )
        BTYPE_REF.append( [ AT_i ,  AT_j ] )
        if(debug): print ' New bond type ',AT_i ,  AT_j
            
        debug = 0
        for b_indx in range(  1,len(BONDS) ):
            new_type = 1
            a_i = BONDS[b_indx][0] 
            a_j =  BONDS[b_indx][1] 
            AT_i = ATYPE[a_i ]
            AT_j = ATYPE[ a_j ]
            for ind in range(len( BTYPE_REF ) ):
                if ( AT_i == BTYPE_REF[ind][0] and  AT_j == BTYPE_REF[ind][1] ):
                    new_type = 0
                    BTYPE_IND.append( ind )
                    break
                if ( AT_j == BTYPE_REF[ind][0] and  AT_i == BTYPE_REF[ind][1] ):
                    new_type = 0
                    BTYPE_IND.append( ind )
                    break
                
                # If new 
            if( new_type ):
                TYPE_CNT = TYPE_CNT + 1
                BTYPE_IND.append( TYPE_CNT )
                BTYPE_REF.append( [ AT_i ,  AT_j ] )
                if(debug): print ' New bond type ',AT_i ,  AT_j,' type ',TYPE_CNT
                
            if(debug): print b_indx,AT_i ,  AT_j,' type ',BTYPE_IND[b_indx],a_i,a_j

    if(debug): sys.exit('BTYPE_IND')
    #
    # Count angles types
    #
    if( len(ANGLES) > 0):
        TYPE_CNT = 0
        AT_i = ATYPE[ ANGLES[0][0] ]
        AT_j = ATYPE[ ANGLES[0][1] ]
        AT_k = ATYPE[ ANGLES[0][2] ]
        ANGTYPE_IND.append( TYPE_CNT )
        ANGTYPE_REF.append( [ AT_i ,  AT_j , AT_k ] )
        # 
        for ang_indx in range( 1,len(ANGLES) ):
            new_type = 1
            AT_i = ATYPE[ ANGLES[ang_indx][0] ]
            AT_j = ATYPE[ ANGLES[ang_indx][1] ]
            AT_k = ATYPE[ ANGLES[ang_indx][2] ]
            for ind in range( len( ANGTYPE_REF ) ):
                if ( AT_i == ANGTYPE_REF[ind][0] and  AT_j == ANGTYPE_REF[ind][1] and  AT_k == ANGTYPE_REF[ind][2] ):
                    new_type = 0
                    ANGTYPE_IND.append( ind )
                    break
    
                if ( AT_k == ANGTYPE_REF[ind][0] and  AT_j == ANGTYPE_REF[ind][1] and  AT_i == ANGTYPE_REF[ind][2]  ):
                    new_type = 0
                    ANGTYPE_IND.append( ind )
                    break
    
            # If new 
            if( new_type ):
                TYPE_CNT = TYPE_CNT + 1
                ANGTYPE_IND.append( TYPE_CNT )
                ANGTYPE_REF.append( [ AT_i ,  AT_j , AT_k ] )
                a_i = ANGLES[ang_indx][0]
                a_j = ANGLES[ang_indx][1]
                a_k = ANGLES[ang_indx][2]
                if(debug):
                    print  ' new angle type ',AT_i ,  AT_j , AT_k
                    print '         ',a_i ,a_j,a_k#,GTYPE[a_i],GTYPE[a_j],GTYPE[a_k]

    #
    # Count dihedrals types
    #
    debug = 0
    if( len(DIH) > 0 ):
        TYPE_CNT = 0
        AT_i = ATYPE[ DIH[0][0] ]
        AT_j = ATYPE[ DIH[0][1] ]
        AT_k = ATYPE[ DIH[0][2] ]
        AT_l = ATYPE[ DIH[0][3] ]

        DTYPE_IND.append( TYPE_CNT )
        DTYPE_REF.append( [ AT_i ,  AT_j , AT_k , AT_l ] )
        
        if(debug): print ' new dihedral type ',TYPE_CNT,AT_i ,  AT_j , AT_k , AT_l
        #
            
        for dih_indx in range( 1,len(DIH) ):
            new_type = 1
            AT_i = ATYPE[ DIH[dih_indx][0] ]
            AT_j = ATYPE[ DIH[dih_indx][1] ]
            AT_k = ATYPE[ DIH[dih_indx][2] ]
            AT_l = ATYPE[ DIH[dih_indx][3] ]
    
            if(debug):
                print "  Checking ",AT_i ,  AT_j , AT_k , AT_l 
    
            for ind in range( len( DTYPE_REF ) ):
                if ( DTYPE_REF[ind][0] == AT_i   and  AT_j == DTYPE_REF[ind][1] and  AT_k == DTYPE_REF[ind][2] and  AT_l == DTYPE_REF[ind][3]   ):
                    new_type = 0
                    DTYPE_IND.append( ind )
                    break
                if ( DTYPE_REF[ind][0] == AT_l and  AT_k == DTYPE_REF[ind][1] and  AT_j == DTYPE_REF[ind][2] and DTYPE_REF[ind][3] == AT_i   ):
                    new_type = 0
                    DTYPE_IND.append( ind )
                    break
                
            if( new_type ):
                for ind in range( len( DTYPE_REF ) ):
                    
                    if ( DTYPE_REF[ind][0] == 'X' and  AT_j == DTYPE_REF[ind][1] and  AT_k == DTYPE_REF[ind][2] and DTYPE_REF[ind][3] == 'X'   ):
                        new_type = 0
                        DTYPE_IND.append( ind )
                        break
                    if ( DTYPE_REF[ind][0] == 'X' and  AT_j == DTYPE_REF[ind][1] and  AT_k == DTYPE_REF[ind][2] and DTYPE_REF[ind][3] == AT_l   ):
                        new_type = 0
                        DTYPE_IND.append( ind )
                        break
                    if ( DTYPE_REF[ind][0] == AT_i and  AT_j == DTYPE_REF[ind][1] and  AT_k == DTYPE_REF[ind][2] and DTYPE_REF[ind][3] == 'X'   ):
                        new_type = 0
                        DTYPE_IND.append( ind )
                        break
    
                    if ( DTYPE_REF[ind][0] == 'X' and  AT_k == DTYPE_REF[ind][1] and  AT_j == DTYPE_REF[ind][2] and DTYPE_REF[ind][3] == 'X'   ):
                        new_type = 0
                        DTYPE_IND.append( ind )
                        break
                    if ( DTYPE_REF[ind][0] == 'X' and  AT_k == DTYPE_REF[ind][1] and  AT_j == DTYPE_REF[ind][2] and DTYPE_REF[ind][3] == AT_i  ):
                        new_type = 0
                        DTYPE_IND.append( ind )
                        break
                    if ( DTYPE_REF[ind][0] == AT_l and  AT_k == DTYPE_REF[ind][1] and  AT_j == DTYPE_REF[ind][2] and DTYPE_REF[ind][3] == 'X'   ):
                        new_type = 0
                        DTYPE_IND.append( ind )
                        break
                        
            # If new 
            if( new_type ):
                TYPE_CNT += 1
                DTYPE_IND.append( TYPE_CNT )
                DTYPE_REF.append( [ AT_i ,  AT_j , AT_k , AT_l ] )
            
                if( debug ):
                    print ' new dihedral type ',TYPE_CNT,AT_i ,  AT_j , AT_k , AT_l #," = ",DTYPE_REF[ind
                    #print '     '
                
            else:
            #   print '      dihedral type ',ind,AT_i ,  AT_j , AT_k , AT_l ," = ",DTYPE_REF[ind]
                if( debug ):print '  type found ',DTYPE_IND[dih_indx],DIH[dih_indx],ind,DTYPE_REF[ind]
             
             
    if( debug ):
        sys.exit( "  checkin dih types ")
            
    debug = 0
    if( debug ):
        
            
        for dih_indx in range( len(DIH) ):
            new_type = 1
            AT_i = ATYPE[ DIH[dih_indx][0] ]
            AT_j = ATYPE[ DIH[dih_indx][1] ]
            AT_k = ATYPE[ DIH[dih_indx][2] ]
            AT_l = ATYPE[ DIH[dih_indx][3] ]
            type_cnt = DTYPE_IND[dih_indx]
            print '  dihedral type ',type_cnt,AT_i ,  AT_j , AT_k , AT_l,DTYPE_REF[ type_cnt ]
        sys.exit('lmp_types')
        
    debug = 0
    if( debug ):
        for i in range(n_atoms):
            print i,ATYPE[i],ATYPE_IND[i]

        print len( ATYPE_REF), len(ATYPE_MASS)
        for i in range( len(ATYPE_MASS)):
            print i,ATYPE_REF[i],ATYPE_MASS[i]

        print ' angle types ',len(ANGTYPE_REF) 
        for i in range( len(ANGTYPE_REF) ):
            print i,ANGTYPE_REF[i]
            
        print ' dihedra types ',len(DTYPE_REF) 
        for i in range( len(DTYPE_REF) ):
            print i,DTYPE_REF[i]
            
        for i in range( len(DTYPE_IND) ):
            print i,DTYPE_IND[i] #,DIH[i]

        sys.exit('lmp_types')
        

    return ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF 



def read_data(data_file):
    """
    Read data file
    """
    import file_io
    import sys, numpy
    
    debug = 0
    verbose = 1

    if( file_io.file_exists(data_file )):
        if( debug ): print "Reading in ",data_file
    else:
        print " Specified file ",data_file," does not exisit "
        sys.exit("Invalid file ")
        

    F = open(data_file , 'r' )
    lines = F.readlines()
    F.close()
    
    # initialize variables
    LAT_CONST = numpy.zeros( (3,3) )
    
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
            if( col[2]  == "xlo"  and col[2]  == "xhi"  ):
                LAT_CONST[0,0] = float( col[1] ) - float( col[0] ) 
            if( col[2]  == "ylo"  and col[2]  == "yhi"  ):
                LAT_CONST[1,1] = float( col[1] ) - float( col[0] ) 
            if( col[2]  == "zlo"  and col[2]  == "zhi"  ):
                LAT_CONST[2,2] = float( col[1] ) - float( col[0] )

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
    
    RESN = n_atoms*[0]
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
                
            RESN[ind] = int(col[1])
            ATYPE_IND[ind]  = int(col[2]) - 1
            CHARGES[ind]  = float(col[3])
            R[ind] = numpy.array( [ float(col[4]),float(col[5]),float(col[6])])
            ATYPE[ind]  = "??"
        
            if( cnt_Atoms >=  n_atoms ):
                read_Atoms = 0

        if(read_Bonds and len(col) >= 4 ):
            cnt_Bonds += 1
            ind = int( col[0]) - 1

            if( ind > n_bonds ):
                print sys.exit(" Error in data file index of bonds exceeds number of bonds specified with bonds ")
                
            BTYPE_IND[ind] = int(col[1] ) - 1
            BONDS[ind] = [int(col[2]) - 1 , int(col[3]) - 1 ]
            if( cnt_Bonds >=  n_bonds ):
                read_Bonds = 0
                
            
        if(read_Angles and len(col) >= 5 ):
            cnt_Angles += 1
            ind = int( col[0]) - 1
            ANGTYPE_IND[ind] = int(col[1] ) - 1
            ANGLES[ind] = [int(col[2]) - 1, int(col[3]) - 1, int(col[4]) - 1 ]
            if( cnt_Angles >=  n_angles ):
                read_Angles = 0
                
            
        if(read_Dihedrals and len(col) >= 6 ):
            cnt_Dihedrals += 1
            ind = int( col[0]) - 1
            
            DTYPE_IND[ind] = int(col[1] ) - 1
            DIH[ind] = [int(col[2]) - 1,int(col[3]) - 1, int(col[4]) - 1,int(col[5]) - 1]
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
    return (ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
            BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
            ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
            DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,RESN,
            ATYPE_IND,CHARGES,R , ATYPE, BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LAT_CONST)
    #

def print_lmp(data_file, ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
              BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
              ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
              DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
              RESN,ATYPE_IND,CHARGES,R , ATYPE,
              BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LAT_CONST):
    """
    Write data file
    """
    import sys, math
    # 
    #
    # ' print lammps data file '
    #

    # Calculate totals
    n_atoms = len( ATYPE )
    n_bonds = int(len(BONDS))
    n_angles = int(len(ANGLES))
    n_dihedrals = int(len(DIH))
    n_impropers = 0
    n_atypes = int(len( ATYPE_REF)) #+ 1
    n_btypes = int(len( BTYPE_REF)) #+ 1
    n_angtypes = int(len( ANGTYPE_REF)) #+ 1
    n_dtypes = int(len( DTYPE_REF)) #+ 1
    n_imptypes = 1
    # Calculate box size
    bmn_x = LAT_CONST[0][0]/-2.0
    bmx_x = LAT_CONST[0][0]/2.0
    bmn_y = LAT_CONST[1][1]/-2.0
    bmx_y = LAT_CONST[1][1]/2.0
    bmn_z = LAT_CONST[2][2]/-2.0
    bmx_z = LAT_CONST[2][2]/2.0

    F = open( data_file, 'w' )
    F.write('  Lammps data file \n')
    F.write('\n')
    F.write( "%10d  atoms \n" % n_atoms )
    F.write( "%10d  bonds \n" %  n_bonds )
    F.write( "%10d  angles \n" % n_angles )
    F.write( "%10d  dihedrals \n" %  n_dihedrals )
    F.write( "%10d  impropers \n" % n_impropers  )
    F.write('\n')
    F.write( "%10d  atom types \n" % n_atypes  )
    F.write( "%10d  bond types \n" % n_btypes )
    F.write( "%10d  angle types \n" % n_angtypes )
    F.write( "%10d  dihedral types \n" % n_dtypes )
    F.write( "%10d  improper types \n" % n_imptypes )
    F.write('\n')
    F.write( "%16.8f %16.8f   xlo xhi \n" %  (bmn_x , bmx_x) )
    F.write( "%16.8f %16.8f   ylo yhi \n" %  (bmn_y , bmx_y ) )
    F.write( "%16.8f %16.8f   zlo zhi \n" %  (bmn_z , bmx_z) )
    F.write('\n')
    F.write( ' Masses \n')
    F.write('\n')
    for ind in range( len(ATYPE_REF) ):
        AT_i = ATYPE_REF[ind]
        mass = ATYPE_MASS[ind]
        n_mass = ind + 1
        print n_mass , mass, ATYPE_REF[ind]
        F.write( "%10d %16.8f   # %5s \n" % ( n_mass , mass , AT_i ) )
    F.write('\n')
    F.write(' Pair Coeffs \n')
    F.write('\n')
    for i_ind in range( len(ATYPE_REF) ):
        eps_i = ATYPE_EP[i_ind]
        sig_i = ATYPE_SIG[i_ind]
        F.write( "%10d %12.6f %12.6f  \n" % (i_ind + 1,eps_i, sig_i  ) )
    F.write('\n')
    F.write(' Bond Coeffs \n')
    F.write('\n')
    for b_ind in range(  len( BTYPE_REF )  ):
        AT_i =  BTYPE_REF[b_ind][0]
        AT_j =  BTYPE_REF[b_ind][1]
        k = BONDTYPE_K[b_ind]
        r_0 = BONDTYPE_R0[b_ind]
        F.write( "%10d %12.6f %12.6f # %5s %5s  \n" % (b_ind + 1, k , r_0, AT_i, AT_j ) )
    F.write('\n')
    F.write(' Angle Coeffs \n')
    F.write('\n')
    for a_ind in range( len( ANGTYPE_REF ) ):
        AT_i =  ANGTYPE_REF[a_ind][0]
        AT_j =  ANGTYPE_REF[a_ind][1]
        AT_k =  ANGTYPE_REF[a_ind][2]
        r_0 = ANGLETYPE_R0[a_ind]
        k = ANGLETYPE_K[a_ind]
        F.write( "%10d %16.8f %16.8f  # %5s %5s  %5s \n" % (a_ind + 1, k , r_0, AT_i, AT_j, AT_k  ) )
    F.write('\n')
    F.write(' Dihedral Coeffs \n')
    F.write('\n')
    for d_ind in range( len( DTYPE_REF ) ):
        AT_i =  DTYPE_REF[d_ind][0]
        AT_j =  DTYPE_REF[d_ind][1]
        AT_k =  DTYPE_REF[d_ind][2]
        AT_l =  DTYPE_REF[d_ind][3]
        func_type = DIHTYPE_F[d_ind]
        if( func_type == 4 or func_type == 9 ):  # Harmonic
            # phase = DIHTYPE_PHASE[d_ind]
            phase = +1
            k = DIHTYPE_K[d_ind]
            pn = DIHTYPE_PN[d_ind]
            F.write( "%10d %16.8f %5d %10d  # %5s %5s  %5s %5s  \n" % (d_ind + 1, k ,phase , pn , AT_i, AT_j, AT_k ,AT_l  ) )
        elif( func_type == 3 ): # Ryckert-Bellemans
            F1 = DIHTYPE_C[d_ind][0]
            F2 = DIHTYPE_C[d_ind][1]
            F3 = DIHTYPE_C[d_ind][2]
            F4 = DIHTYPE_C[d_ind][3]
            F.write( "%10d  %12.6f  %12.6f  %12.6f  %12.6f # %5s %5s  %5s %5s \n" % (d_ind + 1,F1,F2,F3,F4, AT_i, AT_j, AT_k ,AT_l  ) ) # DIHTYPE_C[d_ind][0:3]  ) )
            
    F.write('\n')
    F.write(' Improper Coeffs \n')
    F.write('\n')
    imp_cnt = 0
    for d_ind in range( len( DTYPE_REF ) ):
        func_type = DIHTYPE_F[d_ind]
        if( func_type == 2 ):
            imp_cnt = imp_cnt + 1
            k = DIHTYPE_K[d_ind]
            ki = DIHTYPE_PHASE[d_ind]
            F.write( "%10d %16.8f %12.6  \n" % (d_ind + 1, k , ki  ) )
            
    if( imp_cnt == 0 ):
        F.write( " 1   0.0 0.0 \n")
    F.write('\n')
    F.write(' Atoms \n')
    F.write('\n')
    TOTAL_CHARGE = 0.0
    for i in range(n_atoms):
        F.write( "%9d %9d %8d %12.8f %12.6f %12.6f %12.6f # %5s \n" % (i+1,RESN[i],ATYPE_IND[i]+1,CHARGES[i],R[i][0],R[i][1],R[i][2] , ATYPE[i] )  )
        TOTAL_CHARGE = TOTAL_CHARGE + float( CHARGES[i] )
        
    print ' Total charge = ',TOTAL_CHARGE
    F.write('\n')
    F.write(' Bonds \n')
    F.write('\n')
    for i in range( len(BONDS) ):
        #
        a_i = BONDS[i][0] 
        a_j = BONDS[i][1]
        #
        AT_i =  ATYPE[a_i]
        AT_j =  ATYPE[a_j]
        #
        b_ind = BTYPE_IND[i]  + 1
        
        F.write(  '%9d %8d %9d %9d # %5s %5s \n' % (i+1,b_ind,a_i+1,a_j+1 , AT_i, AT_j ) )
    F.write('\n')
    F.write(' Angles \n')
    F.write('\n')
    for i in range( len(ANGLES) ):
        a_i = ANGLES[i][0] + 1
        a_j = ANGLES[i][1] + 1
        a_k = ANGLES[i][2] + 1
        a_ind = ANGTYPE_IND[i]+ 1 
        F.write(  '%9d %8d %9d %9d %9d \n' % (i+1,a_ind,a_i,a_j,a_k) )
    F.write(  '\n' )
    F.write(' Dihedrals \n')
    F.write('\n')
    for i in range( len(DIH) ):
        a_i = DIH[i][0] + 1
        a_j = DIH[i][1] + 1
        a_k = DIH[i][2] + 1
        a_l = DIH[i][3] + 1
        d_ind = DTYPE_IND[i]+ 1 
        if( DIHTYPE_F[d_ind - 1] != 2 ):
            F.write(  '%9d %8d %9d %9d %9d  %9d \n' % (i+1,d_ind,a_i,a_j,a_k,a_l) )
    F.write( '\n' )
    F.write(' Impropers \n')
    for i in range( len(DIH) ):
        a_i = DIH[i][0] + 1
        a_j = DIH[i][1] + 1
        a_k = DIH[i][2] + 1
        a_l = DIH[i][3] + 1
        d_ind = DTYPE_IND[i]+ 1 
        if( DIHTYPE_F[d_ind - 1] == 2 ):
            F.write(  '%9d %8d %9d %9d %9d  %9d \n' % (i+1,d_ind,a_i,a_j,a_k,a_l) )
            
    F.close()
    

def print_rest(rest_file,data_file,cent_angle,dih_indx,options):
    """
    Write input file for relaxed torsional angle calculation 
    """
    if( len(dih_indx) > 3 ):
        dihindx_line = str( dih_indx[0] ) + " " + str( dih_indx[1] ) + " " + str( dih_indx[2] ) + " " + str( dih_indx[3] )+ " " 
    else:
        print " less than 3 atoms specified in dihedral "
        sys.exit("print_rest ")
        
    # LAMMPS input: minimization w/ soft potential  
    input_lines = ' units 		real ' 
    input_lines = input_lines + '\n'+ ' boundary 	p p p	'
    input_lines = input_lines + '\n'+ ' newton          on  # use off for np>4 '
    input_lines = input_lines + '\n'+ ' atom_style	full	'
    input_lines = input_lines + '\n'+ ' bond_style      harmonic'
    input_lines = input_lines + '\n'+ ' angle_style     harmonic'
    input_lines = input_lines + '\n'+ ' #dihedral_style  harmonic # hybrid multi/harmonic'
    input_lines = input_lines + '\n'+ ' dihedral_style opls'
    input_lines = input_lines + '\n'+ ' improper_style  harmonic	'
    if( options.zero_charges ):
        input_lines = input_lines + '\n'+ ' pair_style 	lj/cut  20.0'
    else:
        input_lines = input_lines + '\n'+ ' pair_style 	lj/cut/coul/long 20.0'
        input_lines = input_lines + '\n'+ ' kspace_style    pppm    1.0e-5'
    input_lines = input_lines + '\n'+ ' pair_modify     mix geometric '
    input_lines = input_lines + '\n'+ ' # special_bonds lj/coul 0.0 0.0 0.5 #angle yes dihedral yes    '
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' neighbor        2.0     bin'
    input_lines = input_lines + '\n'+ ' neigh_modify    delay   5'
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' read_data ' + data_file
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' fix TFIX all langevin 0.0 0.0 100 24601'
    
    
    input_lines = input_lines + '\n'+ ' fix REST all restrain dihedral ' + dihindx_line + ' 2000.0 5000.0  '+str(float(cent_angle))
    input_lines = input_lines + '\n'+ ' fix_modify REST energy no'
    input_lines = input_lines + '\n'+ ' run 10000'
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' dump 1 all xyz 1 fit.xyz '
    input_lines = input_lines + '\n'+ ' compute 1 all pe'
    input_lines = input_lines + '\n'+ ' thermo_style custom step temp pe etotal press vol'
    input_lines = input_lines + '\n'+ ' write_data  out_'  + data_file 
    input_lines = input_lines + '\n'+ ' # sanity check for convergence'
    input_lines = input_lines + '\n'+ ' minimize 1e-6 1e-9 1 100000'
    
    
    f = open(rest_file,'w')
    f.write(input_lines)
    f.close()

    return

def print_rest2(rest_file,data_file,DIH_CONST_ANGLE,DIH_CONST_ATOMS,options):
    """
    Write input file for relaxed torsional angle calculation  version 2
    """
    
    # LAMMPS input: minimization w/ soft potential  
    input_lines = ' units 		real ' 
    input_lines = input_lines + '\n'+ ' boundary 	p p p	'
    input_lines = input_lines + '\n'+ ' newton          on  # use off for np>4 '
    input_lines = input_lines + '\n'+ ' atom_style	full	'
    input_lines = input_lines + '\n'+ ' bond_style      harmonic'
    input_lines = input_lines + '\n'+ ' angle_style     harmonic'
    input_lines = input_lines + '\n'+ ' #dihedral_style  harmonic # hybrid multi/harmonic'
    input_lines = input_lines + '\n'+ ' dihedral_style opls'
    input_lines = input_lines + '\n'+ ' improper_style  harmonic	'
    if( options.zero_charges ):
        input_lines = input_lines + '\n'+ ' pair_style 	lj/cut  20.0'
    else:
        input_lines = input_lines + '\n'+ ' pair_style 	lj/cut/coul/long 20.0'
        input_lines = input_lines + '\n'+ ' kspace_style    pppm    1.0e-5'
    input_lines = input_lines + '\n'+ ' pair_modify     mix geometric '
    input_lines = input_lines + '\n'+ ' # special_bonds lj/coul 0.0 0.0 0.5 #angle yes dihedral yes    '
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' neighbor        2.0     bin'
    input_lines = input_lines + '\n'+ ' neigh_modify    delay   5'
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' read_data ' + data_file
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' fix TFIX all langevin 0.0 0.0 100 24601'
    
    for d_indx in range(len(DIH_CONST_ANGLE)):
        dihindx_line = ""
        dih_atoms =    DIH_CONST_ATOMS[d_indx]
        for a_indx in dih_atoms:
            dihindx_line += " "+str(a_indx)
        dih_angle = DIH_CONST_ANGLE[d_indx]
    
        fix_id = "REST"+str(d_indx)
        input_lines = input_lines + '\n'+ ' fix '+fix_id+' all restrain dihedral ' + dihindx_line + ' 2000.0 5000.0  '+str(float(dih_angle))
        input_lines = input_lines + '\n'+ ' fix_modify '+fix_id+' energy no'
        
        
    input_lines = input_lines + '\n'+ ' run 10000'
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' dump 1 all xyz 1 fit.xyz '
    input_lines = input_lines + '\n'+ ' compute 1 all pe'
    input_lines = input_lines + '\n'+ ' thermo_style custom step temp pe etotal press vol'
    input_lines = input_lines + '\n'+ ' write_data  out_'  + data_file 
    input_lines = input_lines + '\n'+ ' # sanity check for convergence'
    input_lines = input_lines + '\n'+ ' minimize 1e-6 1e-9 1 100000'
    
    
    f = open(rest_file,'w')
    f.write(input_lines)
    f.close()

    return 

def print_sp(rest_file,data_file,options):
    """
    Write input file for single point torsional angle calculation
    """
        
    # LAMMPS input: minimization w/ soft potential  
    input_lines = ' units 		real ' 
    input_lines = input_lines + '\n'+ ' boundary 	p p p	'
    input_lines = input_lines + '\n'+ ' newton          on  # use off for np>4 '
    input_lines = input_lines + '\n'+ ' atom_style	full	'
    input_lines = input_lines + '\n'+ ' bond_style      harmonic'
    input_lines = input_lines + '\n'+ ' angle_style     harmonic'
    input_lines = input_lines + '\n'+ ' #dihedral_style  harmonic # hybrid multi/harmonic'
    input_lines = input_lines + '\n'+ ' dihedral_style opls'
    input_lines = input_lines + '\n'+ ' improper_style  harmonic	'
    
    if( options.zero_charges ):
        input_lines = input_lines + '\n'+ ' pair_style 	lj/cut  20.0'
    else:
        input_lines = input_lines + '\n'+ ' pair_style 	lj/cut/coul/long 20.0'
        input_lines = input_lines + '\n'+ ' kspace_style    pppm    1.0e-5'
        
    input_lines = input_lines + '\n'+ ' pair_modify     mix geometric '
    input_lines = input_lines + '\n'+ ' # special_bonds lj/coul 0.0 0.0 0.5 #angle yes dihedral yes    '
    input_lines = input_lines + '\n'+ ' # kspace_style    pppm    1.0e-5'
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' neighbor        2.0     bin'
    input_lines = input_lines + '\n'+ ' neigh_modify    delay   5'
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' read_data ' + data_file
    input_lines = input_lines + '\n'+ ' '
    input_lines = input_lines + '\n'+ ' fix TFIX all langevin 0.0 0.0 100 24601'
    input_lines = input_lines + '\n'+ ' dump 1 all xyz 1 sp.xyz '    
    input_lines = input_lines + '\n'+ ' thermo_style custom step temp pe etotal press vol'
    input_lines = input_lines + '\n'+ ' run 0'

    f = open(rest_file,'w')
    f.write(input_lines)
    f.close()

    return 

def run_lmp(lammps_dir,lammps_src,lmp_in):
    """
    Run lammps
    """
    
    import sys, os
    from string import replace
    
    run_id = replace(lmp_in,'.in','')

    run_line = lammps_dir + lammps_src + " < " + lmp_in + " > lmp.out "
    print run_line 
    os.system(run_line)
    mv_line = "mv log.lammps " + run_id+".log"
    os.system(mv_line)
    
    return


def get_pe(log_file):
    """
    Get potential energy from log file 
    """
    import sys, os

    KCALtoEV = 0.04336411

    debug = 0

    pe = []
    
    # Open log file and get energy
    print " reading in ",log_file
    f = open(log_file,'r')
    Lines = f.readlines()
    f.close()
    
    
    read_line = 0
    pot_en = "na"
    for line in Lines :
        col = line.split()
        
        if( len(col) > 4 ):
            if( col[0] == "Loop" and col[1] == "time"  ):
                read_line = 0
                
        if( read_line ):
            print " readin in energy ",col[2] 
            pot_en = float( col[2] )*KCALtoEV
            pe.append( pot_en )
            
        if( len(col) > 4 ):
            if( col[2] == "PotEng" ):
                read_line = 1
                
                
    
    if(debug):
        print " potential energies ",pot_en
        for indx in range( len(pe)):
            print indx , pe[indx]
            
        sys.exit("get_pe")

    return pot_en   # in eV 


def last_xmol(xmol_file,options):
    """
    Get file geometry from xmol file 
    """
    debug = 0
    
    F = open(xmol_file,'r')
    Lines = F.readlines()
    F.close()
    
    col = Lines[0].split()
    n_atoms = int( col[0] ) 

    R = []
    line_cnt = 0 
    for line in Lines :
        col = line.split()
        line_cnt = line_cnt + 1
        if( line_cnt > (n_atoms + 2 ) ):
            line_cnt = 1
            R = []
            if( debug ): print " restarting ", line
        if( line_cnt > 2 ):
            R.append( [ float( col[1]), float( col[2]), float( col[3]) ])
            if( debug ): print float( col[1]), float( col[2]), float( col[3]) 
        
    if( debug ):
        for indx in range( len(R)): 
            print indx,R[indx]
        sys.exit("last_xmol")
        
    return R

    
def get_dihangle(R,angle_indx,options):
    """
    Calcuate dihedral angle of 4 atoms
        <-phi->
     I         L
      \      / 
        J - K    
    """
    
    import numpy , prop
    
    debug = 0
    if( debug ):
        for indx in range( len(R)): 
            print indx,R[indx]
            
    i = angle_indx[0]
    a = numpy.array( [ R[i][0],R[i][1],R[i][2] ] ) 
    j = angle_indx[1]
    b= numpy.array( [ R[j][0],R[j][1],R[j][2] ] ) 
    k = angle_indx[2]
    c= numpy.array(  [ R[k][0],R[k][1],R[k][2] ] ) 
    l = angle_indx[3]
    d= numpy.array(  [ R[l][0],R[l][1],R[l][2] ] ) 
    
    
    angle = prop.getDihedral(a,b,c,d)
    
    return angle

    
