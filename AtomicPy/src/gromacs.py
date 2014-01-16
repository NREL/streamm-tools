#! /usr/bin/env python
# subroutines for gromacs file manipulation

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# Email travis.kemper@nrel.gov
# Version 3.00 

kJtoeV = 0.01036427
Kjtokcal = 0.239006 



def read_top(options,top_infile):
    # Read in topology file 
    
    # General python modules 
    import sys, numpy
    
    # atomicpy functions
    import file_io

    debug = 0

    if( file_io.file_exists(top_infile )):
        if( debug  ): print "      Reading in ",top_infile
    else:
        print " Specified .gro file ",top_infile," does not exisit "
        sys.exit("Invalid file ")

    F = open(top_infile , 'r' )
    top_lines = F.readlines()
    F.close()
    
    #
    # Compile included files into single array
    #
    Lines = []
    for line in top_lines:
        col = line.split()
        if ( len(col) > 1 ):
            if ( col[0][:1] != ';' ):
                if( col[0] == "#include" ):
                    include_file = col[1].strip("\"")

                    F_inc = open(include_file , 'r' )
                    include_lines = F_inc.readlines()
                    F_inc.close()
                    
                    if( options.verbose ):
                        print "      Including ",include_file
                        
                    for l_inc in include_lines:
                        c_inc = l_inc.split()
                        if ( len(c_inc) > 1 ):
                            if ( c_inc[0][:1] != ';' ):
                                Lines.append(l_inc)
                                
                                
                else:
                    Lines.append(line)

                    
    #
    # Get moleculetypes and count 
    #
    read_MOLECULETYPE = 0
    read_MOLECULES = 0 

    MOLECULETYPE = []
    MOLECULECNT = []
    MOLECULECNT_ID = []
    MOL_CNT = 0 
    
    for line in Lines:
        col = line.split()
        if ( len(col) > 1 ):
            # Read in molecule types
            if read_MOLECULETYPE :
                if ( len(col) >= 2 ):
                    if ( col[0][:1] != ';' and  col[0][:1] != '[' ):
                        MOLECULETYPE.append( col[0])
                if ( col[0][:1] == '[' ): read_MOLECULETYPE = 0
                        
            # Read in how many times those molecule types are repleated 
            if read_MOLECULES :
                if ( len(col) >= 2 ):
                    if ( col[0][:1] != ';' and  col[0][:1] != '[' ):
                        MOLECULECNT_ID.append( col[0])
                        mol_n = int(col[1])
                        MOLECULECNT.append( mol_n )
                        MOL_CNT += mol_n
                if ( col[0][:1] == '[' ): read_MOLECULES = 0
                        
            if ( col[0][:1] == '[' and col[1] == 'moleculetype' ):
                read_MOLECULETYPE = 1
            if ( col[0][:1] == '[' and col[1] == 'molecules' ):
                read_MOLECULES = 1

    
    
    if( options.verbose ):
        print "      Molecule types found ",len(MOLECULETYPE) 
        for mol_indx in range( len(MOLECULETYPE) ):
            print "         ",mol_indx+1,MOLECULETYPE[mol_indx]
        
        print "      Total molecule types in system ",MOL_CNT
        for mol_indx in range( len(MOLECULECNT_ID) ):
            print "         ",mol_indx+1,MOLECULECNT_ID[mol_indx],MOLECULECNT[mol_indx]
    
 
    
    # Declare molecular list for reading in molecule 
    MOL_ATOMS_INDEX =  [] #numpy.empty( MOL_CNT+1, dtype=int )
    MOL_BONDS_INDEX =  [] #numpy.empty( MOL_CNT+1, dtype=int )
    MOL_ANGLES_INDEX =  [] #numpy.empty( MOL_CNT+1, dtype=int )
    MOL_DIH_INDEX =  [] #numpy.empty( MOL_CNT+1, dtype=int )

    # Atomic lists form molecule i 
    ATYPE_l = []
    RESN_l = []
    RESID_l = []
    GTYPE_l = []
    CHARN_l = []    
    CHARGES_l = []
    AMASS_l = []

    BONDS_l = []
    ANGLES_l = []
    DIH_l = []

    
    # atoms 
    read_ATOMS = 0
    ATOMS_CNT = -1
    # bonds
    read_BONDS = 0
    BONDS_CNT = -1
    # angles
    read_ANGLES = 0
    ANGLES_CNT = -1
    # dihedrals
    read_DIH = 0
    DIH_CNT = -1
    
    # Loop over lines of top file 
    for line in Lines:
        col = line.split()
        if ( len(col) > 2 and col[0][:1] != ';' ):
            # Read in atoms section to local arrays 
            if read_ATOMS :
                if ( col[0][:1] == '[' ):
                    read_ATOMS = 0
                elif ( len(col) >= 7 ):
                    ATYPE_l.append( col[1] )
                    RESN_l.append(  int(col[2]) )
                    RESID_l.append( col[3] )
                    GTYPE_l.append( col[4] )
                    CHARN_l.append(  int(col[5]) )
                    CHARGES_l.append(  float(col[6]) )
                    AMASS_l.append(  float(col[7]) )
                    ATOMS_CNT = ATOMS_CNT + 1 
                        

            # Read in bonds section to local arrays 
            if read_BONDS :
                if ( col[0][:1] == '[' ):
                    read_BONDS = 0
                elif ( len(col) >= 2 ):
                    BONDS_l.append( [ int(col[0])-1,int(col[1])-1 ] )
                    BONDS_CNT = BONDS_CNT + 1 
                    

            # Read in angles section to local arrays 
            if read_ANGLES :
                if ( col[0][:1] == '[' ):
                    read_ANGLES = 0
                elif ( len(col) >= 2 ):
                    ANGLES_l.append( [int(col[0])-1,int(col[1])-1,int(col[2])-1 ] )
                    ANGLES_CNT = ANGLES_CNT + 1 

            # Read in dihedrals section to local arrays 
            if read_DIH :
                if ( col[0][:1] == '[' ):
                    read_DIH = 0
                elif ( len(col) >= 2 ):
                    DIH_l.append( [int(col[0])-1,int(col[1])-1,int(col[2])-1,int(col[3])-1 ] )
                    DIH_CNT = DIH_CNT + 1
            #
            # switch on atoms read in and tack local arrays by molecule number 
            #
            if ( col[0][:1] == '[' and col[1].strip() == 'moleculetype' ):
                MOL_ATOMS_INDEX.append( ATOMS_CNT + 1 )
                MOL_BONDS_INDEX.append( BONDS_CNT + 1 )
                MOL_ANGLES_INDEX.append( ANGLES_CNT + 1 )
                MOL_DIH_INDEX.append( DIH_CNT + 1 )                
                
            if ( col[0][:1] == '[' and col[1].strip() == 'atoms' ):
                read_ATOMS = 1
            # switch on bonds read in and tack local arrays by molecule number 
            if ( col[0][:1] == '[' and col[1].strip() == 'bonds' ):
                read_BONDS = 1
            # switch on angles read in and tack local arrays by molecule number 
            if ( col[0][:1] == '[' and col[1].strip() == 'angles' ):
                read_ANGLES = 1
            # switch on dihedrals read in and tack local arrays by molecule number 
            if ( col[0][:1] == '[' and col[1].strip() == 'dihedrals' ):
                read_DIH = 1

    # index to account for last 
    MOL_ATOMS_INDEX.append( ATOMS_CNT + 1 )
    MOL_BONDS_INDEX.append( BONDS_CNT + 1 )
    MOL_ANGLES_INDEX.append( ANGLES_CNT + 1 )
    MOL_DIH_INDEX.append( DIH_CNT + 1 )
    
    
    if( options.verbose ):
        print "      Total molecules ",MOL_CNT
        for mol_indx in range( len(MOLECULECNT)):
            print "      Molecule id ",MOLECULECNT_ID[mol_indx]
            if( MOL_ATOMS_INDEX[mol_indx] >= 0 ):
                print "         atoms ",MOL_ATOMS_INDEX[mol_indx],MOL_ATOMS_INDEX[mol_indx+1] - 1
            if( MOL_BONDS_INDEX[mol_indx] >= 0 ):
                print "         bonds ",MOL_BONDS_INDEX[mol_indx],MOL_BONDS_INDEX[mol_indx+1] - 1
            if( MOL_ANGLES_INDEX[mol_indx] >= 0 ):
                print "         angles ",MOL_ANGLES_INDEX[mol_indx],MOL_ANGLES_INDEX[mol_indx+1] - 1
            if( MOL_DIH_INDEX[mol_indx] >= 0 ):
                print "         dihedrals ",MOL_DIH_INDEX[mol_indx],MOL_DIH_INDEX[mol_indx+1] - 1
            if( MOLECULECNT[mol_indx] > 0 ):
                print "         will be repeated ",MOLECULECNT[mol_indx]," times "

    # Declare lists to be global 
    ATYPE = []
    RESN = []
    RESID = []
    GTYPE = []
    CHARN = []    
    CHARGES = []
    AMASS = []
    MOLNUMB = []
    
    BONDS = []
    ANGLES = []
    DIH = []
    
    MOLLIST = []
    MOLPNT = []
    
    
    A_CNT = -1
    MOL_CNT = -1 
    
    
    debug = 0
    
    # Loop over molecules
    for repeat_indx in range( len(MOLECULECNT) ):
        mol_id = MOLECULECNT_ID[repeat_indx]
        # Find molecule type
        id_n = -1
        if( debug): print " searching molecule types ",len(MOLECULETYPE)
        for id_indx in range( len(MOLECULETYPE) ):
            id_type = MOLECULETYPE[id_indx]
            if( id_type == mol_id ):
                id_n = id_indx
        if( id_n < 0 ):
            sys.exit(" error in .top read in")
            
        # Repeat specified times
        if( debug ): print " adding molecule ",id_n,MOLECULECNT[repeat_indx]
        
        for mol_repeat in range( MOLECULECNT[repeat_indx]):
            MOL_CNT += 1
            MOLPNT.append( A_CNT+1 )
            
            # Repeat atoms
            N_o = MOL_ATOMS_INDEX[id_n] #- 1
            N_f = MOL_ATOMS_INDEX[id_n+1] - 1
            if( debug ): print " adding atoms",N_o,N_f
            
            # Record number of atoms to shift 
            NA_o = len( ATYPE)
            REF_N = []  # Reference to from initial atom # to current in added molecule 

            for atom_indx in range(N_o,N_f+1):
                # Add atoms to global list
                A_CNT += 1
                ATYPE.append( ATYPE_l[atom_indx] )
                RESN.append(  RESN_l[atom_indx] )
                RESID.append( RESID_l[atom_indx])
                GTYPE.append( GTYPE_l[atom_indx] )
                CHARN.append(  CHARN_l[atom_indx] )
                CHARGES.append(  CHARGES_l[atom_indx] )
                AMASS.append(  AMASS_l[atom_indx] )
                # Reference atom # in global list
                REF_N.append( A_CNT )
                # Molecule information 
                MOLLIST.append( A_CNT )
                MOLNUMB.append( MOL_CNT)
                if( debug ): print " ADDING ",A_CNT,GTYPE_l[atom_indx]
                
            # Repeat bonds
            N_o = MOL_BONDS_INDEX[id_n] #- 1
            N_f = MOL_BONDS_INDEX[id_n+1] - 1
            if( debug ): print " adding bonds ",N_o,N_f
            if( N_f > N_o ):
                for bond_indx in range(N_o,N_f+1):
                    i_o = BONDS_l[bond_indx][0]
                    j_o = BONDS_l[bond_indx][1]
                    i_add = REF_N[i_o]
                    j_add = REF_N[j_o]
                    BONDS.append( [i_add ,j_add] )
                    #if( debug ): print i_o,j_o , " -> ",i_add,j_add
            # Repeat angles
            N_o = MOL_ANGLES_INDEX[id_n] #- 1
            N_f = MOL_ANGLES_INDEX[id_n+1] - 1
            if( debug ): print " adding angles ",N_o,N_f
            if( N_f > N_o ):
                for angle_indx in range(N_o,N_f+1):
                    k_o = ANGLES_l[angle_indx][0]
                    i_o = ANGLES_l[angle_indx][1]
                    j_o = ANGLES_l[angle_indx][2]
                    k_add = REF_N[k_o]
                    i_add = REF_N[i_o]
                    j_add = REF_N[j_o]
                    ANGLES.append( [k_add,i_add,j_add ])
            # Repeat dihedral
            N_o = MOL_DIH_INDEX[id_n] #- 1
            N_f = MOL_DIH_INDEX[id_n+1] - 1
            if( debug ): print " adding dihedral ",N_o,N_f
            if( N_f > N_o ):
                for dih_indx in range(N_o,N_f+1):
                    k_o = DIH_l[dih_indx][0]
                    i_o = DIH_l[dih_indx][1]
                    j_o = DIH_l[dih_indx][2]
                    l_o = DIH_l[dih_indx][3]
                    k_add = REF_N[k_o]
                    i_add = REF_N[i_o]
                    j_add = REF_N[j_o]
                    l_add = REF_N[l_o]
                    DIH.append([k_add,i_add,j_add,l_add])
                
    if(debug): print " Total atoms ",A_CNT
    MOLPNT.append( A_CNT + 1)
        
            
    return ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST


def read_gro(options,in_gro):
    # Read gromacs structure file
    
    import sys, numpy 

    # atomicpy functions
    import file_io

    debug = 0

    if( file_io.file_exists(in_gro )):
        if( debug ): print "Reading in ",in_gro
    else:
        print " Specified .gro file ",in_gro," does not exisit "
        sys.exit("Invalid file ")
    
    # Declare coordinate array R
    R = []        # list of numpy arrays for each atomic vector postion
    VEL = []      # list of numpy arrays for each atomic vector velocity 
    GTYPE = []    # list of gromacs types 
    LV = []       # lattice vectors 

    F = open(in_gro , 'r' )
    Lines = F.readlines()
    F.close()
    
    # Read in .gro file
    #
    line_cnt = 0
    for line in Lines :
        line_cnt = line_cnt + 1
        if( line_cnt > 2 and len(line) >= 44): # skip header
            g = line[10:15]
            GTYPE.append( g )
            x = line[20:28]
            y = line[28:36]
            z = line[36:44]
            # Convert from nm to angstroms 
            r_i =  numpy.array(  [float(x)*10,float(y)*10,float(z)*10] )
            R.append( r_i )
    
    #        
    # Get lattice vector from last line
    #
    
    # Initialize as zero
    LV = numpy.zeros( (3,3) )
    #for d_i in range(3):
    #    LV.append( [0.0,0.0,0.0])
    
    col = line.split()
    n_vec = int( len(col))
    if( n_vec == 3 ):
        LV[0][0] = float( col[0] )*10
        LV[1][1] = float( col[1] )*10
        LV[2][2] = float( col[2] )*10
    if( n_vec == 9 ):
        LV[0][0] = float( col[0] )*10
        LV[1][1] = float( col[1] )*10
        LV[2][2] = float( col[2] )*10
        LV[0][1] = float( col[3] )*10
        LV[0][2] = float( col[4] )*10
        LV[1][0] = float( col[5] )*10
        LV[1][2] = float( col[6] )*10
        LV[2][0] = float( col[7] )*10
        LV[2][1] = float( col[8] )*10
        
    if( options.verbose ):
        print "      Box size ",LV[0,0],LV[1,1],LV[2,2]," angstorms "

    return ( GTYPE, R, VEL, LV )


def read_itp(ff_file):
    import sys

    FF_ATOMTYPES = []
    FF_BONDTYPES = []
    FF_ANGLETYPES = []
    FF_DIHTYPES = []

    F = open(ff_file , 'r' )
    Lines = F.readlines()
    F.close()
    #
    # Get atom types
    # 
    read_line = 0
    ID='atomtypes'
    for line in Lines:
        col = line.split()
        if ( len(col) > 2 ):
            if ( col[0][:1] == '[' and col[1] == 'atomtypes'):
                read_line = 1
            if ( col[0][:1] == '[' and col[1] != 'atomtypes'):
                 read_line = 0
            if read_line :
                if ( len(col) >= 7 ):
                    if ( col[0][:1] != ';' ):
                        FF_ATOMTYPES.append( col )
                        
    #
    # Get bond types
    # 
    read_line = 0
    ID='bondtypes'
    for line in Lines:
        col = line.split()
        if ( len(col) > 2 ):
            if ( col[0][:1] == '[' and col[1] == ID ):
                read_line = 1
            if ( col[0][:1] == '[' and col[1] != ID ):
                 read_line = 0
            if read_line :
                if ( len(col) >= 5 ):
                    if ( col[0][:1] != ';' ):
                        FF_BONDTYPES.append(col)

    #
    # Get angle types
    # 
    read_line = 0
    ID='angletypes'
    for line in Lines:
        col = line.split()
        if ( len(col) > 2 ):
            if ( col[0][:1] == '[' and col[1] == ID ):
                read_line = 1
            if ( col[0][:1] == '[' and col[1] != ID ):
                 read_line = 0
            if read_line :
                if ( len(col) >= 5 ):
                    if ( col[0][:1] != ';' ):
                        FF_ANGLETYPES.append( col )


    #
    # Get dihedral types
    # 
    read_line = 0
    ID='dihedraltypes'
    for line in Lines:
        col = line.split()
        if ( len(col) > 2 ):
            if ( col[0][:1] == '[' and col[1] == ID ):
                read_line = 1
            if ( col[0][:1] == '[' and col[1] != ID ):
                read_line = 0
                
            if read_line :
                if ( len(col) > 6 ):
                    if ( col[0][:1] != ';' ):
                        FF_DIHTYPES.append( col )


    return ( FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES )


def print_gro(gro_file,GTYPE,RESID,RESN,R,LV):
    import sys 
    #
    # printing gro 
    #
    F = open( gro_file, 'w' )
    NA = len(GTYPE)
    F.write( ' log to top ')
    F.write( "\n %-2i "  % int(NA) )
    atom_indx = 0 
    for i in range(NA):
        atom_indx += 1
        if( atom_indx > 10000): atom_indx = 1
        F.write( "\n%5d%-5s%5s%5d%8.3f%8.3f%8.3f"  % ( RESN[i],RESID[i],GTYPE[i],atom_indx,R[i][0]/10.0,R[i][1]/10.0,R[i][2]/10.0 ))
    F.write( "\n %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n" % (LV[0][0]/10.0,LV[1][1]/10.0,LV[2][2]/10.0,LV[0][1]/10.0,LV[0][2]/10.0,LV[1][0]/10.0,LV[1][2]/10.0,LV[2][0]/10.0,LV[2][1]/10.0))
    F.close()


def print_top( top_file,ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
               ,DIH_CONST,DIH_CONST_ANGLE
               ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
               ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LAT_CONST):
    import sys
    
    NA = len( ELN )
    debug = 0 ;
    print
    print ' print top file '
    print
    
    spc = str(' ')
    
    F = open( top_file, 'w' )
    #
    # write atoms 
    #
    F.write('; top file for gromacs \n ')
    F.write(' [  defaults ] \n' )
    F.write(' ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ  \n')
    F.write(' 1               2               yes             0.5     0.8333 \n')
    F.write('\n')
    F.write('; force field file \n')
    F.write('#include "ff-new.itp"  \n')
    F.write('\n')
    F.write( '[ moleculetype ] \n')
    F.write( ' MOL           3 \n')
    F.write('\n')
    F.write( '[ atoms ] \n' )
    TOTAL_CHARGE = 0.0
    for i in range(len(ASYMB)):
        # print i,ATYPE[i] ,RESN[i] ,RESID[i] ,GTYPE[i],CHARN[i],CHARGES[i] ,AMASS[i]
        F.write( "%5d %5s %5d %10s %5s %5d %16.12f %12.6f  \n" % (i+1,ATYPE[i],RESN[i],RESID[i],GTYPE[i],CHARN[i],CHARGES[i],AMASS[i]) )
        TOTAL_CHARGE = TOTAL_CHARGE + float( CHARGES[i] )

    print ' Total charge = ',TOTAL_CHARGE
    F.write( '' + '\n' )
    #
    # write bonds
    #
    F.write( ' [ bonds ] \n' )
    print len(BONDS)-1 , ' number of bonds'
    for i in range( len(BONDS) ):
        a_i = BONDS[i][0] + 1
        a_j = BONDS[i][1] + 1
        b_ind = BTYPE_IND[i] 
        b_type = BONDTYPE_F[b_ind]
        F.write(  '%10d %10d %5d \n' % (a_i,a_j,b_type) )
    F.write(  '\n' )
    #
    # write angles
    #
    F.write( ' [ angles ] \n' )
    print len(ANGLES)-1 , ' number of angles'
    for i in range( len(ANGLES) ):
        a_i = ANGLES[i][0] + 1
        a_j = ANGLES[i][1] + 1
        a_k = ANGLES[i][2] + 1
        ind = ANGTYPE_IND[i] 
        func_type = ANGLETYPE_F[ind]
        F.write(  '%10d %10d %10d %5d \n' % (a_i,a_j,a_k,func_type) )
    F.write(  '\n' )
    #
    # write dihedrals 
    #
    debug = 0
    line = str( '[ dihedrals ]' )
    F.write( line + '\n' )        
    print len(DIH)-1 , ' number of dihedrals '
    for i in range( len(DIH) ):
        a_i = DIH[i][0] + 1
        a_j = DIH[i][1] + 1
        a_k = DIH[i][2] + 1
        a_l = DIH[i][3] + 1
        ind = DTYPE_IND[i] 
        if(debug): print ' dihedral index ',ind, a_i,a_j,a_k,a_l,ATYPE[a_i-1],ATYPE[a_j-1],ATYPE[a_k-1],ATYPE[a_l-1]
        func_type = DIHTYPE_F[ind]
        F.write(  '%10d %10d %10d %10d %5d \n' % (a_i,a_j,a_k,a_l,func_type) )
    if(debug): sys.exit(" print to pdih ")
    
    #
    if( len(DIH_CONST) > 0 ):
        F.write( '\n' )
        F.write( '[ dihedral_restraints ] \n' )
        F.write( '; ai   aj    ak    al  type  label  phi  dphi  kfac  power \n' )
        const_type = 1
        const_phi = DIH_CONST_ANGLE
        const_dphi = 0
        const_kfac = 20000.0
        #for indx_const in range( len( DIH_CONST) ):
        #    print DIH_CONST[indx_const]
        #    print DIH_CONST[indx_const][0],DIH_CONST[indx_const][1],DIH_CONST[indx_const][2],DIH_CONST[indx_const][3]
        #    print const_type,const_phi,const_dphi,const_kfac
        #    F.write( " %8d %8d %8d %8d %8d %f12.6 %8d %16.6 " % (DIH_CONST[indx_const][0],DIH_CONST[indx_const][1],DIH_CONST[indx_const][2],DIH_CONST[indx_const][3],const_type,const_phi,const_dphi,const_kfac))
        #F.write( " %8d %8d %8d %8d %8d %12.6f %8d %16.6f \n" % (DIH_CONST[0],DIH_CONST[1],DIH_CONST[2],DIH_CONST[3],const_type,const_phi,const_dphi,const_kfac))
        F.write( " %8d %8d %8d %8d  %8d 1 %12.6f %8d %16.6f  1 \n" % (DIH_CONST[0],DIH_CONST[1],DIH_CONST[2],DIH_CONST[3],const_type,const_phi,const_dphi,const_kfac))
    #
    # write improper dihedrals 
    #
    debug = 0
    for i in range( len(IMPS) ):
        a_i = IMPS[i][0] + 1
        a_j = IMPS[i][1] + 1
        a_k = IMPS[i][2] + 1
        a_l = IMPS[i][3] + 1
        ind = 1 # IMPTYPE_IND[i] - 1
        if(debug): print ' improper dihedral index ',i, ATYPE[a_i-1],ATYPE[a_j-1],ATYPE[a_k-1],ATYPE[a_l-1],IMPTYPE_F[i]
        #if(debug): print ' dihedral index ',ind, a_i,a_j,a_k,a_l,ATYPE[a_i-1],ATYPE[a_j-1],ATYPE[a_k-1],ATYPE[a_l-1],IMPTYPE_F[i]
        # func_type = IMPTYPE_F[ind]
        # Hack 
        func_type = IMPTYPE_F[i]
        F.write(  '%10d %10d %10d %10d %5d \n' % (a_i,a_j,a_k,a_l,func_type) )
    F.write( '\n' )
    # Print tail
    if(debug): sys.exit("debug imps ")
    
    F.write(' [ system ]  \n ')
    F.write(' Molecular \n')
    F.write('\n')
    F.write('  [ molecules ]  \n ')
    F.write(' MOL    1  \n')
    F.write('\n')
    # close top file 
    F.close()

    return


def print_itp(ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES):
    import sys,top
    #
    # ' checking connections '
    #
    F = open( 'ff-new.itp','w')
    F.write(';  new ff parameters \n')
    F.write(' \n ')
    #
    # Check atom types 
    #
    debug = 0
    AT_LIST = []
    F.write(' [ atomtypes ] \n')
    for i in range( len(ATYPE) ):
        
        check = 1
        for ff_i in range (len(FF_ATOMTYPES)):
            #print ATYPE[i] , FF_ATOMTYPES[ff_i]
            AT_i = ATYPE[i]
            FF_i = FF_ATOMTYPES[ff_i] #.split()
            if( debug):
                print  AT_i , FF_i[0]
            if ( AT_i == FF_i[0] ):
                #
                # check if new atom type
                #
                at_new = 1
                for at_indx in range( len(AT_LIST) ):
                    if( AT_LIST[at_indx] == AT_i ):
                        at_new = 0
                if( at_new ):
                    AT_LIST.append( AT_i )
                    ff_line = '%s \n' % ' '.join(map(str, FF_ATOMTYPES[ff_i]))
                    F.write(  ff_line )
                    
                check = 0
                break
        if check :
            print ' atom type ',i,ASYMB[i],' ',ATYPE[i],'  not in ff'
            print ' grep "',ATYPE[i],'" ','ff.itp'
            sys.exit(' unknow ff specification found in gromacs.print_itp ')
            
    F.write(' \n ')
    #
    # Check bonds  
    #
    NBD_LIST = []*2
    F.write(' [ bondtypes ] \n')
    for i in range( len(BONDS)):        
        check = 1
        for ff_i in range (len(FF_BONDTYPES)):
            # print BONDS[i]
            a_i = BONDS[i][0] 
            a_j = BONDS[i][1] 
            AT_i = ATYPE[ BONDS[i][0] ]
            AT_j = ATYPE[ BONDS[i][1] ]
            FF_l = FF_BONDTYPES[ff_i] #.split() 
            if ( AT_i == FF_l[0] and  AT_j == FF_l[1]   ):
                check = 0
                break
            if (  AT_i == FF_l[1] and  AT_j == FF_l[0]  ):
                check = 0
                break
            
        if check :
            #print ' bond ',i,' ',AT_i,'(',ASYMB[a_i],GTYPE[a_i],')',' - ',AT_j,'(',ASYMB[a_j],GTYPE[a_j],')','  not in ff'
            print ' bond ',i,' ',AT_i,'(',ASYMB[a_i],')',' - ',AT_j,'(',ASYMB[a_j],')','  not in ff'
            print ' grep "',AT_i,'" ','ff.itp',' | grep "',AT_j
            print ' index ',a_i,' or index ',a_j
            sys.exit(' unknow ff specification ')
        else:
            # print AT_i,AT_j,'  good '
            #
            # check if new bond type
            #
            new = 1
            for indx in range( len(NBD_LIST) ):
                if( NBD_LIST[indx][0] == FF_l[0] and NBD_LIST[indx][1] == FF_l[1]   ):
                    new = 0
            if( new ):
                NBD_LIST.append( [ FF_l[0] ,FF_l[1]  ] )
                ff_line = '%s  \n' % ' '.join(map(str, FF_BONDTYPES[ff_i] ))
                F.write(  ff_line )
    F.write(' \n ')
    #
    # Check angles 
    #
    ANG_LIST = []*3
    F.write(' [ angletypes ] \n')
    for i in range( len(ANGLES)):        
        check = 1
        for ff_i in range (len(FF_ANGLETYPES)):
            # print BONDS[i]
            AT_i = ATYPE[ ANGLES[i][0] ]
            AT_j = ATYPE[ ANGLES[i][1] ]
            AT_k = ATYPE[ ANGLES[i][2] ]
            FF_l = FF_ANGLETYPES[ff_i] #.split() 
            if ( AT_i == FF_l[0] and  AT_j == FF_l[1] and  AT_k == FF_l[2]   ):
                check = 0
                break
            if ( AT_k == FF_l[0] and  AT_j == FF_l[1] and  AT_i == FF_l[2]   ):
                check = 0
                break
            
        if check :
            print ' angle ',i,' ',AT_i,' - ',AT_j,' - ',AT_k,'  not in ff'
            print ' grep "',AT_i,'" ','ff.itp',' | grep "',AT_j ,'" | grep "',AT_k,'"'
            print ' index ',ANGLES[i][0] -1 ,' or index ',ANGLES[i][1]-1,' or index ',ANGLES[i][2]-1
            sys.exit(' unknow ff specification ')
        else:
            # print AT_i,AT_j,AT_k,'  good '
            #
            # check if new bond type
            #
            new = 1
            for indx in range( len(ANG_LIST) ):
                if( ANG_LIST[indx][0] == FF_l[0] and ANG_LIST[indx][1] == FF_l[1] and ANG_LIST[indx][2] == FF_l[2]   ):
                    new = 0
            if( new ):
                ANG_LIST.append( [ FF_l[0] ,FF_l[1] ,FF_l[2]  ] )
                ff_line = '%s  \n' % ' '.join(map(str, FF_ANGLETYPES[ff_i] ))
                F.write(  ff_line )
                
    F.write(' \n ')
    #
    # Check dihedrals
    #
    debug = 0
    DIH_LIST = []*4
    F.write(' [ dihedraltypes ] \n')
    for i in range( len(DIH)):        
        check = 1
        for ff_i in range (len(FF_DIHTYPES)):
            # print BONDS[i]
            AT_i = ATYPE[ DIH[i][0] ]
            AT_j = ATYPE[ DIH[i][1] ]
            AT_k = ATYPE[ DIH[i][2] ]
            AT_l = ATYPE[ DIH[i][3] ]
            FF_l = FF_DIHTYPES[ff_i] #.split()

            if ( AT_i == FF_l[0] and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  AT_l == FF_l[3]   ):
                check = 0
                break
            if ( AT_l == FF_l[0] and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  AT_i == FF_l[3]   ):
                check = 0
                break
            
        
        if check :
            #print " searching for alternatives X-i-j-X "
                
            for ff_i in range (len(FF_DIHTYPES)):
                # print BONDS[i]
                AT_i = ATYPE[ DIH[i][0] ]
                AT_j = ATYPE[ DIH[i][1] ]
                AT_k = ATYPE[ DIH[i][2] ]
                AT_l = ATYPE[ DIH[i][3] ]
                FF_l = FF_DIHTYPES[ff_i] #.split()
    
                if ( FF_l[0] == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
                if ( FF_l[0] == 'X' and  AT_k == FF_l[1] and  AT_j == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
                if ( FF_l[0] == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == AT_l   ):
                    check = 0
                    break
                if ( FF_l[0] == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == AT_i   ):
                    check = 0
                    break
                if ( FF_l[0] == AT_i and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
                if ( FF_l[0] == AT_l and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
            
        if check :
            print ' dih ',i,' ',AT_i,' - ',AT_j,' - ',AT_k,' - ',AT_l,'  not in ff'
            # print  FF_DIHTYP[ff_i]
            sys.exit(' unknow ff specification ')
        else:
            # print AT_i,AT_j,AT_k,AT_l,'  good '
            #
            # check if new dih type
            #
            new = 1
            for indx in range( len(DIH_LIST) ):
                if( DIH_LIST[indx][0] == FF_l[0] and DIH_LIST[indx][1] == FF_l[1] and DIH_LIST[indx][2] == FF_l[2] and DIH_LIST[indx][3] == FF_l[3]   ):
                    new = 0
                    
            if( new ):
                DIH_LIST.append( [ FF_l[0] ,FF_l[1] ,FF_l[2] ,FF_l[3] ] )
                
                #
                #atom_i =  DIH[i][1]  
                #atom_j =  DIH[i][2] 
                #NNAB_i = top.calc_nnab(atom_i,NBLIST,NBINDEX) - 1
                #NNAB_j = top.calc_nnab(atom_j,NBLIST,NBINDEX) - 1
                #n_mod = float( NNAB_i + NNAB_j)
                #FF_DIHTYPES_mod = []
                #if(debug):
                #    print ''
                #    print ' pre mod ',FF_DIHTYPES[ff_i]
                #    print n_mod 
                #
                #mod_val = -1
                #ind_max = 12
                #if( debug ): print ' dih type ',FF_DIHTYPES[ff_i][4]
                #if( int(FF_DIHTYPES[ff_i][4].strip()) == 9 ): ind_max = 6
                #FF_DIHTYPES_mod =[]
                #for indx_dihff in range (len(FF_DIHTYPES[ff_i])):
                #    if( debug): print ' index ',indx_dihff,ind_max
                #    if ( indx_dihff > 4 and mod_val < 0 ):
                #        mod_val = 1
                #    if ( FF_DIHTYPES[ff_i][indx_dihff] == ';' or indx_dihff > ind_max ):
                #        mod_val = 0
                #    if( mod_val == 1 ):
                #        if(debug):
                #            print indx_dihff
                #            print " val to mod ", FF_DIHTYPES[ff_i][indx_dihff] 
                #        FF_DIHTYPES_mod.append( str( " %16.8f " % (float( FF_DIHTYPES[ff_i][indx_dihff]  ) /n_mod)))
                #    else:
                #        FF_DIHTYPES_mod.append(  FF_DIHTYPES[ff_i][indx_dihff] ) 
                #         
                #if(debug):  print ' post mod ',FF_DIHTYPES_mod 
                
                #ff_line = '%s  ' % ' '.join(map(str, FF_DIHTYPES[ff_i][0:5] ) )
                #ff_line = '%s  \n' % ' '.join(map(str, FF_DIHTYPES_mod ) )
                ff_line =  '%s ' % ' '.join(map(str, FF_DIHTYPES[ff_i] ) )
                #ff_line = '\n %s  %s  %s  %s  %s ' % ,AT_i,AT_j,AT_k,AT_l
                F.write(  '\n %s  ;  %s - %s - %s - %s ' % (ff_line,AT_i,AT_j,AT_k,AT_l  ) )

    #
    # Check improper dihedrals
    #
    debug = 0
    IMPS_LIST = []*4
    for i in range( len(IMPS)):
        check = 1
        for ff_i in range (len(FF_DIHTYPES)):
            # print BONDS[i]
            AT_i = ATYPE[ IMPS[i][0] ]
            AT_j = ATYPE[ IMPS[i][1] ]
            AT_k = ATYPE[ IMPS[i][2] ]
            AT_l = ATYPE[ IMPS[i][3] ]
            FF_l = FF_DIHTYPES[ff_i] #.split()

            if ( AT_i == FF_l[0] and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  AT_l == FF_l[3]   ):
                check = 0
                break
            if ( AT_l == FF_l[0] and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  AT_i == FF_l[3]   ):
                check = 0
                break
            
        
        if check :
            for ff_i in range (len(FF_DIHTYPES)):
                # print BONDS[i]
                AT_i = ATYPE[ IMPS[i][0] ]
                AT_j = ATYPE[ IMPS[i][1] ]
                AT_k = ATYPE[ IMPS[i][2] ]
                AT_l = ATYPE[ IMPS[i][3] ]
                FF_l = FF_DIHTYPES[ff_i] #.split()
    
                if ( FF_l[0] == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
                if ( FF_l[0] == 'X' and  AT_k == FF_l[1] and  AT_j == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
                if ( FF_l[0] == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == AT_l   ):
                    check = 0
                    break
                if ( FF_l[0] == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == AT_i   ):
                    check = 0
                    break
                if ( FF_l[0] == AT_i and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
                if ( FF_l[0] == AT_l and  AT_j == FF_l[1] and  AT_k == FF_l[2] and FF_l[3] == 'X'   ):
                    check = 0
                    break
            
        if check :
            print ' dih ',i,' ',AT_i,' - ',AT_j,' - ',AT_k,' - ',AT_l,'  not in ff'
            # print  FF_DIHTYP[ff_i]
            sys.exit(' unknow ff specification ')
        else:
            # print AT_i,AT_j,AT_k,AT_l,'  good '
            #
            # check if new dih type
            #
            new = 1
            for indx in range( len(IMPS_LIST) ):
                if( IMPS_LIST[indx][0] == FF_l[0] and IMPS_LIST[indx][1] == FF_l[1] and IMPS_LIST[indx][2] == FF_l[2] and IMPS_LIST[indx][3] == FF_l[3]   ):
                    new = 0
                    
            if( new ):
                IMPS_LIST.append( [ FF_l[0] ,FF_l[1] ,FF_l[2] ,FF_l[3] ] )
                
                #
                #atom_i =  DIH[i][1]  
                #atom_j =  DIH[i][2] 
                #NNAB_i = top.calc_nnab(atom_i,NBLIST,NBINDEX) - 1
                #NNAB_j = top.calc_nnab(atom_j,NBLIST,NBINDEX) - 1
                #n_mod = float( NNAB_i + NNAB_j)
                #FF_DIHTYPES_mod = []
                #if(debug):
                #    print ''
                #    print ' pre mod ',FF_DIHTYPES[ff_i]
                #    print n_mod 
                #
                #mod_val = -1
                #ind_max = 12
                #if( debug ): print ' dih type ',FF_DIHTYPES[ff_i][4]
                #if( int(FF_DIHTYPES[ff_i][4].strip()) == 9 ): ind_max = 6
                #FF_DIHTYPES_mod =[]
                #for indx_dihff in range (len(FF_DIHTYPES[ff_i])):
                #    if( debug): print ' index ',indx_dihff,ind_max
                #    if ( indx_dihff > 4 and mod_val < 0 ):
                #        mod_val = 1
                #    if ( FF_DIHTYPES[ff_i][indx_dihff] == ';' or indx_dihff > ind_max ):
                #        mod_val = 0
                #    if( mod_val == 1 ):
                #        if(debug):
                #            print indx_dihff
                #            print " val to mod ", FF_DIHTYPES[ff_i][indx_dihff] 
                #        FF_DIHTYPES_mod.append( str( " %16.8f " % (float( FF_DIHTYPES[ff_i][indx_dihff]  ) /n_mod)))
                #    else:
                #        FF_DIHTYPES_mod.append(  FF_DIHTYPES[ff_i][indx_dihff] ) 
                #         
                #if(debug):  print ' post mod ',FF_DIHTYPES_mod 
                
                #ff_line = '%s  ' % ' '.join(map(str, FF_DIHTYPES[ff_i][0:5] ) )
                #ff_line = '%s  \n' % ' '.join(map(str, FF_DIHTYPES_mod ) )
                ff_line =  '%s ' % ' '.join(map(str, FF_DIHTYPES[ff_i] ) )
                #ff_line = '\n %s  %s  %s  %s  %s ' % ,AT_i,AT_j,AT_k,AT_l
                F.write(  '\n %s  ;  %s - %s - %s - %s ' % (ff_line,AT_i,AT_j,AT_k,AT_l  ) )
                                                
                                                
    F.write(  '\n ')
    F.close()
    debug = 0
    if(debug): sys.exit('print_itp')

    return  ( AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST )


def print_min_nopbcs(g_mdp):

    min_lines =' title               =  fixed dih calc '
    min_lines = min_lines + '\n'  + 'cpp                 =  /usr/bin/cpp '
    min_lines = min_lines + '\n'  + 'constraints         =  none '
    min_lines = min_lines + '\n'  + 'integrator          =  steep'
    min_lines = min_lines + '\n'  + 'pbc                 =  no'
    min_lines = min_lines + '\n'  + 'periodic_molecules  =  no '
    min_lines = min_lines + '\n'  + 'nsteps              =  1000000 ; total 0.5ns '
    min_lines = min_lines + '\n'  + 'nstlist             =  10'
    min_lines = min_lines + '\n'  + 'ns_type             =  grid'
    min_lines = min_lines + '\n'  + ';Particle-mesh Ewald'
    min_lines = min_lines + '\n'  + 'coulombtype         =  cutoff'
    min_lines = min_lines + '\n'  + 'rvdw                =  2.0'
    min_lines = min_lines + '\n'  + 'rlist               =  2.0'
    min_lines = min_lines + '\n'  + 'rcoulomb            =  2.0'
    min_lines = min_lines + '\n'  + 'emtol  = 0.001'
    min_lines = min_lines + '\n'  + 'emstep = 0.0001 '
    
    f = open(g_mdp,'w')
    f.write(min_lines)
    f.close()

    return 
    

def run_gromacs(g_gro,g_top,g_mdp,s_suf ,options):
    import sys, os
    import file_io
    from string import replace
    

    global kJtoeV,Kjtokcal

    debug = 0

    os.system(options.load_gromacs)

    run_id = replace(g_mdp,'.mdp','')
    g_tpr =  run_id+'.tpr'
    g_log =  run_id+'.log'

    # Determine if calc needs to be ran
    recalc_gro = 0
    run_calc = 1  
    if( file_io.file_exists( g_log )  ):
        print ' log file ',g_log,' already exists will use existing for energies '
        potential_energy = get_logenergy(run_id)
        if( potential_energy != 0.0 ):
            run_calc = 0
    if( recalc_gro ):
        run_calc = 1

    if( run_calc ):
        
        # Remove files so not to rerun previous 
        g_clean = 'rm ' + g_tpr
        os.system(g_clean)

        # Compile gromacs thingy
        g_gromp = options.gromacs_dir+'grompp'+s_suf +' -f '+g_mdp+' -c '+g_gro+' -p '+g_top+' -o '+run_id
        if(debug): print  g_gromp
        os.system(g_gromp)
        # Check to make sure compiled
        # run gromacs 
        g_mdrun = options.gromacs_dir+'mdrun'+s_suf +' -s '+g_tpr+' -o '+run_id+' -x '+run_id+' -c '+run_id+' -g '+run_id +' -e '+run_id
        if(debug): print  g_mdrun
        os.system(g_mdrun)
        
        potential_energy = get_logenergy(run_id)
        
    if(debug):
        print  potential_energy
        sys.exit("run_gromacs")
        
    return potential_energy
    

def get_logenergy(run_id):
    import sys, os

    global kJtoeV,Kjtokcal

    g_log =  run_id+'.log'

    # Open log file and get energy 
    f = open(g_log,'r')
    Lines = f.readlines()
    f.close()

    potential_energy = 0.0
    for line in Lines:
        col = line.split()
        if ( len(col) > 2 ):
            if( col[0] == 'Potential' and col[1] == 'Energy' and col[2] == '=' ):
                potential_energy = float( col[3])*Kjtokcal
    
    return potential_energy

def get_coord(run_id,options):
    import sys, os
    
    R = []

    global kJtoeV,Kjtokcal

    g_gro =  run_id+'.gro'
    g_tpr =  run_id+'.tpr'
    g_gro_w =  run_id+'-W.gro'
    
    # make molecules whole
    if( options.pbcs ): #0] and options.pbcs[1] and options.pbcs[2] ):
        g_mkwhole = "echo -e \" 0 \\n \" | "+options.gromacs_dir + "trjconv"+options.gromacs_sufix +" -f " + g_gro + " -s " + g_tpr + " -o  " +  g_gro_w + " -pbc whole "
        print g_mkwhole
        os.system(g_mkwhole)
        
        g_gro = g_gro_w
    
    # Open log file and get energy 
    f = open(g_gro,'r')
    Lines = f.readlines()
    f.close()
    
    line_cnt = 0
    for line in Lines :
        line_cnt = line_cnt + 1
        if( line_cnt > 2 and len(line) >= 44): # skip header
            x = line[22:28]
            y = line[29:36]
            z = line[37:45]
            R.append( [float(x)*10,float(y)*10,float(z)*10])
            
    return R
    