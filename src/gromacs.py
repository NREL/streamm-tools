#! /usr/bin/env python
"""
subroutines for gromacs file manipulation
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# Email travis.kemper@nrel.gov
# Version 4.00 


from particles  import Particle
from bonds      import Bond,     BondContainer
from angles     import Angle,    AngleContainer
from dihedrals  import Dihedral, DihedralContainer

from parameters import ParameterContainer
from parameters import BondtypesContainer
from parameters import AngletypesContainer
from parameters import DihtypesContainer
from parameters import ljtype
from parameters import bondtype
from parameters import angletype
from parameters import dihtype


def read_gro(strucC,in_gro):
    """
    Read gromacs structure file

    Arguments:
        struc_o (StructureContainer)
        in_gro  (str) GROMACS .gro file

    Returns:
        struc_o (StructureContainer)
    """

    import sys, numpy 

    # atomicpy functions
    import file_io

    debug = False 

    if( file_io.file_exists(in_gro )):
        if( debug ):
            print "Reading in ",in_gro
    else:
        print " Specified .gro file ",in_gro," does not exisit "
        sys.exit("Invalid file ")

    # Open .gro file and read in lines 
    F = open(in_gro , 'r' )
    Lines = F.readlines()
    F.close()

    # Check to see if a previous read has occured
    pt_overwrite = False
    if( len(strucC.ptclC) > 0 ):
        pt_overwrite = True
    # Check of conistent number of atoms
    n_pt = int( Lines[1])
    if( pt_overwrite ):
        if(  len(strucC.ptclC) + 1 != n_pt):
            sys.exit(" Inconsistent number of atoms " )
    #
    # Read in .gro file
    #
    line_cnt = 0
    ptcl_cnt = 0 
    for line in Lines :
        line_cnt = line_cnt + 1
        if( line_cnt > 2 and len(line) >= 44 and ptcl_cnt < n_pt): # skip header
            # Set particle i 
            ptcl_cnt += 1 
            #
            g = line[10:15]
            x = line[20:28]
            y = line[28:36]
            z = line[36:44]
            # Convert from nm to angstroms 
            #r_i =  numpy.array(  [float(x)*10,float(y)*10,float(z)*10] )
            r_i =    [float(x)*10,float(y)*10,float(z)*10]
            if(debug):
                print " particle ",ptcl_cnt,g,r_i
            if( pt_overwrite ):
                pt_i = strucC.ptclC[ptcl_cnt]
                pt_i.tagsDict["gtype"] =  g
                pt_i.position = r_i 
            else:
                pt_i = Particle( r_i ) 
                tagsD = {"gtype":g}
                pt_i.setTagsDict(tagsD)
                strucC.ptclC.put(pt_i)


    #        
    # Get lattice vector from last line
    #
    line = Lines[-1]
    col = line.split()
    n_vec = int( len(col))
    if( n_vec == 3 ):
        strucC.latvec[0][0] = float( col[0] )*10
        strucC.latvec[1][1] = float( col[1] )*10
        strucC.latvec[2][2] = float( col[2] )*10
    if( n_vec == 9 ):
        strucC.latvec[0][0] = float( col[0] )*10
        strucC.latvec[1][1] = float( col[1] )*10
        strucC.latvec[2][2] = float( col[2] )*10
        strucC.latvec[0][1] = float( col[3] )*10
        strucC.latvec[0][2] = float( col[4] )*10
        strucC.latvec[1][0] = float( col[5] )*10
        strucC.latvec[1][2] = float( col[6] )*10
        strucC.latvec[2][0] = float( col[7] )*10
        strucC.latvec[2][1] = float( col[8] )*10

    if( debug ):
        print "      Box size ",strucC.latvec[0][0],strucC.latvec[1][1],strucC.latvec[2][2]," angstorms "
        
    return strucC

def read_top(strucC, top_infile):

    """
    Read in GROMACS topology file

    """
    
    # General python modules 
    import sys, numpy

    # atomicpy functions
    import file_io
    from periodictable import periodictable

    debug = False 
    verbose = False 

    # Load periodic table 
    pt = periodictable()

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

                    if( verbose ):
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
    read_DEFAULTS = 0
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

            if( read_DEFAULTS ):
                if ( len(col) >= 2 ):
                    if ( col[0][:1] != ';' and  col[0][:1] != '[' ):
                        ljmixrule = int( col[1] )
                if ( col[0][:1] == '[' ):
                    read_DEFAULTS = 0
                    
            if ( col[0][:1] == '[' and col[1] == 'moleculetype' ):
                read_MOLECULETYPE = 1
            if ( col[0][:1] == '[' and col[1] == 'molecules' ):
                read_MOLECULES = 1
            if ( col[0][:1] == '[' and col[1] == 'defaults' ):
                read_DEFAULTS = 1



    if( verbose ):
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


    if( verbose ):
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
    BOND_CNT = -1
    ANGLE_CNT = -1
    DIH_CNT = -1

    MOL_CNT = -1 


    # Check to see if a previous read has occured
    pt_overwrite = False
    if( len(strucC.ptclC) > 0 ):
        pt_overwrite = True

    # Find number of atoms in top file
    # Loop over molecules
    pt_cnt_i = 0 
    for repeat_indx in range( len(MOLECULECNT) ):
        mol_id = MOLECULECNT_ID[repeat_indx]
        
        
        # Find molecule type
        id_n = -1
        for id_indx in range( len(MOLECULETYPE) ):
            id_type = MOLECULETYPE[id_indx]
            if( id_type == mol_id ):
                id_n = id_indx
        if( id_n < 0 ):
            sys.exit(" error in .top read in")
            
        # Repeat atoms
        N_o = MOL_ATOMS_INDEX[id_n] #- 1
        N_f = MOL_ATOMS_INDEX[id_n+1] - 1
        np_mol = N_f - N_o + 1
        pt_cnt_i += np_mol*MOLECULECNT[repeat_indx]
        
    # Check of conistent number of atoms
    if( pt_overwrite ):
        pt_cnt = len(strucC.ptclC)
        if( pt_cnt  != pt_cnt_i):
            print " Current structure has %d atoms and %s has %d"%(pt_cnt,top_infile,ATOMS_CNT+1)
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



    #debug = True
    if( debug):
        print "len(MOLECULECNT) ",len(MOLECULECNT)
        #sys.exit(" mol mult debug ")

    # Loop over molecules
    for repeat_indx in range( len(MOLECULECNT) ):
        mol_id = MOLECULECNT_ID[repeat_indx]
        # Find molecule type
        id_n = -1
        if( debug):
            print "mol_id ",mol_id
            print " searching molecule types ",len(MOLECULETYPE)
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
                m_i = AMASS_l[atom_indx]
                q_i = CHARGES_l[atom_indx]
                if( pt_overwrite ):
                    pt_i = strucC.ptclC[A_CNT+1]
                    pt_i.charge = q_i
                    pt_i.mass = m_i
                    if( debug ):
                        print " updating ",A_CNT+1
                else:
                    r_i = [0.0,0.0,0.0]
                    type_i = "??"
                    pt_i = Particle( r_i,type_i,m_i,q_i ) 
                    strucC.ptclC.put(pt_i) 

                pt_i.tagsDict["gtype"] = GTYPE_l[atom_indx]
                pt_i.tagsDict["fftype"] = ATYPE_l[atom_indx]
                pt_i.tagsDict["resname"] = RESN_l[atom_indx]
                pt_i.tagsDict["residue"] = RESID_l[atom_indx]
                pt_i.tagsDict["qgroup"] = CHARN_l[atom_indx]
                pt_i.tagsDict["chain"] = MOL_CNT + 1

                if( debug ):
                    print " particle ",A_CNT+1,pt_i.tagsDict["fftype"], pt_i.tagsDict["chain"] 

                #el = pt.getelementWithSymbol(atomic_symb)
                #el.symbol,el.number,"mass":el.mass,"cov_radii":el.cov_radii,"vdw_radii":el.vdw_radi
                #pt_i.tagsDict["gtype"] = GTYPE_l[atom_indx]

                
                #
                # Reference atom # in global list
                #
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
                    BOND_CNT += 1
                    i_o = BONDS_l[bond_indx][0]
                    j_o = BONDS_l[bond_indx][1]
                    i_add = REF_N[i_o]
                    j_add = REF_N[j_o]

                    if( bonds_overwrite ):
                        bondObj = strucC.bondC[BOND_CNT + 1 ]
                        bondObj.pgid1 = i_o
                        bondObj.pgid2 = j_o
                    else:
                        b_i = Bond( i_o, j_o )            
                        strucC.bondC.put(b_i)
                    if( debug ):
                        print i_o,j_o , " -> ",i_add,j_add

            # Repeat angles
            N_o = MOL_ANGLES_INDEX[id_n] #- 1
            N_f = MOL_ANGLES_INDEX[id_n+1] - 1
            if( debug ): print " adding angles ",N_o,N_f
            if( N_f > N_o ):
                for angle_indx in range(N_o,N_f+1):
                    ANGLE_CNT += 1

                    k_o = ANGLES_l[angle_indx][0]
                    i_o = ANGLES_l[angle_indx][1]
                    j_o = ANGLES_l[angle_indx][2]
                    k_add = REF_N[k_o]
                    i_add = REF_N[i_o]
                    j_add = REF_N[j_o]

                    if( angles_overwrite ):
                        angleObj = strucC.angleC[ANGLE_CNT + 1 ]
                        angleObj.pgid1 = k_o
                        angleObj.pgid2 = i_o
                        angleObj.pgid3 = j_o
                    else:
                        angle_i = Angle( k_o,i_o, j_o )            
                        strucC.angleC.put(angle_i)

            # Repeat dihedral
            N_o = MOL_DIH_INDEX[id_n] #- 1
            N_f = MOL_DIH_INDEX[id_n+1] - 1
            if( debug ): print " adding dihedral ",N_o,N_f
            if( N_f > N_o ):
                for dih_indx in range(N_o,N_f+1):
                    DIH_CNT += 1
                    k_o = DIH_l[dih_indx][0]
                    i_o = DIH_l[dih_indx][1]
                    j_o = DIH_l[dih_indx][2]
                    l_o = DIH_l[dih_indx][3]
                    k_add = REF_N[k_o]
                    i_add = REF_N[i_o]
                    j_add = REF_N[j_o]
                    l_add = REF_N[l_o]

                    if( dih_overwrite ):
                        dObj = strucC.dihC[DIH_CNT + 1 ]
                        dObj.pgid1 = k_o
                        dObj.pgid2 = i_o
                        dObj.pgid3 = j_o
                        dObj.pgid4 = l_o
                    else:
                        dih_i = Dihedral( k_o,i_o, j_o,l_o )            
                        strucC.dihC.put(dih_i)

    if(debug): print " Total atoms ",A_CNT
    MOLPNT.append( A_CNT + 1)


    return (strucC,ljmixrule)


def read_itp( ff_file, ljmixrule):
    """
    Read a gromacs paramter file

    Args:
        ff_file (str) GROMACS itp file
    Return:
        parmC (ParameterContainer) with read in parameters 
    """

    import sys
    
    debug = False

    KJ_KCAL = 0.23901
    GRO_SIG = 5.61230943
    NM_ANG = 10.0

    parmC = ParameterContainer()

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
                        ptype1 = col[0]
                        ljtyp_i = ljtype(ptype1)
                        # Set parameters according to type 
                        g_sig = float(col[5])
                        g_ep  = float(col[6])
                        if( ljmixrule == 2 ):
                            epsilon = g_ep*GRO_SIG
                        elif( ljmixrule == 3 ):
                            epsilon = g_ep*NM_ANG
                        else:
                            print "uknown mixing rule for LJ parameters "
                            sys.exit(" error")
                        
                        parmC.ljtypC.put(ljtyp_i)
                        #print btyp_i
    if( debug):
        print " Bonds have been read in "
        for  ljid, ljtypObj in parmC.ljtypC:
            print ljtypObj
        
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
                        ptype1 = col[0]
                        ptype2 = col[1]
                        # Set parameters according to type 
                        g_type = int( col[2] )                    # Gromacs id 
                        if( g_type == 1 ):
                            r0 = float(col[3])*NM_ANG 
                            kb = float(col[4])*KJ_KCAL/NM_ANG/NM_ANG/2
                            btype = "harmonic"
                        btyp_i = bondtype(ptype1,ptype2,btype)

                        if( g_type == 1 ):
                            btyp_i.setharmonic(r0,kb)

                        parmC.btypC.put(btyp_i)
                        #print btyp_i
    if( debug):
        print " Bonds have been read in "
        for  bid, btypObj in parmC.btypC:
            print btypObj

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
                        ptype1 = col[0]
                        ptype2 = col[1]
                        ptype3 = col[2]
                        # Set parameters according to type 
                        gfunc_type = int( col[3] )                    # Gromacs id 
                        if( gfunc_type == 1 ):
                            theta0 = float( col[4] )
                            kb = float( col[5] )*KJ_KCAL/2
                            atype = "harmonic"
                        atyp_i = angletype(ptype1,ptype2,ptype3,atype)

                        if( gfunc_type == 1 ):
                            atyp_i.setharmonic(theta0,kb)

                        parmC.atypC.put(atyp_i)
                        #print btyp_i
    if( debug):
        print " Angles have been read in "
        for  aid, atypObj in parmC.atypC:
            print atypObj


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
                        ptype1 = col[0]
                        ptype2 = col[1]
                        ptype3 = col[2]
                        ptype4 = col[3]
                        # Set parameters according to type 
                        gfunc_type = int( col[4] )                    # Gromacs id

                        if( gfunc_type == 1 or gfunc_type == 4 ):
                            theat_s = float( col[5] )
                            kb = float( col[6] )*KJ_KCAL
                            mult = float( col[7] )
                            dtype = "harmonic"
                            # Vd(theta) = kb[1 + cos(mult theta - theat_s)]


                        elif( gfunc_type == 2 ):
                            e0 = float( col[5] )
                            ke = float( col[6] )*KJ_KCAL
                            dtype = "improper"
                            # V = 1/2 ke( e-e0)^2

                        elif( gfunc_type == 3 ):
                            C0 = float( col[5] )*KJ_KCAL
                            C1 = float( col[6] )*KJ_KCAL
                            C2 = float( col[7] )*KJ_KCAL
                            C3 = float( col[8] )*KJ_KCAL
                            C4 = float( col[9] )*KJ_KCAL
                            C5 = float( col[10] )*KJ_KCAL
                            dtype = "rb"
                            # Ryckaert-Bellemans function
                            # Vrb(theta) = \sum_n=0^5 C_n [ cos(theata - 180 )^n ]
                            # Translate to opls 
                            k1 = -1.0*( 2.0*C1 + 3.0*C3/2.0)
                            k2 = -1.0*( C2 + C4)
                            k3 = -0.5*C3
                            k4 = -0.25*C4

                        elif( gfunc_type == 5 ):
                            k1 = float( col[5] )*KJ_KCAL
                            k2 = float( col[6] )*KJ_KCAL
                            k3 = float( col[7] )*KJ_KCAL
                            k4 = float( col[8] )*KJ_KCAL
                            dtype = "opls"
                            # opls function
                            # Translate to Ryckaert-Bellemans function
                            C0 = k2 + 0.5*(k1+k3)
                            C1 = 0.5*(-1.0*k1+3.0*k3)
                            C2 = -1.0*k2 + 4.0*k3
                            C3 = -2.0*k3
                            C4 = -4.0*k4
                            C5 = 0.0
                        elif( gfunc_type == 9 ):
                            theat_s = float( col[4] )
                            kb = float( col[5] )
                            mult = float( col[6] )
                            dtype = "multiharmonic"

                        else:
                            print " unknow dihedral type ",gfunc_type
                            raise TypeError

                        dtyp_i = dihtype(ptype1,ptype2,ptype3,ptype4,dtype)

                        if( gfunc_type == 1 or  gfunc_type == 9  ):
                            d = 1.0 
                            dtyp_i.setharmonic(d, mult, kb,theat_s)

                        if( gfunc_type == 3  or gfunc_type == 5  ):
                            dtyp_i.setopls(k1,k2,k3,k4)
                            dtyp_i.setrb(C0,C1,C2,C3,C4,C5)


                        parmC.dtypC.put(dtyp_i)
                        #print btyp_i
    if( debug):
        print " Dihedrals have been read in "
        for  did, dtypObj in parmC.dtypC:
            print dtypObj

    return parmC


