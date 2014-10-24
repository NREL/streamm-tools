#! /usr/bin/env python
"""
subroutines for gromacs file manipulation
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# Email travis.kemper@nrel.gov
# Version 4.00

import sys 

import units 

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
    pt_update = False
    if( len(strucC.ptclC) > 0 ):
        pt_update = True
    # Check of conistent number of atoms
    n_pt = int( Lines[1])
    pt_cnt = len(strucC.ptclC)
    if( pt_update ):
        if( pt_cnt != n_pt):
            print " Current structure has %d atoms and %s has %d"%(pt_cnt,in_gro,n_pt)
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
            resname_i = line[5:10].strip()            
            g = line[10:15]
            x = units.convert_nm_angstroms( float( line[20:28] ))
            y = units.convert_nm_angstroms( float(line[28:36]))
            z = units.convert_nm_angstroms( float(line[36:44]))
            #r_i =  numpy.array(  [float(x)*10,float(y)*10,float(z)*10] )
            r_i =    [x,y,z]
            if(debug):
                print " particle ",ptcl_cnt,g,r_i
            if( pt_update ):
                pt_i = strucC.ptclC[ptcl_cnt]
                pt_i.tagsDict["resname"] =  resname_i
                pt_i.tagsDict["gtype"] =  g
                pt_i.position = r_i 
            else:
                pt_i = Particle( r_i ) 
                tagsD = {"gtype":g,"resname":resname_i}
                pt_i.setTagsDict(tagsD)
                strucC.ptclC.put(pt_i)


    #        
    # Get lattice vector from last line
    #
    line = Lines[-1]
    col = line.split()
    n_vec = int( len(col))
    if( n_vec == 3 ):
        strucC.latvec[0][0] = units.convert_nm_angstroms(float( col[0] ) )
        strucC.latvec[1][1] = units.convert_nm_angstroms(float( col[1] ) )
        strucC.latvec[2][2] = units.convert_nm_angstroms(float( col[2] ) )
    if( n_vec == 9 ):
        strucC.latvec[0][0] = units.convert_nm_angstroms(float( col[0] ) )
        strucC.latvec[1][1] = units.convert_nm_angstroms(float( col[1] ) )
        strucC.latvec[2][2] = units.convert_nm_angstroms(float( col[2] ) )
        strucC.latvec[0][1] = units.convert_nm_angstroms(float( col[3] ) )
        strucC.latvec[0][2] = units.convert_nm_angstroms(float( col[4] ) )
        strucC.latvec[1][0] = units.convert_nm_angstroms(float( col[5] ) )
        strucC.latvec[1][2] = units.convert_nm_angstroms(float( col[6] ) )
        strucC.latvec[2][0] = units.convert_nm_angstroms(float( col[7] ) )
        strucC.latvec[2][1] = units.convert_nm_angstroms(float( col[8] ) )

    if( debug ):
        print "      Box size ",strucC.latvec[0][0],strucC.latvec[1][1],strucC.latvec[2][2]," angstorms "
        
    return strucC

def read_top(strucC, parmC, top_infile):
    """
    Read in GROMACS topology file
    """
    
    # General python modules 
    import sys, numpy

    # atomicpy functions
    import file_io
    from periodictable import periodictable

    debug = False 
    verbose = True 

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
    itp_list = []
    for line in top_lines:
        col = line.split()
        if ( len(col) > 1 ):
            if ( col[0][:1] != ';' ):
                if( col[0] == "#include" ):
                    include_file = col[1].strip("\"")

                    F_inc = open(include_file , 'r' )
                    include_lines = F_inc.readlines()
                    F_inc.close()

                    #
                    # Compile list of itp files
                    #
                    include_sufix = include_file[-3:]
                    if( include_sufix == ".itp" ):
                        itp_list.append(include_file)
                    
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
    pt_update = False
    if( len(strucC.ptclC) > 0 ):
        if(  verbose ): #rank == 0 and
            print "Exisiting particles will updated with top information "
        pt_update = True

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
    if( pt_update ):
        pt_cnt = len(strucC.ptclC)
        if( pt_cnt  != pt_cnt_i):
            print " Current structure has %d atoms and %s has %d"%(pt_cnt,top_infile,pt_cnt_i)
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
            MOLPNT.append( A_CNT + 1 )

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
                if( pt_update ):
                    pt_i = strucC.ptclC[A_CNT+1]
                    if( debug ):
                        print " updating ",A_CNT+1
                else:
                    r_i = [0.0,0.0,0.0]
                    type_i = "??"
                    pt_i = Particle( r_i,type_i,m_i,q_i ) 
                    strucC.ptclC.put(pt_i) 

                pt_i.charge = q_i
                pt_i.mass = m_i
                pt_i.tagsDict["gtype"] = GTYPE_l[atom_indx].strip()
                pt_i.tagsDict["fftype"] = ATYPE_l[atom_indx].strip()
                pt_i.tagsDict["resname"] = RESID_l[atom_indx].strip()
                pt_i.tagsDict["residue"] = RESN_l[atom_indx]
                pt_i.tagsDict["qgroup"] = CHARN_l[atom_indx]
                pt_i.tagsDict["chain"] = MOL_CNT + 1
                
                el = pt.getelementWithMass(m_i)
                if( el.symbol == "VS" ):
                    el.symbol = ATYPE_l[atom_indx].strip()

                # HACK !!
                if(  ATYPE_l[atom_indx].strip() == "LP" ):
                    el.symbol = ATYPE_l[atom_indx].strip()

                pt_i.tagsDict["symbol"] = el.symbol
                pt_i.tagsDict["number"] = el.number
                pt_i.tagsDict["cov_radii"] = el.cov_radii
                pt_i.tagsDict["vdw_radi"] = el.vdw_radii
                #pt_i.tagsDict[""] = el.

                if( debug ):
                    print " particle ",A_CNT+1,pt_i.tagsDict["symbol"],pt_i.tagsDict["fftype"], pt_i.tagsDict["chain"] 
                
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


    #
    # Read in parameters from itp files 
    # 
    for itp_file in  itp_list:
        parmC = read_itp( parmC, ff_file, ljmixrule)

    return (strucC,parmC,ljmixrule)


def read_itp( parmC, ff_file, ljmixrule):
    """
    Read a gromacs paramter file

    Args:
        ff_file (str) GROMACS itp file
    Return:
        parmC (ParameterContainer) with read in parameters 
    """

    import sys
    
    debug = False

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
                        mass_i = float( col[2] )
                        g_sig = float(col[5])
                        g_ep  = float(col[6])
                        if( ljmixrule == 2 or ljmixrule == 3 ):
                            sigma = units.convert_nm_angstroms(g_sig)
                            epsilon = units.convert_kJmol_kcalmol(g_ep)
                        else:
                            print "uknown mixing rule for LJ parameters "
                            sys.exit(" error")

                        ljtyp_i.setmass(mass_i)
                        ljtyp_i.setparam(epsilon,sigma)
                        parmC.ljtypC.put(ljtyp_i)
                        #print btyp_i
    if( debug):
        print " LJ paramters have been read in "
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
                            r0 = units.convert_nm_angstroms( float(col[3]) )
                            kb = units.convert_g_bond_kb( float(col[4]) )
                            btype = "harmonic"
                        btyp_i = bondtype(ptype1,ptype2,btype)
                        btyp_i.set_g_indx(g_type)
                        
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
                            theta0 = float( col[4] )        # degrees 
                            kb = units.convert_g_angle_kb( float( col[5] ) )
                            atype = "harmonic"
                        atyp_i = angletype(ptype1,ptype2,ptype3,atype)
                        atyp_i.set_g_indx(gfunc_type)

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
                            kb = units.convert_kJmol_kcalmol( float( col[6] ) )
                            mult = float( col[7] )
                            dtype = "harmonic"
                            # Vd(theta) = kb[1 + cos(mult theta - theat_s)]


                        elif( gfunc_type == 2 ):
                            e0 = float( col[5] )
                            ke = units.convert_kJmol_kcalmol( float( col[6] ) )
                            dtype = "improper"
                            # V = 1/2 ke( e-e0)^2

                        elif( gfunc_type == 3 ):
                            C0 = units.convert_kJmol_kcalmol( float( col[5] ) )
                            C1 = units.convert_kJmol_kcalmol( float( col[6] ) )
                            C2 = units.convert_kJmol_kcalmol( float( col[7] ) )
                            C3 = units.convert_kJmol_kcalmol( float( col[8] ) )
                            C4 = units.convert_kJmol_kcalmol( float( col[9] ) )
                            C5 = units.convert_kJmol_kcalmol( float( col[10] ) )
                            dtype = "rb"
                            # Ryckaert-Bellemans function
                            # Vrb(theta) = \sum_n=0^5 C_n [ cos(theata - 180 )^n ]

                        elif( gfunc_type == 5 ):
                            k1 = units.convert_kJmol_kcalmol( float( col[5] ) )
                            k2 = units.convert_kJmol_kcalmol( float( col[6] ) )
                            k3 = units.convert_kJmol_kcalmol( float( col[7] ) )
                            k4 = units.convert_kJmol_kcalmol( float( col[8] ) )
                            dtype = "opls"
                            # opls function
                            
                        elif( gfunc_type == 9 ):
                            theat_s = float( col[5] )
                            kb = units.convert_kJmol_kcalmol( float( col[6] ) )
                            mult = float( col[7] )
                            dtype = "multiharmonic"

                        else:
                            print " unknow dihedral type ",gfunc_type
                            raise TypeError

                        dtyp_i = dihtype(ptype1,ptype2,ptype3,ptype4,dtype)
                        dtyp_i.set_g_indx(gfunc_type)

                        if( gfunc_type == 1 or  gfunc_type == 9  ):
                            d = 1.0 
                            dtyp_i.setharmonic(d, mult, kb,theat_s)

                        if( gfunc_type == 2 ):
                            dtyp_i.setimp(e0,ke) 

                            
                        if( gfunc_type == 3 ):
                            dtyp_i.setrb(C0,C1,C2,C3,C4,C5)  # Sets oplsa as well since they are equivalent 
                        if(  gfunc_type == 5  ):
                            dtyp_i.setopls(k1,k2,k3,k4)   # Sets rb as well since they are equivalent 

                        parmC.dtypC.put(dtyp_i)

                        #print btyp_i
    debug = False
    if( debug):
        print " Dihedrals have been read in "
        for  did, dtypObj in parmC.dtypC:
            print dtypObj
        sys.exit(" itp dih readin debug ")

    return parmC



def print_gro(strucC,gro_file):
    """
    Write gromacs structure file 
    """

    import sys 
    #
    # printing gro 
    #
    F = open( gro_file, 'w' )
    F.write( ' log to top ')
    F.write( "\n %-2i "  % int(len(strucC.ptclC)) )
    atom_indx = 0 
    for pid, ptclObj  in strucC.ptclC:
        atom_indx += 1
        if( atom_indx > 10000): atom_indx = 1
        r_i = ptclObj.position
        r_i_nm = [units.convert_angstroms_nm(r_i[0]) ,units.convert_angstroms_nm(r_i[1]) ,units.convert_angstroms_nm(r_i[2]) ]
        F.write( "\n%5d%-5s%5s%5d%8.3f%8.3f%8.3f"  % ( ptclObj.tagsDict["residue"],ptclObj.tagsDict["resname"][:5],ptclObj.tagsDict["gtype"][:5],atom_indx,r_i_nm[0],r_i_nm[1],r_i_nm[2] ))

    latvec = strucC.get_latvec()
    F.write( "\n %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n" % ( latvec[0][0]/10.0,latvec[1][1]/10.0,latvec[2][2]/10.0,latvec[0][1]/10.0,latvec[0][2]/10.0,latvec[1][0]/10.0,latvec[1][2]/10.0,latvec[2][0]/10.0,latvec[2][1]/10.0))
    F.close()

def print_top(strucC,top_file):
    """
    Write gromacs topology file
    """

    spc = str(' ')
    
    F = open( top_file, 'w' )
    #
    # write atoms 
    #
    F.write('; top file for gromacs \n ')
    F.write('\n')
    F.write('; force field file \n')
    F.write('#include "ff-new.itp"  \n')
    F.write('\n')
    F.write( '[ moleculetype ] \n')
    F.write( ' MOL           3 \n')
    F.write('\n')
    F.write( '[ atoms ] \n' )
    TOTAL_CHARGE = 0.0

    for pid, ptclObj  in strucC.ptclC:
        # print i,ATYPE[i] ,RESN[i] ,RESID[i] ,GTYPE[i],CHARN[i],CHARGES[i] ,AMASS[i]
        F.write( "%5d %5s %5d %10s %5s %5d %16.12f %12.6f  \n" % (pid,ptclObj.tagsDict["fftype"],ptclObj.tagsDict["residue"],ptclObj.tagsDict["resname"],ptclObj.tagsDict["gtype"],ptclObj.tagsDict["qgroup"],ptclObj.charge,ptclObj.mass) )
        TOTAL_CHARGE = TOTAL_CHARGE + float(ptclObj.charge)

    print ' Total charge = ',TOTAL_CHARGE
    F.write( '' + '\n' )
    #
    # write bonds
    #
    F.write( ' [ bonds ] \n' )
    print len(strucC.bondC) , ' number of bonds'
    for b_o, bondObj_o  in strucC.bondC:
        pid_i = bondObj_o.pgid1
        pid_j = bondObj_o.pgid2
        gfunc_type = bondObj_o.get_g_indx()
        F.write(  '%10d %10d %5d \n' % (pid_i,pid_j,gfunc_type) )
    F.write(  '\n' )
    #
    # write angles
    #
    F.write( ' [ angles ] \n' )
    print len(strucC.angleC)-1 , ' number of angles'
    for a_o,angleObj_o in strucC.angleC:
        pid_k = angleObj_o.pgid1
        pid_i = angleObj_o.pgid2
        pid_j = angleObj_o.pgid3
        gfunc_type = angleObj_o.get_g_indx()
        F.write(  '%10d %10d %10d %5d \n' % (pid_k,pid_i,pid_j,gfunc_type) )
    F.write( '\n' )
    #
    # write dihedrals
    #
    debug = 0
    line = str( '[ dihedrals ]' )
    F.write( line + '\n' )        
    print len(strucC.dihC), ' number of dihedrals '
    for d_o,dihObj_o in strucC.dihC:
        pid_k = dihObj_o.pgid1
        pid_i = dihObj_o.pgid2
        pid_j = dihObj_o.pgid3
        pid_l = dihObj_o.pgid4
        #if(debug): print ' dihedral index ',ind, pid_k,pid_i,pid_j,pid_l,ATYPE[a_i-1],ATYPE[a_j-1],ATYPE[a_k-1],ATYPE[a_l-1]
        gfunc_type = dihObj_o.get_g_indx()
        F.write(  '%10d %10d %10d %10d %5d \n' % (pid_k,pid_i,pid_j,pid_l,gfunc_type) )
        
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

    F.close()


def print_itp(paramC,itp_file,ff_type):
    """"
    Write gromacs parameter file
    """
    
    F = open( itp_file,'w')
    F.write(';  new ff parameters \n')
    F.write(' \n ')
    
    if( ff_type == "oplsaa" ):
            
        nbfunc = 1
        ljmixrule = 3
        genpairs = "yes"
        fudgeLJ = 0.5
        fudgeQQ = 0.5
        
    elif( ff_type == "amber" ):

        nbfunc = 1
        ljmixrule = 2
        genpairs = "yes"
        fudgeLJ = 0.5
        fudgeQQ = 0.8333
        
    else:
        print " force-field type unknown  "

    
    F.write(' \n [ defaults ] ')
    F.write(' \n ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ ')
    F.write(' \n  %d %d %s  %f %f ' % ( nbfunc,ljmixrule,genpairs,fudgeLJ,fudgeQQ  ))
    F.write(' \n ')
    F.write(' \n ')

    ljtypC_p =  paramC.ljtypC
    btypC_p =  paramC.btypC
    atypC_p =  paramC.atypC
    dtypC_p =  paramC.dtypC

    #
    # Write particle types   
    #
    F.write('\n [ atomtypes ] ')
    for lj_p, ljObj_p  in ljtypC_p:
        sigma = units.convert_angstroms_nm( ljObj_p.sigma )
        epsilon = units.convert_kcalmol_kJmol( ljObj_p.epsilon )
        out_line = "\n %s   %d  %f %f %s %f %f  "%(ljObj_p.ptype1,ljObj_p.pid,ljObj_p.mass,ljObj_p.charge,ljObj_p.ptype,sigma,epsilon)
        F.write(out_line)

    F.write(' \n ')
    #
    # Write bond types
    #
    F.write('\n [ bondtypes ] ')
    for b_p, btypeObj_p  in btypC_p:
        g_type = btypeObj_p.get_g_indx()
        if( g_type == 1 ):
            r0 = units.convert_angstroms_nm(  btypeObj_p.get_r0() )
            kb = units.convert_kb_g_bond( btypeObj_p.get_kb()) 
            out_line = "\n %s   %s  %d   %f  %f  "%(btypeObj_p.get_ptype1(),btypeObj_p.get_ptype2(),g_type,r0,kb)
        else:
            print "  unknown gromacs bond type index  ",g_type
            sys.exit(" error in printing itp file ")
            
        F.write(out_line)
    F.write(' \n ')
    #
    # Write angle  types
    #
    F.write('\n [ angletypes ] ')
    for a_p,atypeObj_p in atypC_p:        
        g_type = atypeObj_p.get_g_indx()
        if( g_type == 1  ):
            theta0 = atypeObj_p.get_theta0()
            kb = units.convert_kb_g_angle( atypeObj_p.get_kb()) 
            out_line = "\n %s   %s   %s  %d   %f  %f  "%(atypeObj_p.get_ptype1(),atypeObj_p.get_ptype2(),atypeObj_p.get_ptype3(),g_type,theta0,kb)
        else:
            print "  unknown gromacs angle type index  ",g_type
            sys.exit(" error in printing itp file ")
        F.write(out_line)
    F.write(' \n ')
    #
    # Write dihedral  types
    #
    F.write('\n [ dihedraltypes ] ')
    for d_p,dtypObj_p in dtypC_p:
        g_type = dtypObj_p.get_g_indx()
        if( g_type == 1 or g_type == 4  or g_type == 9 ):
            theat_s = dtypObj_p.get_theat_s()
            kb = units.convert_kcalmol_kJmol( dtypObj_p.get_kb() )
            mult =  dtypObj_p.get_mult()
            #print "  theat_s , kb ,mult ",theat_s , kb ,mult
            out_line = "\n %s   %s   %s   %s  %d   %f  %f  %d "%(dtypObj_p.get_ptype1(),dtypObj_p.get_ptype2(),dtypObj_p.get_ptype3(),dtypObj_p.get_ptype4(),g_type,theat_s,kb,mult)
        elif( g_type == 2  ):
            e0, ke_kcalmol  = dtypObj_p.getimp()
            ke = units.convert_kcalmol_kJmol( ke_kcalmol )
            out_line = "\n %s   %s   %s   %s  %d   %f  %f  "%(dtypObj_p.get_ptype1(),dtypObj_p.get_ptype2(),dtypObj_p.get_ptype3(),dtypObj_p.get_ptype4(),g_type,e0, ke)            
        elif( g_type == 3  ):
            Clist_kcalmol  = dtypObj_p.get_rbClist()
            Clist = []
            for Cindx in Clist_kcalmol:
                Clist.append( units.convert_kcalmol_kJmol(Cindx))
                
            out_line = "\n %s   %s   %s   %s  %d  %f  %f  %f  %f  %f "%(dtypObj_p.get_ptype1(),dtypObj_p.get_ptype2(),dtypObj_p.get_ptype3(),dtypObj_p.get_ptype4(),g_type,Cindx[0],Cindx[1],Cindx[2],Cindx[3],Cindx[4],Cindx[5])            
        else:
            print "  unknown gromacs angle type index  ",g_type
            sys.exit(" error in printing itp file ")
        F.write(out_line)

    F.write(' \n ')


    F.close()
    

        
"""

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
    
def print_itp(new_itp,norm_dihparam,ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES):
    Write gromacs parameter file
    
    import sys,top
    #
    # ' checking connections '
    #
    F = open( new_itp,'w')
    F.write(';  new ff parameters \n')
    F.write(' \n ')
    
    ff_type = "oplsaa"
    if( ff_type == "oplsaa" ):
            
        nbfunc = 1
        combrule = 3
        genpairs = "yes"
        fudgeLJ = 0.5
        fudgeQQ = 0.5
        
    elif( ff_type == "amber" ):

        nbfunc = 1
        combrule = 2
        genpairs = "yes"
        fudgeLJ = 0.5
        fudgeQQ = 0.8333
        
    else:
        print " force-field type unknown  "
        

    F.write(' \n [ defaults ] ')
    F.write(' \n ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ ')
    F.write(' \n  %d %d %s  %f %f ' % ( nbfunc,combrule,genpairs,fudgeLJ,fudgeQQ  ))
    F.write(' \n ')
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
            if ( AT_l == FF_l[0] and  AT_k == FF_l[1] and  AT_j == FF_l[2] and  AT_i == FF_l[3]   ):
                check = 0
                break
        
        if check :
            
            for ff_i in range ( len(FF_DIHTYPES) ):
                
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
            
            
        #if( not check ):

        # print AT_i,AT_j,AT_k,AT_l,'  good '
        #
        # check if new dih type
        #
        new = 1
        for indx in range( len(DIH_LIST) ):
            if( DIH_LIST[indx][0] == AT_i and DIH_LIST[indx][1] == AT_j and DIH_LIST[indx][2] == AT_k and DIH_LIST[indx][3] == AT_l   ):
                new = 0
                
        if( new ):
            DIH_LIST.append( [ AT_i,AT_j,AT_k,AT_l] )
            #
            # If dihedrals in reference itp file are for single set of atoms 
            #
            FF_DIHTYPES_mod = [ AT_i,AT_j,AT_k,AT_l ]
            #
            if( norm_dihparam ):
                atom_i =  DIH[i][1]  
                atom_j =  DIH[i][2] 
                NNAB_i = top.calc_nnab(atom_i,NBLIST,NBINDEX) - 1
                NNAB_j = top.calc_nnab(atom_j,NBLIST,NBINDEX) - 1
                n_mod = float( NNAB_i + NNAB_j )/ 2.0
                #FF_DIHTYPES_mod = []
                if(debug):
                    print ''
                    print ' pre mod ',FF_DIHTYPES[ff_i]
                    print "     nb_i nb_j",NNAB_i,NNAB_j
                    print '     n norm ',n_mod
                
                mod_val = -1
                ind_max = 12
                if( debug ): print ' dih type ',FF_DIHTYPES[ff_i][4]
                if( int(FF_DIHTYPES[ff_i][4].strip()) == 9 ): ind_max = 6
                for indx_dihff in range (4,len(FF_DIHTYPES[ff_i])):
                    if( debug): print ' index ',indx_dihff,ind_max, FF_DIHTYPES[ff_i][indx_dihff] 
                    if ( indx_dihff > 4 and mod_val < 0 ):
                        mod_val = 1
                    if ( FF_DIHTYPES[ff_i][indx_dihff] == ';' or indx_dihff > ind_max ):
                        mod_val = 0
                        
                    if( mod_val == 1 ):
                        if(debug):
                            print indx_dihff
                            print " val to mod ", FF_DIHTYPES[ff_i][indx_dihff] 
                        FF_DIHTYPES_mod.append( str( " %16.8f " % (float( FF_DIHTYPES[ff_i][indx_dihff]  ) /n_mod)))
                    else:
                        FF_DIHTYPES_mod.append(  FF_DIHTYPES[ff_i][indx_dihff] ) 
                    
                if(debug):  print ' post mod ',FF_DIHTYPES_mod
                ff_line = '  %s  ' % ' '.join(map(str, FF_DIHTYPES_mod ) )
            else:
                for indx_dihff in range (4,len(FF_DIHTYPES[ff_i])):
                    FF_DIHTYPES_mod.append(  FF_DIHTYPES[ff_i][indx_dihff] )
                    
                ff_line =  '  %s ' % ' ' .join(map(str, FF_DIHTYPES_mod ) )
            
            #ff_line = '%s  ' % ' '.join(map(str, FF_DIHTYPES[ff_i][0:5] ) )
            #ff_line =  '%s ' % ' '.join(map(str, FF_DIHTYPES[ff_i] ) )
            #ff_line = '\n %s  %s  %s  %s  %s ' % ,AT_i,AT_j,AT_k,AT_l
            F.write(  '\n %s  ;  %s - %s - %s - %s ' % (ff_line,AT_i,AT_j,AT_k,AT_l  ) )
            
    if( debug):
        sys.exit("norm ")
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

    def write_gro(self,dir_id,output_id ): # Move out of class
        Write out gromacs gro file
        # Version 1 will be dependent on Atomicpy
        import gromacs 

        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        
        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        ATYPE = []
        RING_NUMB = []
        GTYPE = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( ptclObj.position)  )
            AMASS.append( float(ptclObj.mass)  )
            CHARGES.append( float(ptclObj.charge)  )
            MOLNUMB.append( int(ptclObj.tagsDict["chain"])  )
            RESID.append( ptclObj.tagsDict["resname"][0:5]  )
            RESN.append( int(ptclObj.tagsDict["residue"])  )
            ATYPE.append( ptclObj.tagsDict["fftype"]  )
            RING_NUMB.append( int(ptclObj.tagsDict["ring"])  )
            GTYPE.append( ptclObj.tagsDict["gtype"]  )


            print " GTYPE ", ptclObj.tagsDict["gtype"] 
            print " RESID ", ptclObj.tagsDict["resname"] 
            print " RESN ", ptclObj.tagsDict["residue"] 
       
        # Set cubic lattice constant to 5 nm arbitrary 
        LV = np.zeros( (3,3) )
            
        LV[0][0] = self.latvec[0][0]
        LV[1][1] = self.latvec[1][1]
        LV[2][2] = self.latvec[2][2]
        
        out_gro = dir_id+"/"+output_id + ".gro"
        gromacs.print_gro(out_gro,GTYPE,RESID,RESN,R,LV)

    def write_top(self,dir_id,output_id,norm_dihparam,itp_file ): # Move out of class
        Write out gromacs gro file
        # Version 1 will be dependent on Atomicpy
        import gromacs , elements, top , lammps, groups , atom_types

        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        
        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        ATYPE = []
        RING_NUMB = []
        GTYPE = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( ptclObj.position)  )
            AMASS.append( float(ptclObj.mass)  )
            CHARGES.append( float(ptclObj.charge)  )
            MOLNUMB.append( int(ptclObj.tagsDict["chain"])  )
            RESID.append( ptclObj.tagsDict["resname"][0:5]  )
            RESN.append( int(ptclObj.tagsDict["residue"])  )
            ATYPE.append( ptclObj.tagsDict["fftype"]  )
            RING_NUMB.append( int(ptclObj.tagsDict["ring"])  )
            GTYPE.append( ptclObj.tagsDict["gtype"]  )
       
        BONDS = []
        for b_i,bondObj in  self.bondC:
            BONDS.append( [bondObj.pgid1 - 1, bondObj.pgid2 -1])

            print " make_top bonds ",bondObj.pgid1 , bondObj.pgid2

        # Set cubic lattice constant to 5 nm arbitrary 
        LV = np.zeros( (3,3) )
            
        LV[0][0] = self.latvec[0][0]
        LV[1][1] = self.latvec[1][1]
        LV[2][2] = self.latvec[2][2]
        
        # Find atomic number based on atomic symbol 
        ELN = elements.asymb_eln(ASYMB)
        NA = len(ELN)
        
        # Create neighbor list form bonds
        NBLIST,NBINDEX = groups.build_nablist_bonds(ELN,BONDS)
        #NBLIST,NBINDEX = self.bonded_nblist() #groups.build_nablist_bonds(ELN,BONDS)

        print_nb = True
        # Make compatable with 0-(N-1) index of atoms 
        #for n_indx in range(len(NBLIST)):
        #    if( print_nb):
        #        print " changing NBLIST ",NBLIST[n_indx] ," to ",NBLIST[n_indx] -1 
        #    NBLIST[n_indx] =NBLIST[n_indx] -1
        #for n_indx in range(len(NBINDEX)-1):
        #    if( print_nb):
        #        print " changing NBINDEX ",NBINDEX[n_indx] ," to ",NBINDEX[n_indx+1]
        #    NBINDEX[n_indx] =NBINDEX[n_indx+1]
        #
        if( print_nb):

            for p_i in range(len(self.ptclC)):
                N_i_o = NBINDEX[p_i]
                N_i_f = NBINDEX[p_i+1]
                print " atom ",p_i+1, ELN[p_i]," has ",N_i_f - N_i_o
                for indx_j in range( N_i_o,N_i_f):
                    atom_j = NBLIST[indx_j]
                    print "      nb ",atom_j+1," atomic # ", ELN[atom_j]
                    
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)

        for a_indx in range(len(ANGLES)):
            print " angles ",ANGLES[a_indx][0]+1,ANGLES[a_indx][1]+1,ANGLES[a_indx][2]+1
        
        #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)


        #
        # Set charge groups
        #
        verbose = False 
        CG_SET = []
        CHARN = []
        one = 1
        for i in range( len(ELN) ):
            CG_SET.append(one)
            CHARN.append(one)
        CHARN = top.set_chargegroups(verbose,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LV)

        # Read in parameter files 
        
        FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(itp_file)

        # Identify total number of atom types for lammps output 
        ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH)

        # Check atom types to be sure each atom of the same type has the same number of neighbors 
        ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)
    
        ATYPE_EP, ATYPE_SIG = top.atom_parameters(itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
        BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(itp_file,BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
        ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(itp_file,ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
        DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(itp_file, norm_dihparam, DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )

        IMPTYPE_F  = top.imp_parameters(itp_file)

        IMPS,IMPTYPE_F = atom_types.set_ptma_imps(NA,NBLIST, NBINDEX,ELN,ASYMB,IMPS,IMPTYPE_F)


        top_file = dir_id+"/"+output_id + ".top"
        DIH_CONST = []
        DIH_CONST_ANGLE = []
        gromacs.print_top( top_file,ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
               ,DIH_CONST,DIH_CONST_ANGLE
               ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
               ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
"""
