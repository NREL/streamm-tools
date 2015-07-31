#! /usr/bin/env python

# Dependency modules 
import numpy as np 
import sys
import copy

# Streamm toolkit modules 
from structureContainer import StructureContainer
from particles     import Particle, ParticleContainer

from parameters    import ParameterContainer
from parameters    import ljtype,   LJtypesContainer
from parameters    import bondtype, BondtypesContainer
from parameters    import angletype,AngletypesContainer
from parameters    import dihtype,  DihtypesContainer
from parameters    import imptype,  ImptypesContainer

from bonds         import Bond,     BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer
from impropers     import Improper, ImproperContainer


from periodictable import periodictable

import units

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("--in_com", dest="in_com", type="string", default="", help="Input gaussain .com file ")
    parser.add_option("--in_itp", dest="in_itp", type="string", default="", help="Input  .itp  GROMACS file ")
    parser.add_option("--out_gro", dest="out_gro", type="string", default="out.gro", help="Output .gro GROMACS file")
    parser.add_option("--out_top", dest="out_top", type="string", default="out.top", help="Output .top GROMACS file")
    parser.add_option("--out_itp", dest="out_itp", type="string", default="out.itp", help="Output .itp GROMACS file")
    #
    (options, args) = parser.parse_args()
        
    return options, args

def check_int(var):
    """
    check if variable var is an integer
    """
    
    try:
        x = int(var)
    except ValueError:
        return False
    
    return True
    
def read_com(strucC,data_file):
    """
    Read in structure information from gaussian com file

    Args:
        strucC (StructureContainer) 
        data_file (str) gaussian com file
    Return:
        strucC (StructureContainer) 
        
    """

    # Check to see if a previous read has occured
    pt_update = False
    #if( len(strucC.ptclC) > 0 ):
    #    pt_update = True

    F = open(data_file,'r')
    Lines = F.readlines()
    F.close()

    line_cnt = 0

    read_r = False
    
    for line in Lines :
        col = line.split()
        line_cnt += 1

        if( read_r ):
            if (  len(col) < 4  ):
                read_r = False
            else:
                p_i += 1
                symbol = str(col[0])
                r_i = [float(col[1]),float(col[2]),float(col[3]) ]
                
                if( pt_update ):
                    strucC.ptclC[p_i].position = r_i
                    strucC.ptclC[p_i].setTagsDict({"symbol":symbol})                    
                else:
                    pt_i = Particle( r_i ) 
                    pt_i.setTagsDict({"symbol":symbol})                    
                    strucC.ptclC.put(pt_i)

        if( len(col) == 2 ):
            if( check_int( col[0] ) and check_int( col[1] )  ):
                read_r = True
                p_i = 0
        
    return (strucC)

def write_gro(strucC,data_file):
    """
    Write gromacs structure file 
    """
    #
    latvec = strucC.getLatVec()

    gro_lines =  " com2gro \n"
    gro_lines += " %-2i \n"  %( int(len(strucC.ptclC)) )
    atom_indx = 0 
    for pid, pt_i  in strucC.ptclC:
        atom_indx += 1
        if( atom_indx > 10000): atom_indx = 1
        r_i = pt_i.position
        r_i_nm = [units.convert_angstroms_nm(r_i[0]) ,units.convert_angstroms_nm(r_i[1]) ,units.convert_angstroms_nm(r_i[2]) ]
        gro_lines += "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"  % ( pt_i.tagsDict["residue"],pt_i.tagsDict["resname"][:5],pt_i.tagsDict["gtype"][:5],atom_indx,r_i_nm[0],r_i_nm[1],r_i_nm[2] )
    gro_lines += " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (units.convert_angstroms_nm(latvec[0][0]),units.convert_angstroms_nm(latvec[1][1]),units.convert_angstroms_nm(latvec[2][2]),units.convert_angstroms_nm(latvec[0][1]),units.convert_angstroms_nm(latvec[0][2]),units.convert_angstroms_nm(latvec[1][0]),units.convert_angstroms_nm(latvec[1][2]),units.convert_angstroms_nm(latvec[2][0]),units.convert_angstroms_nm(latvec[2][1])) 

    F = open( data_file, 'w' )
    F.write(gro_lines)
    F.close()

def delta_r_c(r_i,r_j,latticevec):
    """
    Find magnitude of dr_ij using cubic (_c) periodic boundry conditions 

    Args:
         r_i (vector numpy[3] ) position of particle i 
         r_j (vector numpy[3] ) position of particle j
         latticevec (numpy[3][3] ) 3 lattice vectors 
    Returns:
         dr_pbc (numpy[3]) vector from r_i to r_j 
    """
    r_ij  = r_j - r_i
    
    
    r_x = r_ij[0] - latticevec[0][0] * round( r_ij[0]/  latticevec[0][0] )
    r_y = r_ij[1] - latticevec[1][1] * round( r_ij[1]/  latticevec[1][1] )
    r_z = r_ij[2] - latticevec[2][2] * round( r_ij[2]/  latticevec[2][2] )
    
    dr_pbc = np.array( [r_x,r_y,r_z] )

    return dr_pbc


def getAngle(r_i,r_j):
    """
    Calcuate angle
      k - i - j 
      r_i = r_ik
      r_j = r_ij
      cos( \theta ) = ( a dot b ) / ( |a| |b| )
      
    Args:
         r_i (numpy[3] ) vector from particle i to particle k
         r_j (numpy[3] ) vector from particle i to particle j 
    Returns:
         ang_deg (float) angle_{kij}
      
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


def build_covnablist(strucC):
    """
    Build covalent neighbor list from elements and positions 
    """

    debug = False 
    p_time = False
    
    n_ptcl = len( strucC.ptclC )
    maxnnab = n_ptcl*12            # Assume max nieghbors as fcc 

    cov_buffer = 1.25
    
    radi_cov =  []
    cov_nblist = np.empty( maxnnab,  dtype=int )
    cov_nbindx = np.empty( maxnnab,  dtype=int )

    NNAB = 0

    if( p_time ): t_i = datetime.datetime.now()

    if( len( strucC.ptclC) <= 0 ):
        error_line = " empyt particle container passed to  build_covnablist "
        sys.exit(error_line)

    for pid_i, ptclObj_i  in strucC.ptclC:
        cov_nbindx[pid_i] = NNAB + 1
        el_symb = ptclObj_i.tagsDict["symbol"]
        el_i = ptclObj_i.tagsDict["number"]
        r_i = np.array( [ float( ptclObj_i.position[0] ),float( ptclObj_i.position[1] ),float( ptclObj_i.position[2] )] )
	rc_i = ptclObj_i.tagsDict["cov_radii"]*cov_buffer
	
	NNAB_i = NNAB 
        for pid_j, ptclObj_j  in strucC.ptclC:
            if( pid_j != pid_i ):
                el_j =  ptclObj_j.tagsDict["number"]
                r_j = np.array( [ float( ptclObj_j.position[0] ),float( ptclObj_j.position[1] ),float( ptclObj_j.position[2] )] )
		r_ij = r_j - r_i
		mag_dr =  np.linalg.norm(r_ij)
                #r_ij = delta_r(r_i,r_j)
                rc_j = ptclObj_j.tagsDict["cov_radii"]*cov_buffer
                r_cov = rc_i + rc_j
    
                if( mag_dr <= r_cov ):
                    NNAB = NNAB + 1
                    cov_nblist[NNAB] =  pid_j	
		    if( debug ):
			print '    atom i/j ', pid_i,pid_j,el_i,el_j
			print '       cov radi ',rc_i , rc_j
			print '       r_i ',r_i
			print '       r_j ',r_j
			print '       r_ij ',r_ij 
			print '       |r_ij| ',mag_dr 
			print '       r_cov ',r_cov
		    
	if( debug ):
	    print "  atom ",pid_i ,el_i, " has ",NNAB - NNAB_i  ," bonded nieghbors "
	
    if( p_time ):
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	print "  build_nablist dt ",dt_min,dt_sec

    
    # Account for final atom position
    cov_nbindx[pid_i+1] =  NNAB + 1

    if ( debug ):
        for pid_i, ptclObj_i  in strucC.ptclC:
            N_o = cov_nbindx[ pid_i  ]
            N_f = cov_nbindx[ pid_i + 1 ] - 1
            NNAB = N_f - N_o + 1
            print ' atom ', pid_i ,' has ',NNAB," neighbors ",N_o    
            #
            # Find number of elements
            #
            for indx in range( N_o,N_f+1):
                j = cov_nblist[indx]
                print  "     j  ",j
                #    el_j = ELN[j]
                #    ELCNT[j] = ELCNT[j] + 1

        sys.exit('debug build_covnablist')

    return (cov_nblist, cov_nbindx)


def nblist_bonds(strucC,cov_nblist, cov_nbindx):
    import sys
    """
    Generate bonds from neighbor list
    """
    debug = False

    if( debug):
        print " Initial # of bonds ",len(strucC.bondC)
    
    for pid_i, ptclObj_i  in strucC.ptclC:    
        N_o = cov_nbindx[pid_i ]
        N_f = cov_nbindx[pid_i+ 1 ] - 1
        for indx_j in range( N_o,N_f+1):
            pid_j = cov_nblist[indx_j]
            if( pid_j > pid_i):
                b_i = Bond( pid_i, pid_j )            
                strucC.bondC.put(b_i)

    if( debug):
        print " Final # of bonds ",len(strucC.bondC)
        sys.exit(" debug nblist_bonds")


def nblist_angles(strucC,cov_nblist, cov_nbindx):
    import sys
    """
    Generate angles from neighbor list 
    """

    for pid_i, ptclObj_i  in strucC.ptclC:    
        N_o = cov_nbindx[pid_i ]
        N_f = cov_nbindx[pid_i+ 1 ] - 1
        NNAB = N_f - N_o + 1
        if( NNAB >= 2 ):
            for indx_j in range( N_o,N_f):
                pid_j = cov_nblist[indx_j]
                for indx_k in range( indx_j+1,N_f+1):
                    pid_k = cov_nblist[indx_k]
                    if ( pid_j != pid_k ):
                        a_i = Angle( pid_k,pid_i, pid_j )            
                        strucC.angleC.put(a_i)

def nblist_dih(strucC,cov_nblist, cov_nbindx):
    import sys
    """
    Generate dihedrals from neighbor list 
    """

    limdih = False 
    limitdih_n = 0
 
    for pid_i, ptclObj_i  in strucC.ptclC:    
        N_o = cov_nbindx[pid_i ]
        N_f = cov_nbindx[pid_i+ 1 ] - 1
        NNAB = N_f - N_o + 1
        for indx_j in range( N_o,N_f+1):
            pid_j = cov_nblist[indx_j]
            if( pid_j > pid_i):
		dih_ij_cnt = 0
		atom_k_i = -1
		atom_l_i = -1
                for indx_k in range( N_o,N_f+1 ):
                    pid_k = cov_nblist[indx_k]
                    if ( pid_k != pid_j ):
                        No_j =  cov_nbindx[ pid_j  ]
                        Nf_j = cov_nbindx[pid_j+ 1 ] - 1
                        for indx_l in range( No_j,Nf_j+1):
                            pid_l = cov_nblist[indx_l]
                            if ( pid_l != pid_i and pid_l != pid_k ):
				if( limdih ):
				    if( dih_ij_cnt < limitdih_n ):
					if(  limitdih_n ==  2 and atom_k != atom_k_i and atom_l != atom_l_i ):

                                            d_i = Dihedral( pid_k, pid_i, pid_j, pid_l )            
                                            strucC.dihC.put(d_i)
                                            dih_ij_cnt += 1
					    atom_k_i = pid_k
					    atom_l_i = pid_l
					elif( limitdih_n !=  2 ):
                                            d_i = Dihedral( pid_k, pid_i, pid_j, pid_l )            
                                            strucC.dihC.put(d_i)
					    dih_ij_cnt += 1
					
				else:
                                    d_i = Dihedral( pid_k, pid_i, pid_j, pid_l )            
                                    strucC.dihC.put(d_i)

def calc_nnab(i,NBINDEX):
    """
    Return number of nieghbors for a given atom 
    """
    #
    # Find number of elements 
    #
    N_o = NBINDEX[ i  ]
    N_f = NBINDEX[  i+1  ] - 1 
    NNAB = N_f - N_o + 1
    return NNAB


def calc_elcnt(i,strucC,NBLIST,NBINDEX):
    """
    Return 
    """
    
    import numpy
    #
    # Find number of elements 
    #
    ELCNT = numpy.zeros(120, dtype =int )
    N_o = NBINDEX[ i  ]
    N_f = NBINDEX[  i+1  ] - 1 
    for indx in range( N_o,N_f+1):
        j = NBLIST[indx]
        el_j = int( strucC.ptclC[j].tagsDict["number"] )
	if( el_j >= 0 ):
	    ELCNT[el_j] = ELCNT[el_j] + 1

    return ELCNT 



def add_ff_prop(strucC):
    """
     Add in particle properties need for .gro file
     """
    for pid, pt_i  in strucC.ptclC:
        add_dict = pt_i.tagsDict
        add_dict["residue"] = 1
        add_dict["resname"] = "MOLRES"
        add_dict["fftype"] = "??"
        add_dict["lmptype"] = -1
        add_dict["qgroup"] = 1
        add_dict["ffmass"] = pt_i.tagsDict["mass"] 
        add_dict["chain"] = "??"
        add_dict["gtype"] = pt_i.tagsDict["symbol"]
        pt_i.setTagsDict(add_dict)                    
        
    return (strucC)


def add_atomic_prop(strucC):
    """
     Add in particle properties based on atomic symbol
     """

    # Load periodic table 
    elements = periodictable()
        
    for pid, pt_i  in strucC.ptclC:
        try:
            symbol_i = pt_i.tagsDict["symbol"]
        except KeyError:
            error_line = "Atomic symbol not set"
            sys.exit(error_line)

        el_i = elements.getelementWithSymbol(symbol_i)
        add_dict = pt_i.tagsDict
        add_dict["number"] = el_i.number
        add_dict["mass"] = el_i.mass
        add_dict["cov_radii"] = el_i.cov_radii
        add_dict["vdw_radii"] = el_i.vdw_radii
        pt_i.setTagsDict(add_dict)                    
        
    return (strucC)


def read_itp(parmC,ff_file):
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
    # Get DEFAULTS
    # 
    read_DEFAULTS = False
    for line in Lines:
        col = line.split()
        if ( len(col) > 1 ):
            if( read_DEFAULTS ):
                if ( len(col) >= 2 ):
                    if ( col[0][:1] != ';' and  col[0][:1] != '[' ):
                        col_int = 0
                            
                        parmC.set_nbfunc( int( col[col_int] ) )
                        parmC.set_combmixrule( int( col[col_int+1] ) )
                        parmC.set_genpairs( str( col[col_int+2] ) )
                        parmC.set_fudgeLJ( float( col[col_int+3] ) )
                        parmC.set_fudgeQQ( float( col[col_int+4] ) )
                        # Need to check with ljmixrule found in top
                        #  if one is found...
                if ( col[0][:1] == '[' ):
                    read_DEFAULTS = False
                    
            if ( col[0][:1] == '[' and col[1] == 'defaults' ):
                read_DEFAULTS = True

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

                        col_int = 0 
                        if(  len(col[0]) > 2 ):
                            col_int = 1
                        
                        ptype1 = col[col_int]
                        ljtyp_i = ljtype(ptype1)
                        # Set parameters according to type
                        mass_i = float( col[col_int+2] )
                        g_sig = float(col[col_int+5])
                        g_ep  = float(col[col_int+6])
                        if(  parmC.get_combmixrule() == 2 or parmC.get_combmixrule() == 3 ):
                            sigma = units.convert_nm_angstroms(g_sig)
                            epsilon = units.convert_kJmol_kcalmol(g_ep)
                        else:
                            print "uknown mixing rule for combmixrule "
                            sys.exit(" error ")

                        ljtyp_i.setmass(mass_i)
                        ljtyp_i.setparam(epsilon,sigma)
                        parmC.ljtypC.put(ljtyp_i)
                        # print btyp_i
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
                    if ( col[0][:1] != ';' and col[0][:1] != '#' ):
                        ptype1 = col[0]
                        ptype2 = col[1]
                        ptype3 = col[2]
                        ptype4 = col[3]
                        # Set parameters according to type
                        gfunc_type = int( col[4] )                    # Gromacs id

                        if( gfunc_type == 1 or gfunc_type == 9 ):
                            theat_s = float( col[5] )
                            kb = units.convert_kJmol_kcalmol( float( col[6] ) )
                            mult = float( col[7] )
                            dtype = "multiharmonic"
                            # Vd(theta) = kb[1 + cos(mult theta - theat_s)]


                        elif( gfunc_type == 2 ):
                            e0 = float( col[5] )
                            ke = units.convert_g_angle_kb( float( col[6] ) )
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
                            
                        elif( gfunc_type == 4 ):
                            e0 = float( col[5] )
                            ke = units.convert_kJmol_kcalmol( float( col[6] ) )
                            pn = int( col[7] )
                            dtype = "periodicimproper"
                        else:
                            print " unknow dihedral type ",gfunc_type
                            raise TypeError

                        if( gfunc_type == 2 ):
                            imptyp_i = imptype(ptype1,ptype2,ptype3,ptype4,dtype)
                            imptyp_i.set_g_indx(gfunc_type)
                            imptyp_i.setimp(e0,ke)
                            parmC.imptypC.put(imptyp_i)
                        elif( gfunc_type == 4 ):
                            imptyp_i = imptype(ptype1,ptype2,ptype3,ptype4,dtype)
                            imptyp_i.set_g_indx(gfunc_type)
                            imptyp_i.setimp(e0,ke)
                            imptyp_i.set_pn(pn)
                            parmC.imptypC.put(imptyp_i)
                        elif( gfunc_type == 1 or  gfunc_type == 9  ):
                            dtyp_i = dihtype(ptype1,ptype2,ptype3,ptype4,dtype)
                            dtyp_i.set_g_indx(gfunc_type)
                            dtyp_i.setharmonic( mult, kb,theat_s)
                            parmC.dtypC.put(dtyp_i)
                        elif( gfunc_type == 3 ):
                            dtyp_i = dihtype(ptype1,ptype2,ptype3,ptype4,dtype)
                            dtyp_i.set_g_indx(gfunc_type)
                            dtyp_i.setrb(C0,C1,C2,C3,C4,C5)  # Sets oplsa as well since they are equivalent 
                            parmC.dtypC.put(dtyp_i)
                        elif(  gfunc_type == 5  ):
                            dtyp_i = dihtype(ptype1,ptype2,ptype3,ptype4,dtype)
                            dtyp_i.set_g_indx(gfunc_type)
                            dtyp_i.setopls(k1,k2,k3,k4)   # Sets rb as well since they are equivalent 
                            parmC.dtypC.put(dtyp_i)
                        else:
                            print " tyring to set unknow dihedral type ",gfunc_type
                            raise TypeError

                        #print btyp_i
    debug = False
    if( debug):
        print " Dihedrals have been read in "
        for  did, dtypObj in parmC.dtypC:
            print dtypObj
        sys.exit(" itp dih readin debug ")

    return parmC


def oplsaa_atomtypes(strucC,cov_nblist, cov_nbindx):
    """
    Guess OPLS-aa atom types based on coordination and nearest nieghbors 
    """
    
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = calc_nnab(pid_i,cov_nbindx)
        ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        #
        # label carbons 
        #
        if ptclObj_i.tagsDict["number"] == 6 :
            if int(NNAB) == 4 :
                ptclObj_i.tagsDict["fftype"] = 'CT' # Alkane
            if  int(NNAB) == 3 :
                ptclObj_i.tagsDict["fftype"] = 'CA'  # Conjugated 
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'C:'   # Allene

                    
            if int(NNAB) == 1 :
                ptclObj_i.tagsDict["fftype"] = '' # Aromatic C
                error_line =  " WARNING!!! carbon index ",pid_i," bonded to single atom "
                sys.exit(error_line)
        #
        # label oxygens
        #
        if( ptclObj_i.tagsDict["number"] == 8 ):
            if int(NNAB) == 1 :
                ptclObj_i.tagsDict["fftype"] = 'O' # double bonded
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'OS' # ether

        #
        # label nitrogens 
        #
        if ptclObj_i.tagsDict["number"] == 7 :
            if int(NNAB) == 3 :      # amide
                ptclObj_i.tagsDict["fftype"] = 'N' 

        #
        # label sulfurs
        #
        if( ptclObj_i.tagsDict["number"] == 16 ):
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'S'   #  Thioether RSR (UA)


    #
    # label hydrogens
    #
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = calc_nnab(pid_i,cov_nbindx)
        ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        if( ptclObj_i.tagsDict["number"] == 1 ):
            if ( NNAB > 1  ):
                sys.exit(' over coordinated H')
            pid_j = cov_nblist[ cov_nbindx[pid_i] ]
            ptclObj_j = strucC.ptclC[pid_j]
            ELCNT_j =  calc_elcnt(pid_j,strucC,cov_nblist,cov_nbindx)

            if ( ptclObj_j.tagsDict["fftype"]== 'CA' ):
                ptclObj_i.tagsDict["fftype"] = 'HA' #
            if ( ptclObj_j.tagsDict["fftype"]== 'CT' ):
                ptclObj_i.tagsDict["fftype"] = 'HC' #


    # Check for unidentified atoms
    for pid_i, ptclObj_i  in strucC.ptclC:
        if ( ptclObj_i.tagsDict["fftype"] == '?' ):
            print ' atom ', pid_i , ptclObj_i.tagsDict["number"],' unknow '
            sys.exit(' Unknow atom ')



def set_param(struc_o,param_all,norm_dihparam,cov_nblist, cov_nbindx):
    """
    Set force-field parameters
    """

    use_last = True 

    log_file = "param.log"
    param_out = open(log_file,"w")

    ptclC_o =  struc_o.ptclC
    bondC_o  = struc_o.bondC
    angleC_o  = struc_o.angleC
    dihC_o  = struc_o.dihC
    impC_o  = struc_o.impC
    
    ljtypC_all =  param_all.ljtypC
    btypC_all =  param_all.btypC
    atypC_all =  param_all.atypC
    dtypC_all =  param_all.dtypC
    imptypC_all =  param_all.imptypC
    
    #
    # Create new parameter container to hold each unique type
    #   which exsists in the structure container struc_i
    #   this is necessary for parameter outputs to no have
    #   redundent values 
    #
    paramC_p = ParameterContainer()
    
    ljtypC_p =  paramC_p.ljtypC
    btypC_p =  paramC_p.btypC
    atypC_p =  paramC_p.atypC
    dtypC_p =  paramC_p.dtypC
    imptypC_p =  paramC_p.imptypC
    paramC_p.set_nbfunc( param_all.get_nbfunc() )
    paramC_p.set_combmixrule( param_all.get_combmixrule() )
    paramC_p.set_genpairs( param_all.get_genpairs() )
    paramC_p.set_fudgeLJ( param_all.get_fudgeLJ() )
    paramC_p.set_fudgeQQ( param_all.get_fudgeQQ() )
    #
    # Count atom types
    #
    debug_lj = False 
    for pid_o, ptclObj_o  in ptclC_o:    
        new_type = True
        lj_p = 0
        fftype_i = ptclObj_o.tagsDict["fftype"]
        #mass_i = ptclObj_o.mass

        if( debug_lj ):
            print " Checking type ",fftype_i
        
        for lj_p, ljObj_p  in ljtypC_p:    
            if( fftype_i == ljObj_p.ptype1 ):
                new_type = False
                ptclObj_o.type = str(lj_p)
                ptclObj_o.tagsDict["lmptype"] = lj_p
                if(debug_lj):
                    print ' maches previous ',pid_o, lj_p,ljObj_p.ptype1
                 
        if( new_type ):

            if( debug_lj ):
                print " New type ",fftype_i
                    
            # Find type in ljtypC_all
            cnt_check = 0
            type_found = False 
            ptclObj_o.type = str(lj_p+1)
            ptclObj_o.tagsDict["lmptype"] = lj_p+1
            for lj_all, ljObj_all  in ljtypC_all:
                if( fftype_i == ljObj_all.ptype1):
                    cnt_check += 1
                    type_found = True 
                    set_obj = ljObj_all
                if( type_found and not use_last  ):
                    # Set mass of particle type
                    #   This is need for some force-field input files (LAMMPS)
                    #mass_i =
                    set_obj.setpid( ptclObj_o.tagsDict["number"] )
                    #ljObj_all.setmass( mass_i )
                    ljtypC_p.put(set_obj)
                    type_found = False
                    
            if(  use_last and type_found ):
                set_obj.setpid( ptclObj_o.tagsDict["number"] )
                #ljObj_all.setmass( mass_i )
                ljtypC_p.put(set_obj)
                type_found = False 
                
                    
            if( cnt_check < 1 ):
                print " No LJ parameters were found for atom type %s ",fftype_i
                #raise TypeErrorb
            elif( cnt_check > 1 ):
                print " Multiple LJ parameters (%d) were found for atom type %s "%(cnt_check,fftype_i)
                if( not use_last  ):
                    raise TypeError
    if( debug_lj ):

        for lj_all, ljObj_all  in ljtypC_p:
            print ljObj_all
        sys.exit("LJ debug 89798")
        
    #
    # Count bonds types
    #
    debug = 0
    for b_o, bondObj_o  in bondC_o:
        new_type = True
        btyp_p = 0
        pid_i = bondObj_o.pgid1
        pid_j = bondObj_o.pgid2
        fftype_i =  ptclC_o[ pid_i ].tagsDict["fftype"]
        fftype_j =  ptclC_o[ pid_j ].tagsDict["fftype"]
        r_i = np.array( ptclC_o[ bondObj_o.pgid1 ].position  )
        r_j = np.array( ptclC_o[ bondObj_o.pgid2 ].position  )
        bond_len = np.linalg.norm(delta_r_c(r_i,r_j,struc_o.getLatVec() ))
        
        for btyp_p, btypObj_p  in btypC_p:
            p_i = btypObj_p.ptype1 
            p_j = btypObj_p.ptype2
            if( fftype_i == p_i  and  fftype_j == p_j ):
                new_type = False
                bondObj_o.set_lmpindx(btyp_p)
                bondObj_o.set_g_indx(btypObj_p.get_g_indx())
                log_line=" Setting bond atoms %s - %s numbers %d - %d wiht bond length %f to type %d with r_o %f  delta %f \n"%(fftype_i,fftype_j,pid_i,pid_j,bond_len,btyp_p,btypObj_p.get_r0(),bond_len-btypObj_p.get_r0() )
                param_out.write(log_line)
                break 
            elif( fftype_i == p_j  and  fftype_j == p_i ):
                new_type = False
                bondObj_o.set_lmpindx(btyp_p)
                bondObj_o.set_g_indx(btypObj_p.get_g_indx())
                log_line=" Setting bond atoms %s - %s numbers %d - %d wiht bond length %f to type %d with r_o %f  delta %f \n"%(fftype_i,fftype_j,pid_i,pid_j,bond_len,btyp_p,btypObj_p.get_r0(),bond_len-btypObj_p.get_r0() )
                param_out.write(log_line)
                break 
                
        if( new_type ):
            # Find type in btypC_all
            cnt_check = 0
            type_found = False 
            bondObj_o.set_lmpindx(btyp_p+1)
            
            for b_all, btypObj_all  in btypC_all:
                all_i = btypObj_all.ptype1 
                all_j = btypObj_all.ptype2
                if( fftype_i == all_i  and  fftype_j == all_j ):
                    type_found = True 
                if( fftype_j == all_i  and  fftype_i == all_j ):
                    type_found = True
                    
                if( type_found ):
                    cnt_check += 1
                    bondObj_o.set_g_indx(btypObj_all.get_g_indx())
                    btypC_p.put(btypObj_all)
                    if( debug ):
                        print " %d  Bond parameters were found for bond type %s-%s "%(cnt_check,fftype_i,fftype_j)
                        
                    type_found = False 
                    log_line=" Setting bond atoms %s - %s numbers %d - %d wiht bond length %f to type %d with r_o %f  delta %f \n"%(fftype_i,fftype_j,pid_i,pid_j,bond_len,btyp_p,btypObj_all.get_r0(),bond_len-btypObj_all.get_r0() )
                    param_out.write(log_line)
                    
            if( cnt_check < 1 ):
                print " No Bond parameters were found for bond type %s-%s "%(fftype_i,fftype_j)
                raise TypeError
            elif( cnt_check > 1 ):
                print " %d  Bond parameters were found for bond type %s-%s "%(cnt_check,fftype_i,fftype_j)
                for btyp_p, btypObj_p  in btypC_p:
                    print btyp_p ,btypObj_p.ptype1 ,btypObj_p.ptype2                    
                if( not use_last  ):
                    raise TypeError
                    

    #
    # Count angles types
    #
    debug = 0
    for a_o,angleObj_o in angleC_o:
        new_type = True
        atyp_p = 0
        pid_k = angleObj_o.pgid1
        pid_i = angleObj_o.pgid2
        pid_j = angleObj_o.pgid3
        fftype_k =  ptclC_o[ angleObj_o.pgid1 ].tagsDict["fftype"]
        fftype_i =  ptclC_o[ angleObj_o.pgid2 ].tagsDict["fftype"]
        fftype_j =  ptclC_o[ angleObj_o.pgid3 ].tagsDict["fftype"]
        r_k = np.array( ptclC_o[ pid_k ].position  )
        r_i = np.array( ptclC_o[ pid_i ].position  )
        r_j = np.array( ptclC_o[ pid_j ].position  )
        r_ik = delta_r_c(r_i,r_k,struc_o.getLatVec() )
        r_ij = delta_r_c(r_i,r_j,struc_o.getLatVec() )
        angle_kij = getAngle(r_ik,r_ij)
        for atyp_p, atypObj_p  in atypC_p:
            p_k = atypObj_p.ptype1 
            p_i = atypObj_p.ptype2 
            p_j = atypObj_p.ptype3
            theta0_kij = atypObj_p.get_theta0()
            delta_theta = angle_kij - theta0_kij
            if( fftype_k == p_k  and  fftype_i == p_i  and  fftype_j == p_j ):
                new_type = False
                angleObj_o.set_lmpindx(atyp_p)
                angleObj_o.set_g_indx(atypObj_p.get_g_indx())
                
                log_line=" Setting angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f to type %d with theta_o %f  delta %f \n"%(fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij,atyp_p,theta0_kij,delta_theta )
                param_out.write(log_line)
                    
                break 
            if( fftype_j == p_k  and  fftype_i == p_i  and  fftype_k == p_j ):
                new_type = False
                angleObj_o.set_lmpindx(atyp_p)
                angleObj_o.set_g_indx(atypObj_p.get_g_indx())
                log_line=" Setting angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f to type %d with theta_o %f  delta %f \n"%(fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij,atyp_p,theta0_kij,delta_theta )
                param_out.write(log_line)
                break 
                
        if( new_type ):
            # Find type in btypC_all
            cnt_check = 0
            type_found = False 
            angleObj_o.set_lmpindx(atyp_p+1)
            for a_all, atypObj_all  in atypC_all:
                all_k = atypObj_all.ptype1 
                all_i = atypObj_all.ptype2 
                all_j = atypObj_all.ptype3
                theta0_kij =  atypObj_all.get_theta0()
                delta_theta = angle_kij-theta0_kij
                if( fftype_k == all_k  and fftype_i == all_i  and  fftype_j == all_j ):
                    type_found = True 
                if( fftype_j == all_k  and  fftype_i == all_i  and  fftype_k == all_j ):
                    type_found = True 

                if( type_found ):
                    cnt_check += 1
                    angleObj_o.set_g_indx(atypObj_all.get_g_indx())
                    atypC_p.put(atypObj_all)
                    type_found = False
                    if( debug ):
                        print " %d Angles parameters were found for bond type %s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j)
                    log_line=" Setting angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f to type %d with theta_o %f  delta %f \n"%(fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij,atyp_p+1,theta0_kij,delta_theta )
                    param_out.write(log_line)
                    
            if( cnt_check < 1 ):
                print " No Angles parameters were found for bond type %s-%s-%s "%(fftype_k,fftype_i,fftype_j)
                raise TypeError
            elif( cnt_check > 1 ):
                log_line=" %d Angles parameters were found for angle atoms %s - %s - %s numbers %d - %d - %d  wiht angle %f  \n"%(cnt_check,fftype_k,fftype_i,fftype_j,pid_k,pid_i,pid_j,angle_kij )
                #param_out.write(log_line)
                print log_line
                atypC_p.findtype(fftype_k,fftype_i,fftype_j)

                if( not use_last  ):
                    raise TypeError
               
 
    #
    # Count dihedrals types
    #
    debug = False
    if( debug):
        for d_all, dtypObj_all  in dtypC_all:
            all_k = dtypObj_all.ptype1 
            all_i = dtypObj_all.ptype2 
            all_j = dtypObj_all.ptype3
            all_l = dtypObj_all.ptype4
            print " all types in parameters ",all_k,all_i,all_j,all_l,dtypObj_all.get_type()
            #sys.exit("  type check debug ")
            
    imp_cnt = 0 
    for d_o,dihObj_o in dihC_o:
        new_type = True
        dtyp_p = 0
        pid_i = dihObj_o.pgid2 
        pid_j = dihObj_o.pgid3
        fftype_k =  ptclC_o[ dihObj_o.pgid1 ].tagsDict["fftype"]
        fftype_i =  ptclC_o[ pid_i ].tagsDict["fftype"]
        fftype_j =  ptclC_o[pid_j ].tagsDict["fftype"]
        fftype_l =  ptclC_o[ dihObj_o.pgid4 ].tagsDict["fftype"]

        if( debug):
            print " checking ",fftype_k, fftype_i,  fftype_j , fftype_l
        # Check to see if dihedral type is already in parameter set for the structure container
        for dtyp_p, dtypObj_p  in dtypC_p:
            p_k = dtypObj_p.ptype1 
            p_i = dtypObj_p.ptype2 
            p_j = dtypObj_p.ptype3
            p_l = dtypObj_p.ptype4
            if( fftype_k == p_k  and  fftype_i == p_i  and  fftype_j == p_j  and  fftype_l == p_l ):
                new_type = False
                dihObj_o.set_lmpindx(dtyp_p)
                dihObj_o.set_g_indx(dtypObj_p.get_g_indx())

                if( debug):
                    print " dihObj_o.get_lmpindx()  ",dihObj_o.get_lmpindx() 
                    print " dihObj_o.get_g_indx()  ",dihObj_o.get_g_indx() 
                    print "  previous type ",dtyp_p,p_k,p_i,p_j,p_l,dihObj_o.get_g_indx()
                break
            if( fftype_l == p_k  and  fftype_j == p_i  and  fftype_i == p_j  and  fftype_k == p_l ):
                new_type = False
                dihObj_o.set_lmpindx(dtyp_p)
                dihObj_o.set_g_indx(dtypObj_p.get_g_indx())
                if( debug):
                    print " dihObj_o.get_lmpindx()  ",dihObj_o.get_lmpindx() 
                    print " dihObj_o.get_g_indx()  ",dihObj_o.get_g_indx() 
                    print "  previous type ",dtyp_p,p_k,p_i,p_j,p_l,dihObj_o.get_g_indx()
                break

        # If it is not in the parameter set for the struture container
        #  find it in the parameters from the reference parameter file 
        if( new_type ):
            # Find type in btypC_all
            cnt_check = 0
            type_found = False 
            # Set type to new type = last type+1
            dihObj_o.set_lmpindx(dtyp_p+1)

            if( debug):
                print "  new type checking against %d read in parameters "%len(dtypC_all)

            copy_type = False 
            for d_all, dtypObj_all  in dtypC_all:
                all_k = dtypObj_all.ptype1 
                all_i = dtypObj_all.ptype2 
                all_j = dtypObj_all.ptype3
                all_l = dtypObj_all.ptype4

                if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                    copy_type = True    
                if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l   ):
                    copy_type = True

                if( copy_type ):
                    cnt_check += 1
                    dtypObj_temp = copy.deepcopy(dtypObj_all)
                    #dtypObj_temp.set_g_indx(dtypObj_all.get_g_indx())
                    type_found = True
                    copy_type = False 
                    if( debug ):
                        print " %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                        print "     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l
                    
            if( not type_found ):
                if(debug):
                    print " checking  X - FF - FF - FF "
                copy_type = False 
                for d_all, dtypObj_all  in dtypC_all:
                    all_k = dtypObj_all.ptype1 
                    all_i = dtypObj_all.ptype2 
                    all_j = dtypObj_all.ptype3
                    all_l = dtypObj_all.ptype4

                    if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                        copy_type = True
                    if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                        copy_type = True
                    if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l  ):
                        copy_type = True
                    if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                        copy_type = True
                        
                    if( copy_type ):
                        cnt_check += 1
                        dtypObj_temp = copy.deepcopy(dtypObj_all)
                        #dtypObj_temp.set_g_indx(dtypObj_all.get_g_indx())
                        type_found = True 
                        copy_type = False 
                        if( debug ):
                            print " %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                            print "     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l
                            
            if( not type_found ):
                if(debug):
                    print " checking  X - FF - FF - X "
                copy_type = False 
                for d_all, dtypObj_all  in dtypC_all:
                    all_k = dtypObj_all.ptype1 
                    all_i = dtypObj_all.ptype2 
                    all_j = dtypObj_all.ptype3
                    all_l = dtypObj_all.ptype4

                    if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                        copy_type = True
                    if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                        copy_type = True

                    if( copy_type ):
                        cnt_check += 1
                        dtypObj_temp = copy.deepcopy(dtypObj_all)
                        #dtypObj_temp.set_g_indx(dtypObj_all.get_g_indx())
                        type_found = True 
                        copy_type = False 
                        if( debug ):
                            print " %d Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                            print "     from  type to dtypC_p  from ",all_k,all_i,all_j,all_l

            if( cnt_check < 1 ):
                print " No Dih parameters were found for dih type %s-%s-%s-%s "%(fftype_k,fftype_i,fftype_j,fftype_l)
                raise TypeError
            elif( cnt_check > 1 ):
                print " %d Dih parameters were found for dih type %s-%s-%s-%s please check parameter file  "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                print dtypObj_temp
                #dtypObj_temp_list.findtype(fftype_k,fftype_i,fftype_j,fftype_l)
                if( not use_last  ):
                    raise TypeError

            if( type_found ):
                if( debug ):

                    print " adding new type to dtypC_p  from ",dtypObj_temp.type, dtypObj_temp.ptype1,dtypObj_temp.ptype2,dtypObj_temp.ptype3,dtypObj_temp.ptype4
                
                # Set FF types to read in bond to remove X's 
                dtypObj_temp.ptype1 = fftype_k
                dtypObj_temp.ptype2 = fftype_i
                dtypObj_temp.ptype3 = fftype_j
                dtypObj_temp.ptype4 = fftype_l

                if( norm_dihparam ):
                    # normalize by number of nieghbors
                    dihen_norm = 1.0
                    if( debug):
                        print " Normalizing dihedral potential "
                        print " finding types for ",pid_i,pid_j
                    NNAB_i = calc_nnab(pid_i,cov_nbindx) - 1
                    NNAB_j = calc_nnab(pid_j,cov_nbindx) - 1
                    
                    dihen_norm = float( NNAB_i + NNAB_j)/2.0

                    if(debug): print " dihen_norm ",dihen_norm
                        
                    dtypObj_temp.normforceconstants(dihen_norm)


                dihObj_o.set_lmpindx(dtyp_p+1)
                dihObj_o.set_g_indx(dtypObj_temp.get_g_indx())
                dtypC_p.put(dtypObj_temp)                
                if( debug ):
                    print " %d Dih parameters were found for dih type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                    print " len(dtypC_p) ",len(dtypC_p) 
                    print " dtypObj_temp.get_lmpindx()  ",dtypObj_temp.get_lmpindx() 
                    print " dtypObj_temp.get_g_indx()  ",dtypObj_temp.get_g_indx() 
    #
    # Count improper dihedrals types
    #
    debug = False
    if( debug):
        for d_all, imptypObj_all  in imptypC_all:
            all_k = imptypObj_all.ptype1 
            all_i = imptypObj_all.ptype2 
            all_j = imptypObj_all.ptype3
            all_l = imptypObj_all.ptype4
            print " all types in parameters ",all_k,all_i,all_j,all_l,imptypObj_all.get_type()
            #sys.exit("  type check debug ")
            
    imp_cnt = 0 
    for imp_o,impObj_o in impC_o:
        new_type = True
        imptyp_p = 0
        pid_k = impObj_o.pgid1
        pid_i = impObj_o.pgid2 
        pid_j = impObj_o.pgid3
        pid_l = impObj_o.pgid4 
        fftype_k =  ptclC_o[ impObj_o.pgid1 ].tagsDict["fftype"]
        fftype_i =  ptclC_o[ pid_i ].tagsDict["fftype"]
        fftype_j =  ptclC_o[pid_j ].tagsDict["fftype"]
        fftype_l =  ptclC_o[ impObj_o.pgid4 ].tagsDict["fftype"]

        if( debug):
            print " checking ",fftype_k, fftype_i,  fftype_j , fftype_l
        # Check to see if impedral type is already in parameter set for the structure container
        for imptyp_p, imptypObj_p  in imptypC_p:
            p_k = imptypObj_p.ptype1 
            p_i = imptypObj_p.ptype2 
            p_j = imptypObj_p.ptype3
            p_l = imptypObj_p.ptype4
            if( fftype_k == p_k  and  fftype_i == p_i  and  fftype_j == p_j  and  fftype_l == p_l ):
                new_type = False
                impObj_o.set_lmpindx(imptyp_p)
                impObj_o.set_g_indx(imptypObj_p.get_g_indx())

                if( debug):
                    print " impObj_o.get_lmpindx()  ",impObj_o.get_lmpindx() 
                    print " impObj_o.get_g_indx()  ",impObj_o.get_g_indx() 
                    print "  previous type ",imptyp_p,p_k,p_i,p_j,p_l,impObj_o.get_g_indx()
                break
            if( fftype_l == p_k  and  fftype_j == p_i  and  fftype_i == p_j  and  fftype_k == p_l ):
                new_type = False
                impObj_o.set_lmpindx(imptyp_p)
                impObj_o.set_g_indx(imptypObj_p.get_g_indx())
                if( debug):
                    print " impObj_o.get_lmpindx()  ",impObj_o.get_lmpindx() 
                    print " impObj_o.get_g_indx()  ",impObj_o.get_g_indx() 
                    print "  previous type ",imptyp_p,p_k,p_i,p_j,p_l,impObj_o.get_g_indx()
                break

        # If it is not in the parameter set for the struture container
        #  find it in the parameters from the reference parameter file 
        if( new_type ):
            # Find type in btypC_all
            cnt_check = 0
            type_found = False 
            # Set type to new type = last type+1
            impObj_o.set_lmpindx(imptyp_p+1)

            if( debug):
                print "  new type checking against %d read in parameters "%len(imptypC_all)

            copy_type = False 
            for d_all, imptypObj_all  in imptypC_all:
                all_k = imptypObj_all.ptype1 
                all_i = imptypObj_all.ptype2 
                all_j = imptypObj_all.ptype3
                all_l = imptypObj_all.ptype4

                if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                    copy_type = True    
                if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l   ):
                    copy_type = True

                if( copy_type ):
                    cnt_check += 1
                    imptypObj_temp = copy.deepcopy(imptypObj_all)
                    #imptypObj_temp.set_g_indx(imptypObj_all.get_g_indx())
                    type_found = True
                    copy_type = False 
                    if( debug ):
                        print " %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                        print "     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l
                    
            if( not type_found ):
                if(debug):
                    print " checking  X - FF - FF - FF "
                copy_type = False 
                for d_all, imptypObj_all  in imptypC_all:
                    all_k = imptypObj_all.ptype1 
                    all_i = imptypObj_all.ptype2 
                    all_j = imptypObj_all.ptype3
                    all_l = imptypObj_all.ptype4

                    if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                        copy_type = True
                    if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                        copy_type = True
                    if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l  ):
                        copy_type = True
                    if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                        copy_type = True
                        
                    if( copy_type ):
                        cnt_check += 1
                        imptypObj_temp = copy.deepcopy(imptypObj_all)
                        #imptypObj_temp.set_g_indx(imptypObj_all.get_g_indx())
                        type_found = True 
                        copy_type = False 
                        if( debug ):
                            print " %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                            print "     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l
                            
            if( not type_found ):
                if(debug):
                    print " checking  X - FF - FF - X "
                copy_type = False 
                for d_all, imptypObj_all  in imptypC_all:
                    all_k = imptypObj_all.ptype1 
                    all_i = imptypObj_all.ptype2 
                    all_j = imptypObj_all.ptype3
                    all_l = imptypObj_all.ptype4

                    if ( all_k == 'X' and  all_i == fftype_i and  all_j == fftype_j and all_l == 'X'   ):
                        copy_type = True
                    if ( all_l == 'X' and  all_j == fftype_i and   all_i == fftype_j  and all_k == 'X'   ):
                        copy_type = True

                    if( copy_type ):
                        cnt_check += 1
                        imptypObj_temp = copy.deepcopy(imptypObj_all)
                        #imptypObj_temp.set_g_indx(imptypObj_all.get_g_indx())
                        type_found = True 
                        copy_type = False 
                        if( debug ):
                            print " %d Imp Dih parameters were found for bond type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                            print "     from  type to imptypC_p  from ",all_k,all_i,all_j,all_l

            if( cnt_check < 1 ):
                print " No Dih parameters were found for dih type %s-%s-%s-%s "%(fftype_k,fftype_i,fftype_j,fftype_l)
                raise TypeError
            elif( cnt_check > 1 ):
                print " %d Dih parameters were found for dih type %s-%s-%s-%s please check parameter file  "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                print imptypObj_temp
                #imptypObj_temp_list.findtype(fftype_k,fftype_i,fftype_j,fftype_l)
                if( not use_last  ):
                    raise TypeError

            if( type_found ):
                if( debug ):

                    print " adding new type to imptypC_p  from ",imptypObj_temp.type, imptypObj_temp.ptype1,imptypObj_temp.ptype2,imptypObj_temp.ptype3,imptypObj_temp.ptype4
                
                # Set FF types to read in bond to remove X's 
                imptypObj_temp.ptype1 = fftype_k
                imptypObj_temp.ptype2 = fftype_i
                imptypObj_temp.ptype3 = fftype_j
                imptypObj_temp.ptype4 = fftype_l

                norm_impdihparam = False 
                if( norm_impdihparam ):
                    # normalize by number of nieghbors
                    dihen_norm = 1.0
                    if( debug):
                        print " Normalizing dihedral potential "
                        print " finding types for ",pid_i,pid_j
                    NNAB_i = calc_nnab(pid_i,cov_nbindx) - 1
                    NNAB_j = calc_nnab(pid_j,cov_nbindx) - 1
                    
                    dihen_norm = float( NNAB_i + NNAB_j)/2.0

                    if(debug): print " dihen_norm ",dihen_norm
                    print " dihen_norm ",dihen_norm

                    imptypObj_temp.normforceconstants(dihen_norm)

                impObj_o.set_lmpindx(imptyp_p+1)
                impObj_o.set_g_indx(imptypObj_temp.get_g_indx())
                imptypC_p.put(imptypObj_temp)

                if( debug ):
                    print " %d Dih parameters were found for dih type %s-%s-%s-%s "%(cnt_check,fftype_k,fftype_i,fftype_j,fftype_l)
                    print " len(imptypC_p) ",len(imptypC_p) 
                    print " imptypObj_temp.get_lmpindx()  ",imptypObj_temp.get_lmpindx() 
                    print " imptypObj_temp.get_g_indx()  ",imptypObj_temp.get_g_indx() 

    debug = False 
    if(debug):
        print " LJ atom types found %d "%(len(ljtypC_p))
        for lj_p, ljObj_p  in ljtypC_p: 
            print lj_p,ljObj_p.ptype1,ljObj_p.mass,ljObj_p.epsilon,ljObj_p.sigma
        print " Bond types found %d "%(len(btypC_p))
        for btyp_p, btypObj_p  in btypC_p:
            print btyp_p ,btypObj_p.ptype1 ,btypObj_p.ptype2,btypObj_p.get_lmpindx(),btypObj_p.get_g_indx()
        print " Angle types found %d "%(len(atypC_p))
        for atyp_p, atypObj_p  in atypC_p:
            print atyp_p ,atypObj_p.ptype1 ,atypObj_p.ptype2,atypObj_p.ptype3,atypObj_p.get_lmpindx(),atypObj_p.get_g_indx()
        print " Dih types found %d "%(len(dtypC_p))
        for dtyp_p, dtypObj_p  in dtypC_p:
            print dtyp_p ,dtypObj_p.ptype1 ,dtypObj_p.ptype2,dtypObj_p.ptype3,dtypObj_p.ptype4,dtypObj_p.get_lmpindx(),dtypObj_p.get_g_indx()
        print " imp Dih types found %d "%(len(paramC_p.imptypC))
        for imptyp_p, dtypObj_p  in imptypC_p:
            print imptyp_p ,dtypObj_p.ptype1 ,dtypObj_p.ptype2,dtypObj_p.ptype3,dtypObj_p.ptype4,dtypObj_p.get_lmpindx(),dtypObj_p.get_g_indx()
        sys.exit('find_types')

    debug = False 
    if(debug):
        print "  All particles should have new type labeled as interger stored as a string "
        for pid_o, ptclObj_o  in ptclC_o:
            print ptclObj_o.tagsDict["fftype"],ptclObj_o.type
        for d_o,dihObj_o in dihC_o:
            print " lmpindx() g_indx()  ",d_o,dihObj_o.pgid1,dihObj_o.pgid2,dihObj_o.pgid3,dihObj_o.pgid4, dihObj_o.get_lmpindx() ,dihObj_o.get_g_indx() 

    param_out.close()

    return paramC_p,struc_o


def write_top(strucC,top_file,out_itp):
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
    F.write('#include "%s"  \n'%(out_itp))
    F.write('\n')
    F.write( '[ moleculetype ] \n')
    F.write( ' MOL           3 \n')
    F.write('\n')
    F.write( '[ atoms ] \n' )
    F.write( '; nr  :  fftype : residue #  : residue name  : atom name  :  charge group #:  charge : mass  \n')
    TOTAL_CHARGE = 0.0

    for pid, ptclObj  in strucC.ptclC:
        F.write( "%5d %5s %5d %10s %5s %5d %16.12f %12.6f  \n" % (pid,ptclObj.tagsDict["fftype"],ptclObj.tagsDict["residue"],ptclObj.tagsDict["resname"],ptclObj.tagsDict["gtype"],ptclObj.tagsDict["qgroup"],ptclObj.charge,ptclObj.tagsDict["ffmass"]) )
        TOTAL_CHARGE = TOTAL_CHARGE + float(ptclObj.charge)

    print ' Total charge = ',TOTAL_CHARGE
    F.write( '' + '\n' )
    #
    # write bonds
    #
    F.write( ' [ bonds ] \n' )
    F.write( '; i    j  func       b0          kb \n')    
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
    F.write( ';  i    j    k  func       th0       cth \n' )
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
    F.write(  ' [ dihedrals ] \n' )        
    F.write(  ';i   j   k   l	   func	 \n' )
    print len(strucC.dihC), ' number of dihedrals '
    for d_o,dihObj_o in strucC.dihC:
        pid_k = dihObj_o.pgid1
        pid_i = dihObj_o.pgid2
        pid_j = dihObj_o.pgid3
        pid_l = dihObj_o.pgid4
        #if(debug): print ' dihedral index ',ind, pid_k,pid_i,pid_j,pid_l,ATYPE[a_i-1],ATYPE[a_j-1],ATYPE[a_k-1],ATYPE[a_l-1]
        gfunc_type = dihObj_o.get_g_indx()
        F.write(  '%10d %10d %10d %10d %5d \n' % (pid_k,pid_i,pid_j,pid_l,gfunc_type) )

    print len(strucC.impC), ' number of improper dihedrals '
    for d_o,impObj_o in strucC.impC:
        pid_k = impObj_o.pgid1
        pid_i = impObj_o.pgid2
        pid_j = impObj_o.pgid3
        pid_l = impObj_o.pgid4
        #if(debug): print ' impedral index ',ind, pid_k,pid_i,pid_j,pid_l,ATYPE[a_i-1],ATYPE[a_j-1],ATYPE[a_k-1],ATYPE[a_l-1]
        gfunc_type = impObj_o.get_g_indx()
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


def write_itp(paramC,out_itp):
    """"
    Write gromacs parameter file
    """
    
    F = open( out_itp,'w')
    F.write(';  new ff parameters \n')
    F.write(' \n ')
    
    
    F.write(' \n [ defaults ] ')
    F.write(' \n ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ ')
    nbfunc = paramC.get_nbfunc()
    combmixrule =  paramC.get_combmixrule()
    genpairs =  paramC.get_genpairs()
    fudgeLJ =  paramC.get_fudgeLJ()
    fudgeQQ =  paramC.get_fudgeQQ()

    
    F.write(' \n  %d %d %s  %f %f ' % ( nbfunc,combmixrule,genpairs,fudgeLJ,fudgeQQ  ))
    F.write(' \n ')
    F.write(' \n ')

    ljtypC_p =  paramC.ljtypC
    btypC_p =  paramC.btypC
    atypC_p =  paramC.atypC
    dtypC_p =  paramC.dtypC
    imptypC_p =  paramC.imptypC

    #
    # Write particle types   
    #
    F.write('\n [ atomtypes ] ')
    for lj_p, ljObj_p  in ljtypC_p:
        sigma = units.convert_angstroms_nm( ljObj_p.sigma )
        epsilon = units.convert_kcalmol_kJmol( ljObj_p.epsilon )
        out_line = "\n %s   %d  %f %f %s %f %f  "%(ljObj_p.ptype1,ljObj_p.pid,ljObj_p.get_mass(),ljObj_p.charge,ljObj_p.ptype,sigma,epsilon)
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
            error_line = "No impropers should be in dihedral type container "
            sys.exit(error_line)
        elif( g_type == 3  ):
            Clist_kcalmol  = dtypObj_p.get_rbClist()
            Clist = []
            for Cindx in Clist_kcalmol:
                Clist.append( units.convert_kcalmol_kJmol(Cindx))
            out_line = "\n %s   %s   %s   %s  %d  %f  %f  %f  %f  %f %f "%(dtypObj_p.get_ptype1(),dtypObj_p.get_ptype2(),dtypObj_p.get_ptype3(),dtypObj_p.get_ptype4(),g_type,Clist[0],Clist[1],Clist[2],Clist[3],Clist[4],Clist[5])            
        else:
            print "  unknown gromacs dihedral type index  ",g_type
            sys.exit(" error in printing itp file ")
        F.write(out_line)

    for d_p,imptypObj_p in imptypC_p:
        g_type = imptypObj_p.get_g_indx()
        if( g_type == 2  ):
            e0, ke_kcalmol  = imptypObj_p.getimp()
            ke = units.convert_kb_g_angle( ke_kcalmol )
            out_line = "\n %s   %s   %s   %s  %d   %f  %f  "%(imptypObj_p.get_ptype1(),imptypObj_p.get_ptype2(),imptypObj_p.get_ptype3(),imptypObj_p.get_ptype4(),g_type,e0, ke)            
        elif( g_type == 3  ):
            e0, ke_kcalmol  = imptypObj_p.getimp()
            pn  = imptypObj_p.get_pn()
            ke = units.convert_kb_g_angle( ke_kcalmol )
            out_line = "\n %s   %s   %s   %s  %d   %f  %f  %d "%(imptypObj_p.get_ptype1(),imptypObj_p.get_ptype2(),imptypObj_p.get_ptype3(),imptypObj_p.get_ptype4(),g_type,e0, ke, np)            
        else:
            print "  unknown gromacs improper dihedral type index  ",g_type
            sys.exit(" error in printing itp file ")
        F.write(out_line)


    F.write(' \n ')


    F.close()
    
def main():
    """
    Read in gaussian fchk file and create an xyz file 
    """

    #        
    # Read options 
    #
    options, args = get_options()
    #
    #  Initialize blank system 
    # 
    struc_o = StructureContainer()
    #struc_o.unset_verbose()
    param_o = ParameterContainer()

    struc_o = read_com(struc_o,options.in_com)
    struc_o = add_atomic_prop(struc_o)
    struc_o = add_ff_prop(struc_o)

    cov_nblist, cov_nbindx = build_covnablist(struc_o)
    nblist_bonds(struc_o,cov_nblist, cov_nbindx)
    nblist_angles(struc_o,cov_nblist, cov_nbindx)
    nblist_dih(struc_o,cov_nblist, cov_nbindx)

    oplsaa_atomtypes(struc_o,cov_nblist, cov_nbindx)

    param_o = read_itp(param_o,options.in_itp)
    norm_dihparam = False
    param_i,struc_o  = set_param(struc_o,param_o,norm_dihparam,cov_nblist, cov_nbindx)
    write_gro(struc_o,options.out_gro)
    write_top(struc_o,options.out_top,options.out_itp)
    write_itp(param_i,options.out_itp)

    
if __name__=="__main__":
    main()
