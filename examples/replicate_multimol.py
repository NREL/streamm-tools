#! /usr/bin/env python

import numpy as np
import random , math, sys

# Streamm toolkit modules 
import mpiBase
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

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)

    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("--mol1_data", dest="mol1_data", type="string", default="", help="Input for molecule 1 .data LAMMPS file")
    parser.add_option("--mol1_n", dest="mol1_n", type=int, default=0, help="Numnber times molecule 1 is to be replicated  ")
    parser.add_option("--mol2_data", dest="mol2_data", type="string", default="", help="Input for molecule 1 .data LAMMPS file")
    parser.add_option("--mol2_n", dest="mol2_n", type=int, default=0, help="Numnber times molecule 1 is to be replicated  ")
    parser.add_option("--mol3_data", dest="mol3_data", type="string", default="", help="Input for molecule 1 .data LAMMPS file")
    parser.add_option("--mol3_n", dest="mol3_n", type=int, default=0, help="Numnber times molecule 1 is to be replicated  ")
    parser.add_option("--out_data", dest="out_data", type="string", default="out.data", help="Output .data LAMMPS file")
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="out.xyz", help="Output xyz (xmol) file")
    #
    #
    parser.add_option("--seed", dest="seed", type="int", help="Seed for random")
    parser.add_option("--atomic_cut", dest="atomic_cut", type=float, default=3.5, help="Minimum distance between atoms of molecules ")
    parser.add_option("--lc_expand", dest="lc_expand", type=float, default=0.100, help="Fraction of the box size to increase system size after max_sys is excieded ")
    parser.add_option("--max_sys", dest="max_sys", type=float, default=3, help="Maximum system recreations at a certain lattice constant ")
    parser.add_option("--max_mol_place", dest="max_mol_place", type=float, default=50, help="Maximum attempts to place a molecule  ")
    parser.add_option("--calc_overlap", dest="calc_overlap",  default=True,action="store_false",help=" Turn off calculation of molecular overlap, handy to build a quick input file" )
    
    (options, args) = parser.parse_args()
        
    return options, args



def vec_shift(strucC,r_shift):
    """
    Shift structure by vector
    NOTE: only called in pbcs.py

    Arguments
        r_shift (numpy vector) to shift all the cordinates by

    """

    for pid, ptclObj in strucC.ptclC :
        r_i = np.array( ptclObj.position ) + r_shift
        ptclObj.position = [r_i[0],r_i[1],r_i[2]]

def center_mass(strucC):

    """
    Find center of mass of a structure
    NOTE: needs to consider units

    Return
        r_mass (numpy array) position of the center of mass
    """
    import numpy as np

    total_mass_i = strucC.getTotMass()
    r_mass = np.array( [0.0,0.0,0.0] )

    for pid, ptclObj in strucC.ptclC :
        r_mass += ptclObj.mass*np.array( ptclObj.position )

    r_mass = r_mass/total_mass_i

    return r_mass
    
def shift_center_mass(strucC,r_shift):

    """
    Translate center of mass of a structure to a location 
    NOTE: keep here (maybe re-name method) only in pbcs.py

    Return
        r_shift (numpy array) position of the center of mass

    """
    r_mass = center_mass(strucC)
    r_m_s = r_shift - r_mass
    vec_shift(strucC,r_m_s)

def sq_drij_c(r_i,r_j,latticevec):
    """
    Reutrn square magnitude of distance between two vectors
    using cubic periodic boundry conditions 
    Args:
         r_i (vector numpy[3] ) position of particle i 
         r_j (vector numpy[3] ) position of particle j
         latticevec (numpy[3][3] ) 3 lattice vectors 
    Returns:
         sq_dr (float) square of vector from r_i to r_j 

    """

    dr_pbc = delta_r_c(r_i,r_j,latticevec)
    
    sq_dr = np.dot( dr_pbc,dr_pbc)
    
    return sq_dr
    
    
def rotate(strucC, rot_angle_i, rot_angle_j):

    """
    Rotate particles in particle container. Rotation angle i around y axis
    and rotation angle j around z-axis

    Arguments
        rot_angle_i (float)  0 - pi 
        rot_angle_j (float)  0 - pi 

    """
    import numpy as np 
    import math

    # set variables of rotation matrix
    #   for rotation i around y axis 
    cy = math.cos(rot_angle_i)
    sy = math.sin(rot_angle_i)
    #   for rotation j around z axis 
    cz = math.cos(rot_angle_j)
    sz = math.sin(rot_angle_j)

    # loop over each particle 
    for pid, ptclObj in strucC.ptclC :
        xd = ptclObj.position[0]
        yd = ptclObj.position[1]
        zd = ptclObj.position[2]
        # Apply rotation matrix i and j to get new postion 
        r_x =  cy*cz*xd - sz*cy*yd + sy*zd 
        r_y =  sz*xd    + cz*yd            
        r_z = -sy*cz*xd + sy*sz*yd + cy*zd

        # ptclObj.position = numpy.array( [r_x,r_y,r_z] )
        ptclObj.position = [r_x,r_y,r_z]




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

def getlength(strucC):

    """
    Calculate length from maximum seperation between particles
    NOTE: pbcs also calls methods from StructureContainer and this should be moved.
    NOTE: pbcs functionality should be moved here since lattice vectors are here

    Return
        struc_len (float) 
    """

    sq_maxdr = -1000000.0 
    for p_i, ptclObj_i in strucC.ptclC :
        r_i = np.array( ptclObj_i.position )
        for p_j, ptclObj_j in strucC.ptclC :
            r_j = np.array( ptclObj_j.position )
            dr_ij = delta_r_c(r_i,r_j,strucC.latvec)
            r_ij_sq = np.dot( dr_ij,dr_ij)
            if( r_ij_sq > sq_maxdr):
                sq_maxdr = r_ij_sq


    struc_len = np.sqrt(sq_maxdr)

    return struc_len
            
def read_data( strucC , parmC , data_file):
    """
    Read Lammps data file


    Arguments:
        strucC     (StructureContainer)
        parmC      (ParameterContainer)
        data_file  (str) data file
    
    ReturnS:
        strucC     (StructureContainer)
        parmC      (ParameterContainer)
    
    """
        
    debug = 0
    verbose = True

    set_chain_numbers = True

    if( not set_chain_numbers ):
        print " Warning not reading in chain numbers!!! "

    # Load periodic table 
    elements = periodictable()
    
    try:
        with open(data_file,'r') as F:
            lines = F.readlines()
            F.close()
    except IOError:
        sys.exit("Invalid file ")
    
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
        print " n_impropers ",n_impropers
        print ""
        print "n_atom_types",n_atypes
        print "n_bond_types",n_btypes
        print "n_angle_types",n_angtypes
        print "n_dihedral_types",n_dtypes
        print "n_imp_dihedral_types",n_imptypes



    # Check to see if a previous read has occured
    pt_overwrite = False
    if( len(strucC.ptclC) > 0 ):
        pt_overwrite = True
    # Check of conistent number of atoms
    if( pt_overwrite ):
        if(  len(strucC.ptclC)  != n_atoms):
            print "  %d atoms in passed structure "%(len(strucC.ptclC))
            print "  %d atoms in data file "%(n_atoms)
            sys.exit(" Inconsistent number of atoms " )

    else:
        #
        # Initialize particle container
        #
        for pid_i in range(n_atoms):
            pt_i = Particle( )
            strucC.ptclC.put(pt_i)

 
    bonds_overwrite = False
    if( len(strucC.bondC) > 0 ):
        bonds_overwrite = True

        if(  len(strucC.bondC)  != n_bonds):
            print "  %d bonds in passed structure "%(len(strucC.bondC))
            print "  %d bonds in data file "%(n_bonds)
            sys.exit(" Inconsistent number of bonds " )
    
    angles_overwrite = False
    if( len(strucC.angleC) > 0 ):
        angles_overwrite = True
        if(  len(strucC.angleC)  != n_angles):
            print "  %d angles in passed structure "%(len(strucC.angleC))
            print "  %d angles in data file "%(n_angles)
            sys.exit(" Inconsistent number of angles " )
    
    dih_overwrite = False
    if( len(strucC.dihC) > 0 ):
        dih_overwrite = True
        if(  len(strucC.dihC)  != n_dihedrals):
            print "  %d dihedrals in passed structure "%(len(strucC.dihC))
            print "  %d dihedrals in data file "%(n_dihedrals)
            sys.exit(" Inconsistent number of dihedrals " )
              
    imp_overwrite = False
    if( len(strucC.impC) > 0 ):
        imp_overwrite = True
        if(  len(strucC.impC)  != n_impropers):
            print "  %d impropers in passed structure "%(len(strucC.impC))
            print "  %d impropers in data file "%(n_impropers)
            sys.exit(" Inconsistent number of impropers " )          
            
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
    ATYPE_MASS = np.zeros(n_atypes)
    ATYPE_EP = np.zeros(n_atypes)
    ATYPE_SIG = np.zeros(n_atypes)
    
    BTYPE_REF = n_btypes*[2*[""]]
    BONDTYPE_R0 = np.zeros(n_btypes)
    BONDTYPE_K = np.zeros(n_btypes)
    
    ANGTYPE_REF = n_angtypes*[3*[""]]
    ANGLETYPE_R0 = np.zeros(n_angtypes)
    ANGLETYPE_K = np.zeros(n_angtypes)
    
    DTYPE_REF = n_dtypes*[4*[""]]
    DIHTYPE_C = np.zeros((n_dtypes,4))
    DIHTYPE_F = np.zeros(n_dtypes)
    DIHTYPE_K = np.zeros(n_dtypes)
    DIHTYPE_PN = np.zeros(n_dtypes)
    DIHTYPE_PHASE = np.zeros(n_dtypes)

    IMPTYPE_REF = n_imptypes*[4*[""]]
    IMPTYPE_F = np.zeros(n_imptypes)
    IMPTYPE_E0 = np.zeros(n_imptypes)
    IMPTYPE_K = np.zeros(n_imptypes)

    MOLNUMB = n_atoms*[0]
    ATYPE_IND  = n_atoms*[0]
    CHARGES  = np.zeros(n_atoms)
    R = n_atoms*[np.zeros(3)]
    ATYPE  = n_atoms*[""]
    
    BONDS = n_bonds*[[0,0]]
    BTYPE_IND = n_bonds*[0]
    
    ANGLES = n_angles*[[0,0,0]]
    ANGTYPE_IND = n_angles*[0]
    
    DIH = n_dihedrals*[[0,0,0,0]]
    DTYPE_IND = n_dihedrals*[0]


    

    #
    # Check if values exist and need to be updated or don't and need to be created 
    #
    ljtyp_update = False
    ljtyp_cnt = 0 
    if( len(parmC.ljtypC) > 0 ):
        print " LJ types will be updated "
        ljtyp_update = True
    btyp_update = False
    btyp_cnt = 0 
    if( len(parmC.btypC) > 0 ):
        print " Bond types will be updated "
        btyp_update = True
    atyp_update = False
    atyp_cnt = 0 
    if( len(parmC.atypC) > 0 ):
        print " Angle types will be updated "
        atyp_update = True
    dtyp_update = False
    dtyp_cnt = 0 
    if( len(parmC.dtypC) > 0 ):
        print " Dihedral types will be updated "
        dtyp_update = True
    imptyp_update = False
    imptyp_cnt = 0 
    if( len(parmC.imptypC) > 0 ):
        print " Improper dihedrals types will be updated "
        imptyp_update = True
    #
    # Read in data parameters 
    #
    
    for line in lines:
        col = line.split()
        if( read_Masses and  len(col) >= 2 ):
            
            cnt_Masses += 1
            ind = int(col[0]) - 1
            ATYPE_MASS[ind] = float(col[1])
                
            if( len(col) >= 4 ):
                ATYPE_REF[ind] = col[3]
                ptype1 = col[3]
            else:
                ATYPE_REF[ind] = "??"
                ptype1 = "??"
            
            mass_i = float(col[1])

        
            if( ljtyp_update ):
                ljtyp_cnt = ind  + 1
                if( ljtyp_cnt > len(parmC.ljtypC) ):
                    print "Mass index %d larger then length of previously read ljtypC %d"%(ind,len(parmC.ljtypC))
                ljtyp_i = parmC.ljtypC[ljtyp_cnt]
                ljtyp_i.setmass(mass_i)
            else:
                ljtyp_i = ljtype(ptype1)
                ljtyp_i.setmass(mass_i)
                parmC.ljtypC.put(ljtyp_i)
            
            # Turn of mass read
            if(cnt_Masses ==  n_atypes ):
                read_Masses = 0

        if( read_Pair and  len(col) >= 3 ):
            cnt_Pair += 1
            
            ind = int(col[0]) - 1
            ATYPE_EP[ind] = float(col[1])
            ATYPE_SIG[ind] = float(col[2])

            epsilon = float(col[1])
            sigma = float(col[2])
            ljtyp_ind = int(col[0])

            ljtyp_i = parmC.ljtypC[ljtyp_ind]
            ljtyp_i.setparam(epsilon,sigma)
            # Turn pair parameter read off 
            if(cnt_Pair >=  n_atypes ):
                read_Pair = 0

        
        if( read_Bond_coeff and  len(col) >= 3 ):
            cnt_Bond_coeff += 1
            #AT_i = int(col[0])
            #AT_j = int(col[1])
            b_ind = int( col[0]) - 1
            if( b_ind > n_btypes ):
                error_line = " Error in data file index of bond parameter exceeds number of bond parameters specified with bond types "
                sys.exit(error_line)
                
            BTYPE_REF[b_ind][0] = "??"
            BTYPE_REF[b_ind][1] = "??"
            BONDTYPE_K[b_ind] = float(col[1])
            BONDTYPE_R0[b_ind] = float(col[2])

            ptype1 = "??"
            ptype2 = "??"
            lmpindx = int( col[0])
            kb = float(col[1]) 
            r0 = float(col[2]) 
            btype = "harmonic"
            g_type = 1

            if( btyp_update ):
                btyp_cnt = b_ind + 1 
                btyp_i = parmC.btypC[btyp_cnt]
                btyp_i.setharmonic(r0,kb)
                btyp_i.set_g_indx(g_type)
                btyp_i.set_lmpindx(lmpindx)
            else:
                btyp_i = bondtype(ptype1,ptype2,btype)
                btyp_i.setharmonic(r0,kb)
                btyp_i.set_g_indx(g_type)
                btyp_i.set_lmpindx(lmpindx)
                parmC.btypC.put(btyp_i)
            
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

            ptype1 = "??"
            ptype2 = "??"
            ptype3 = "??"
            lmpindx = int( col[0])
            theta0 = float( col[2] )        # degrees 
            kb =  float( col[1] ) 
            atype = "harmonic"
            gfunc_type = 1


            if( atyp_update ):
                atyp_cnt = a_ind + 1 
                atyp_i = parmC.atypC[atyp_cnt]
                atyp_i.set_g_indx(gfunc_type)
                atyp_i.set_lmpindx(lmpindx)
                atyp_i.setharmonic(theta0,kb)
            else:
                atyp_i = angletype(ptype1,ptype2,ptype3,atype)
                atyp_i.set_g_indx(gfunc_type)
                atyp_i.set_lmpindx(lmpindx)
                atyp_i.setharmonic(theta0,kb)
                parmC.atypC.put(atyp_i)
            
            if( cnt_Angle_coeff >=  n_angtypes ):
                read_Angle_coeff = 0
                
        
        if( read_Dihedral_coeff and  len(col) >= 3 ):
            cnt_Dihedral_coeff += 1
            #AT_i = int(col[0])
            #AT_j = int(col[1])
            d_ind = int( col[0]) - 1
            
            if( debug): print " reading dih type ",d_ind," cnt ",cnt_Dihedral_coeff," of ",n_dtypes
            
            if( d_ind > n_dtypes ):
                error_line =  " Error in data file index of dihedral parameter %d exceeds number of dihedral parameters %d "%(d_ind , n_dtypes)
                error_line += " specified with dihedral types "
                print sys.exit(error_line)
                
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

            ptype1 = "??"
            ptype2 = "??"
            ptype3 = "??"
            ptype4 = "??"
            # Set parameters according to type 
            gfunc_type = 3 
            dtype = "opls"

            lmpindx = int( col[0] )
            k1 =  float( col[1] )
            k2 =  float( col[2] ) 
            k3 =  float( col[3] ) 
            k4 =  float( col[4] )


            if( dtyp_update ):
                dtyp_cnt = d_ind + 1 
                dtyp_i = parmC.dtypC[dtyp_cnt]
                dtyp_i.set_g_indx(gfunc_type)
                dtyp_i.set_lmpindx(lmpindx)
                dtyp_i.setopls(k1,k2,k3,k4)
            else:
                dtyp_i = dihtype(ptype1,ptype2,ptype3,ptype4,dtype)
                dtyp_i.set_g_indx(gfunc_type)
                dtyp_i.set_lmpindx(lmpindx)
                dtyp_i.setopls(k1,k2,k3,k4)
                parmC.dtypC.put(dtyp_i)
            
            if( cnt_Dihedral_coeff >=  n_dtypes ):
                read_Dihedral_coeff = 0
                
        
        if( read_Improper_coeff and  len(col) >= 3 ):
            cnt_Improper_coeff += 1
            #AT_i = int(col[0])
            #AT_j = int(col[1])
            imp_ind = int( col[0]) - 1
            
            if( debug): print " reading imp dih type ",imp_ind," cnt ",cnt_Improper_coeff," of ",n_imptypes
            
            if( imp_ind > n_imptypes ):
                error_line =  " Error in data file index of improper parameter %d exceeds number of improper parameters %d "%(imp_ind , n_imptypes)
                error_line += " specified with dihedral types "
                print sys.exit(error_line)
                
            IMPTYPE_REF[imp_ind][0] = "??"
            IMPTYPE_REF[imp_ind][1] = "??"
            IMPTYPE_REF[imp_ind][2] = "??"
            IMPTYPE_REF[imp_ind][3] = "??"
            
            # Assume OPLS dihedral type
            IMPTYPE_F[imp_ind] = 2
            KE = float(col[1])
            Eo = float(col[2])
            IMPTYPE_E0[imp_ind] = Eo
            IMPTYPE_K[imp_ind] = KE

            ptype1 = "??"
            ptype2 = "??"
            ptype3 = "??"
            ptype4 = "??"
            # Set parameters according to type 
            g_indx = 2
            dtype = "improper"

            lmpindx = int( col[0] )


            if( imptyp_update ):
                imptyp_cnt = imp_ind + 1 
                imptyp_i = parmC.imptypC[imptyp_cnt]
                imptyp_i.set_g_indx(g_indx)
                imptyp_i.setimp(Eo,KE)
                imptyp_i.set_lmpindx(lmpindx)                
            else:
                imptyp_i = imptype(ptype1,ptype2,ptype3,ptype4,dtype)
                imptyp_i.set_g_indx(g_indx)
                imptyp_i.setimp(Eo,KE)
                imptyp_i.set_lmpindx(lmpindx)
                parmC.imptypC.put(imptyp_i)
            
            if( cnt_Improper_coeff >=  n_imptypes ):
                read_Improper_coeff = 0

            
        if( read_Atoms and len(col) >= 7 ):
            cnt_Atoms += 1
            ind = int( col[0]) - 1
            if( ind > n_atoms ):
                print sys.exit(" Error in data file index of atoms exceeds number of atoms specified with atoms ")
            
            chain_i = int(col[1]) 
            lmptype_i = int(col[2]) #- 1
            indx = int(col[2]) - 1
            q_i = float(col[3])
            m_i = ATYPE_MASS[indx]
            fftype_i = ATYPE_REF[indx]
            el = elements.getelementWithMass(m_i)
            
            if( el.symbol == "VS" ):
                el.symbol = ATYPE_l[atom_indx].strip()
                fftype_i = "VS"
                m_i = 0.0 

            # HACK !!
            if( ATYPE_MASS[indx] == 9.0  ):
                el.symbol = "LP"
                fftype_i = "LP"
                m_i = 0.0
            
            r_i =  [ float(col[4]),float(col[5]),float(col[6])] 
            type_i = str(lmptype_i)

            #tagsD = {"chain":chain_i,"symbol":el.symbol,"number":el.number,"mass":el.mass,"cov_radii":el.cov_radii,"vdw_radii":el.vdw_radii}
            #if( pt_overwrite ):

            pt_i = strucC.ptclC[ind+1]
            pt_i.position = r_i
            pt_i.charge = q_i
            pt_i.mass = m_i

            add_dict = pt_i.tagsDict
            # Set properties read in data file 
            if( set_chain_numbers ): add_dict["chain"] = chain_i
            add_dict["symbol"] = el.symbol
            add_dict["number"] = el.number
            add_dict["mass"] = el.mass
            add_dict["cov_radii"] = el.cov_radii
            add_dict["vdw_radii"] = el.vdw_radii
            add_dict["lmptype"] = lmptype_i
            if( fftype_i != "??" ):
                add_dict["fftype"] = fftype_i
            # 
            add_dict["ffmass"] = ATYPE_MASS[indx]
            pt_i.setTagsDict(add_dict)
            

        
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
                bondObj.set_lmpindx(int(col[1] ))
            else:
                bondObj = Bond( i_o, j_o )
                bondObj.set_lmpindx(int(col[1] ))
                strucC.bondC.put(bondObj)
                
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
                angleObj.set_lmpindx(int(col[1] ))
            else:
                angleObj = Angle( k_o,i_o, j_o )
                angleObj.set_lmpindx(int(col[1] ))
                strucC.angleC.put(angleObj)
            
                
            
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
                dObj.set_lmpindx(int(col[1] ))
            else:
                dObj = Dihedral( k_o,i_o, j_o,l_o )
                dObj.set_lmpindx(int(col[1] ))
                strucC.dihC.put(dObj)
                
            if( cnt_Dihedrals >=  n_dihedrals ):
                read_Dihedrals = 0
                
            
        if(read_Impropers and len(col) >= 2 ):
            cnt_Impropers += 1
            ind = int( col[0]) - 1
            
            k_o = int(col[2])
            i_o = int(col[3])
            j_o = int(col[4])
            l_o = int(col[5])
            
            if( imp_overwrite ):
                impObj = strucC.impC[cnt_Impropers]
                impObj.pgid1 = k_o
                impObj.pgid2 = i_o
                impObj.pgid3 = j_o
                impObj.pgid4 = l_o
                impObj.set_lmpindx(int(col[1] ))
                impObj.set_type("improper")
            else:
                impObj = Improper( k_o,i_o, j_o,l_o )
                impObj.set_lmpindx(int(col[1] ))
                impObj.set_type("improper")
                strucC.impC.put(impObj)
                

            if( cnt_Impropers >=  n_impropers ):
                read_Impropers = 0

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
                            
        #    cnt_Bonds += 1
        #    ind = int( col[0]) - 1
        #    BTYPE_IND[ind] = int(col[1] ) - 1
        #    BONDS[ind][0] = int(col[2])
        #    if( cnt_Bonds >=  n_atoms ):
        #        read_Bonds = 0
        #        
            

                #
    if( ljtyp_update ):
        if( ljtyp_cnt != len(parmC.ljtypC) ):
            print " Number of LJ types read in %d does not match previously read %d "%(ljtyp_cnt,len(parmC.ljtypC))

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
    return (strucC,parmC)


    """
            else:
                pt_i = Particle( r_i,type_i,q_i,m_i)
                # Set properties read in data file 
                pt_i.tagsDict["chain"] = chain_i
                pt_i.tagsDict["symbol"] = el.symbol
                pt_i.tagsDict["number"] = el.number
                pt_i.tagsDict["mass"] = el.mass
                pt_i.tagsDict["cov_radii"] = el.cov_radii
                pt_i.tagsDict["vdw_radii"] = el.vdw_radii
                pt_i.tagsDict["lmptype"] = lmptype_i
                if( fftype_i != "??" ):
                    pt_i.tagsDict["fftype"] = fftype_i
                # 
                pt_i.tagsDict["ffmass"] = ATYPE_MASS[indx]
                #pt_i.setTagsDict(tagsD)
                strucC.ptclC.put(pt_i)
    """


def add_struc(verbose,p,molecule_rep,struc_i,n_struc_i,atomic_cut,lc_expand,max_sys,max_mol_place,calc_overlap,seed):
    """
    Add structure struc_i to molecule_rep n_struc_i times via random placement
    """

    random.seed(seed)
    
    debug = False 

    cut_ij_sq = atomic_cut*atomic_cut

    max_mol_residue_number = 0
    for pid, ptclObj in struc_i.ptclC :
        if( ptclObj.tagsDict["residue"] > max_mol_residue_number ):
            max_mol_residue_number  = ptclObj.tagsDict["residue"]
    
    
    org = np.array( [0.0,0.0,0.0] ) # struc_i.get_origin()
    
    
    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    n_dim = 3
    ang_acc = 1000  # number of digets in random angle
    # Initilize
    sys_mol_n = 0     # Number of oligomers add to the system
    sys_attempts = 0    # Number of times the system has been reset
    struc_add_cnt = 0   # Total number of structures added to the final structure

    #
    # Place the atomic indices into list 
    # 
    pointIndices = range( len(struc_i.ptclC)  )
    if( debug ):
        print rank, size," splitOnProcs "
    # Create a list of atomic indices for each processor 
    myChunk  = p.splitListOnProcs(pointIndices)
    
    #
    # Start adding molecules to the system
    #
    add_mol = True
    if( rank == 0  ):
        print " Adding %d  molecules  "%n_struc_i
        if( seed == None ):
            print "   random seed will be taken from clock"
        else:
            print "   with random seed %s "%(seed)
        
    while ( add_mol ):
        #
        # Initialize 
        #
        add_mol = True
        overlap_found = True
        strucadd_atempts = 0


        # For each structure add to 
        while ( overlap_found ):
            strucadd_atempts += 1

            rot_angle_i_o = 0.0 
            rot_angle_j_o = 0.0 
            r_random_o  = np.zeros(n_dim)

            if ( rank == 0 ):
                #
                #  Get random rotation angles from single processor  
                #
                rot_angle_i_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                rot_angle_j_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                #
                #  Get random translation from single processor 
                #
                r_random_o = np.zeros(n_dim)
                for x_indx in range( n_dim ):
                    r_random_o[x_indx] = random.randrange(0,ang_acc*int(molecule_rep.latvec[x_indx][x_indx]) )/float(ang_acc)

                debug = 0
                if( debug ):
                    print " ran ",x_indx,(molecule_rep.latvec[x_indx][x_indx])
                    print rot_angle_i_o,rot_angle_j_o,r_random_o
                    #sys.exit(" Random # 's test 1")

            p.barrier() # Barrier for MPI_COMM_WORLD
            #
            # Broadcast random rotation angles and translations to all processors 
            #
            rot_angle_i = p.bcast(rot_angle_i_o)

            rot_angle_j = p.bcast(rot_angle_j_o)
            r_random = p.bcast(r_random_o)
            p.barrier() # Barrier for MPI_COMM_WORLD
            #
            # Get coordinates of randomly rotated and shifted 
            #
            shift_center_mass(struc_i,org)
            rotate(struc_i,rot_angle_i,rot_angle_j)
            vec_shift(struc_i,r_random)

            overlap = 0
            if( len(molecule_rep.ptclC) > 0 ):
                #
                # If there are particles in the system check atoms do not overlap
                #   
                if( calc_overlap ):

                    for p_i, ptclObj_i in struc_i.ptclC(myChunk):
                        r_i = np.array( ptclObj_i.position )
                        for p_sys, ptclObj_sys in molecule_rep.ptclC :
                            r_sys = np.array( ptclObj_sys.position )
                            r_ij_sq = sq_drij_c(r_i,r_sys,molecule_rep.getLatVec() )
                            if( r_ij_sq < cut_ij_sq ):
                                overlap = 1

            p.barrier() # Barrier for MPI_COMM_WORLD
            #
            # Reduce sum the overlap variable from all the processors
            #   if it is zero everywhere there was no overlap detected 
            #
            overlap_sum = p.allReduceSum(overlap)
            p.barrier() # Barrier for MPI_COMM_WORLD

            if( overlap_sum ==  0 ):

                # If no overlap detected add molecule to the system 
                sys_mol_n += 1
                struc_add_cnt += 1
                # Rest molecule numbers
                for pid, ptclObj in struc_i.ptclC :
                    ptclObj.tagsDict["chain"] = struc_add_cnt
                    if( sys_mol_n > 1 ):
                        res_numb_i = max_mol_residue_number + ptclObj.tagsDict["residue"]

                        if(debug):
                            print " residue number updat %d -> %d "%(ptclObj.tagsDict["residue"],res_numb_i)
                        ptclObj.tagsDict["residue"] = res_numb_i

                # add molecule structure to system structure
                struc_i.setLatVec(molecule_rep.getLatVec())
                molecule_rep += struc_i

                if( verbose ):
                    if( rank == 0  ):
                        print "      -  Molecule ",sys_mol_n," has been added to the system after ",strucadd_atempts," placment attempts "
                        print "         system has %d atoms and %d bonds "%(len(molecule_rep.ptclC),len(molecule_rep.bondC))
                        #print " Printing  molecule_rep bonds "
                        #molecule_rep.printbondlengths()

                overlap_found = False
                        
            else:
                overlap_found = True

            if( strucadd_atempts >= max_mol_place ):
                # If attempts to place molecule into the system exceed max set by max_mol_place

                #   reset system and star over 
                if(  rank == 0  ):
                    if(  verbose ):


                        print " add0 molecule_rep. .getLatVec() ",molecule_rep.getLatVec()


                        print "        -  Attempts to add molecule ",sys_mol_n+1," has exceeded max attempts ",max_mol_place," system will be reset for the ",sys_attempts," time "

                sys_mol_n = 0

                struc_add_cnt = 0 
                strucadd_atempts = 0
                sys_attempts += 1

                # Save lattice vectors as to no loose any expansions 
                latvec_i = molecule_rep.getLatVec()

                print " saving s1 latvec_i ",latvec_i

                # Delete system 
                del molecule_rep
                molecule_rep = StructureContainer()  # Output replicated structure
                # Set lattice vectors 
                molecule_rep.setLatVec(latvec_i) 


            if( sys_attempts >= max_sys  ):


                print " exp0 molecule_rep. .getLatVec() ",molecule_rep.getLatVec()


                # If the system has been reset over max_sys times expand the box size by lc_expand
                molecule_rep.expandLatVec(lc_expand)

                # Save lattice vectors as to no loose any expansions 
                latvec_i = molecule_rep.getLatVec()


                print " exp1 molecule_rep. .getLatVec() ",molecule_rep.getLatVec()

                # Delete system 
                del molecule_rep
                molecule_rep = StructureContainer()  # Output replicated structure
                # Set lattice vectors 
                molecule_rep.setLatVec(latvec_i) 

                sys_attempts = 0

                if( verbose ):                
                    if( rank == 0  ):

                        print '          - Number of system resets has exceeded the maximum  (option max_sys) ',max_sys
                        print '          - Lattice vectors will be expanded by (option lc_expand)',lc_expand
                        print '             v_1 ',latvec_i[0]
                        print '             v_2 ',latvec_i[1]
                        print '             v_3 ',latvec_i[2]

        p.barrier() # Barrier for MPI_COMM_WORLD


        if( sys_mol_n ==  n_struc_i  ):
            # If all the molecule have been added exit while loop and print system 
            add_mol = False
            latvec_mol = molecule_rep.getLatVec()
            p.barrier() # Barrier for MPI_COMM_WORLD
            if( verbose and rank == 0  ):
                print " All molecules  have been added "

    return molecule_rep
    


def add_struc_ongrid(verbose,p,molecule_rep,struc_i,n_sol_l,atomic_cut,lc_expand,max_sys,max_mol_place,calc_overlap):
    """
    Add structure struc_i to molecule_rep n_struc_i times onto a grid 
    """

    molecule_rep_mass = molecule_rep.getTotMass()
    sol_mass = struc_i.getTotMass()

    cut_ij_sq = atomic_cut*atomic_cut
    
    sol_buf = atomic_cut
    
    max_residue_number_i = 0
    max_chain_number_i = 0
    for pid, ptclObj in molecule_rep.ptclC :
        if( ptclObj.tagsDict["residue"] > max_residue_number_i ): max_residue_number_i  = ptclObj.tagsDict["residue"]
        if( ptclObj.tagsDict["chain"] > max_chain_number_i ): max_chain_number_i  = ptclObj.tagsDict["chain"]

    min_sol_residue = 100000000
    max_sol_residue_number = -100000000
    min_sol_chain = 100000000
    max_sol_chain_number = -100000000
    for pid, ptclObj in struc_i.ptclC :
        if( ptclObj.tagsDict["residue"] < min_sol_residue ): min_sol_residue =ptclObj.tagsDict["residue"]
        if( ptclObj.tagsDict["residue"] > max_sol_residue_number ): max_sol_residue_number =ptclObj.tagsDict["residue"]
        if( ptclObj.tagsDict["chain"] < min_sol_chain ): min_sol_chain =ptclObj.tagsDict["chain"]
        if( ptclObj.tagsDict["chain"] > max_sol_chain_number ): max_sol_chain_number =ptclObj.tagsDict["chain"]

    delt_reside_number = max_sol_residue_number - min_sol_residue +1
    delt_chain_number = max_sol_chain_number - min_sol_chain +1


    print "max_residue_number_i ",max_residue_number_i
    print "min_sol_residue ",min_sol_residue
    print "max_sol_residue_number ",max_sol_residue_number
    print "delt_reside_number ",delt_reside_number
    print ""    
    print "max_chain_number_i ",max_chain_number_i
    print "delt_chain_number ",delt_chain_number
    
    for pid, ptclObj in struc_i.ptclC :
        ptclObj.tagsDict["residue"] = max_residue_number_i 
        ptclObj.tagsDict["chain"] = max_chain_number_i
        
        # print pid,"residue", ptclObj.tagsDict["residue"],ptclObj.tagsDict["chain"]
    
    sol_maxlength = getlength(struc_i)
    org =  np.array( [0.0,0.0,0.0] ) #get_origin(molecule_rep)

    debug = False 
    
    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    n_dim = 3
    ang_acc = 1000  # number of digets in random angle
    # Initilize
    sys_mol_n = 0     # Number of oligomers add to the system
    sys_attempts = 0    # Number of times the system has been reset
    struc_add_cnt = 0   # Total number of structures added to the final structure

    latvec = molecule_rep.latvec
    len_target_ang = latvec[0][0]
    # Check to be sure the grid is large enough
    sol_box_side = int(math.ceil(n_sol_l**(1.0/3.0) ) )    # Number of solvents per box side 
    sol_length = sol_maxlength + sol_buf                   # length of solvent 
    sol_length_sq = sol_length*sol_length                  # length of solvent squared for overlap calculati

    len_sol_box_ang = sol_length*sol_box_side
    if( len_sol_box_ang > len_target_ang): len_target_ang = len_sol_box_ang
    sol_grid_len = len_target_ang/float(sol_box_side)

    sol_total_n = sol_box_side*sol_box_side*sol_box_side
    sol_extra_n = sol_total_n - n_sol_l 
    
    if( rank == 0  ):

        print "n_sol_l**(1.0/3.0)  ",n_sol_l**(1.0/3.0) 
        print "n_sol_l ",n_sol_l 
        print "sol_maxlength ", sol_maxlength
        print "sol_buf ", sol_buf
        print "sol_box_side ", sol_box_side
        print ""
        print "sol_length ", sol_length
        print "len_sol_box_ang ", len_sol_box_ang
        print "sol_grid_len ", sol_grid_len
        print "sol_total_n  ",sol_total_n
        print "sol_extra_n ",sol_extra_n
        print "len_target_ang ",len_target_ang 
        print ""
        print "Length of original target box size ",len_target_ang
        print "Number of solvent molecules per lengths ",sol_box_side
        print "Gid spaceing    ",sol_grid_len
        

    # Set lattice vector to new box size
    latvec_list = [np.array([len_target_ang,0.0,0.0]),np.array( [0.0,len_target_ang,0.0]),np.array( [0.0,0.0,len_target_ang]) ]

    # Solvent molecule container
    sol_rep = StructureContainer()  
    sol_rep.setLatVec(latvec_list)
    molecule_rep.setLatVec(latvec_list)

    # If options set for fixed then set random seed (for regression testing)
    #if fixed_rnd_seed:
    #    random.seed(0)

    sys_sol_n = 0     # Number of solvent add to the system
    sys_attempts = 0  # Number of times the system has been reset
    #
    # Start adding molecules to the system
    #
    if( rank == 0  ):
        print " Adding %d  solvent  "%n_sol_l
        
    add_sol = True
    while ( add_sol ):
        #
        # Initialize 
        #
        strucadd_atempts = 0

        bad_grid_point = False 
        #
        # 
        #
        if( rank == 0  ):
            print " starting %d of %d "%(sys_sol_n,n_sol_l)

        # loop over the number of times each solvent in the solvent list needs to be replicated
        if( bad_grid_point ): break
        if( not add_sol): break
        # Check overlap with other particles in original structure container
                
        overlap_found = True

        if( bad_grid_point ): break
        if( not add_sol): break 
        while ( overlap_found ):

            strucadd_atempts += 1
            for x_indx in range(sol_box_side):
                if( not add_sol): break 
                for y_indx in range(sol_box_side):
                    if( not add_sol): break 
                    for z_indx in range(sol_box_side):

                        if( sys_sol_n == n_sol_l ):
                            add_sol = False
                            break

                        l_x =  float(x_indx)*sol_grid_len
                        l_y =  float(y_indx)*sol_grid_len
                        l_z =  float(z_indx)*sol_grid_len

                        if( debug and  rank == 0 and  calc_overlap != 0 ):
                            print " Checking overlap for solvent %d at lattice point %f %f %f "%(sys_sol_n,l_x,l_y,l_z)

                        if( l_x > sol_rep.latvec[0][0] or l_y > sol_rep.latvec[1][1] or l_z > sol_rep.latvec[2][2] ):
                            if( rank == 0 ):
                                print " Lattic point beyond box %f %f %f "%(sol_rep.latvec[0][0],sol_rep.latvec[1][1], sol_rep.latvec[2][2])
                            bad_grid_point = True
                            break 

                        lat_pos = np.array( [l_x,l_y,l_z] )

                        # Make sure there is no overlap with the added molecules
                        overlap = 0
                        if( calc_overlap == 1  ):

                            for p_i, ptclObj_i in molecule_rep.ptclC :
                                r_i = np.array( ptclObj_i.position )
                                r_ij_sq = sq_drij_c(r_i,lat_pos,molecule_rep.getLatVec() )
                                if( r_ij_sq < sol_length_sq ):
                                    overlap = 1

                        p.barrier() # Barrier for MPI_COMM_WORLD
                        #
                        # Reduce sum the overlap variable from all the processors
                        #   if it is zero everywhere there was no overlap detected 
                        #
                        overlap_sum = p.allReduceSum(overlap)

                        if( overlap_sum ==  0 ):
                            # If no overlap detected add molecule to the system 
                            sys_sol_n += 1
                            overlap_found = False 
                            # Shift molecule to lattice point
                            shift_center_mass(struc_i,org)
                            vec_shift(struc_i,lat_pos)

                            struc_i.setLatVec(sol_rep.getLatVec())

                            if( debug):
                                print sys_sol_n,lat_pos,struc_i.center_mass()

                            debug_resid = False 
                            # Rest molecule numbers
                            for pid, ptclObj in struc_i.ptclC :
                                #ptclObj.tagsDict["chain"] = sys_sol_n
                                #if( sys_sol_n > 1 ):
                                res_numb_i = delt_reside_number + ptclObj.tagsDict["residue"]
                                chain_numb_i = delt_chain_number + ptclObj.tagsDict["chain"]
                                if(debug_resid):
                                        print " solvent residue number updat %d -> %d residue %d "%(ptclObj.tagsDict["residue"],res_numb_i,sys_sol_n)
                                        print " solvent residue number updat %d -> %d chain %d "%(ptclObj.tagsDict["chain"],chain_numb_i,sys_sol_n)

                                ptclObj.tagsDict["residue"] = res_numb_i
                                ptclObj.tagsDict["chain"] = chain_numb_i
                                #else:
                                #    if(debug):
                                #        print " solvent residue number updat %d -> %d chain %d "%(ptclObj.tagsDict["residue"],initial_sol_residue_number+1,sys_sol_n)
                                #
                                #   ptclObj.tagsDict["residue"] = max_residue_number_i + ptclObj.tagsDict["residue"] - min_sol_residue
                                    #                                    #
                            sol_rep += struc_i

                            if( verbose ):
                                if( rank == 0  ):
                                    print "      -  Molecule %d  has been added to the system at lattice point %f %f %f  "%(sys_sol_n,l_x,l_y,l_z)

                        else:
                            overlap_found = True
                            if( debug): print "  lattice point was found to overlap "



        if( sys_sol_n == n_sol_l ):
            add_sol = False
            p.barrier() # Barrier for MPI_COMM_WORLD
            if( verbose and rank == 0  ):
                print " All solvents  have been added "

        else:
            #
            # If attempts to place solvent molecule into the system failed
            #
            sys_attempts += 1
            if( rank == 0  ):
                print '        - Placment failed  after adding %d particles and %d attempts   '%(sys_sol_n,sys_attempts)
            
            sys_sol_n = 0                
            #
            # Save lattice vectors as to not loose any expansions
            #
            latvec_i = sol_rep.getLatVec()
            #
            # Delete system
            #
            del sol_rep
            sol_rep = StructureContainer()  # Output replicated structure
            #
            # Set lattice vectors
            #
            sol_rep.setLatVec(latvec_i)
            #
            #
            shrink_grid  = False
            if( shrink_grid ):

                if(  (sol_grid_len*0.99) > sol_length):
                    # If grid spacing is larger than the solvent length shrink grid
                    sol_grid_len = sol_grid_len*0.99

                    #if( verbose ):                
                    if( rank == 0  ):
                        print '          - Reducing grid space to  %f which is still larger than the molecular length %f  '%(sol_grid_len,sol_length)

                else:
                    # Otherwise increase volume and set grid spacing to solvent length
                    sol_grid_len = sol_length
                    #
                    # Expand the box size by a single solvent length 
                    #
                    #sol_rep.expandLatVec(options.lc_expand)
                    sol_rep.latvec[0][0] = sol_rep.latvec[0][0] + sol_length
                    sol_rep.latvec[1][1] = sol_rep.latvec[1][1] + sol_length
                    sol_rep.latvec[2][2] = sol_rep.latvec[2][2] + sol_length
                    #
                    # Increase the number of solvent molecules along the box by 1
                    #
                    sol_box_side += 1
                    

                    #if( verbose ):                
                    if( rank == 0  ):
                        print '          - Lattice vectors will be expanded by solvent length %f ',sol_length
            else:

                #
                # Expand the box size by a Grid spacing 
                #
                expand_len  = sol_grid_len
                #sol_rep.expandLatVec(options.lc_expand)
                sol_rep.latvec[0][0] = sol_rep.latvec[0][0] + expand_len
                sol_rep.latvec[1][1] = sol_rep.latvec[1][1] + expand_len
                sol_rep.latvec[2][2] = sol_rep.latvec[2][2] + expand_len
                # sol_grid_len += sol_grid_len + float(int(sol_length/sol_box_side))
                #
                # Increase the number of solvent molecules along the box by 1
                #
                sol_box_side += 1

                max_sol_residue_number = 0
                max_sol_chain_number = 0 
                for pid, ptclObj in struc_i.ptclC :
                    if( ptclObj.tagsDict["residue"] > max_sol_residue_number ): max_sol_residue_number =ptclObj.tagsDict["residue"]
                    if( ptclObj.tagsDict["chain"] > max_sol_chain_number ): max_sol_chain_number =ptclObj.tagsDict["chain"]

                for pid, ptclObj in struc_i.ptclC :
                    ptclObj.tagsDict["residue"] = ptclObj.tagsDict["residue"] - max_sol_residue_number + max_residue_number_i
                    ptclObj.tagsDict["chain"] = ptclObj.tagsDict["chain"]  - max_sol_chain_number  + max_chain_number_i

                
                #if( verbose ):                
                if( rank == 0  ):
                    print '          - Lattice vectors will be expanded by solvent length %f ',sol_length
                        
        p.barrier() # Barrier for MPI_COMM_WORLD


    # Add replicated solvents to final structure
    molecule_rep += sol_rep
    molecule_rep.setLatVec(sol_rep.latvec)

    return molecule_rep

def write_data(strucC,parmC,data_file):

    """
    Write data file
    """
    import sys, math
    # 
    #
    # ' print lammps data file '
    #

    debug = False 
    if( debug ):
        print "int(len(strucC.dihC)) ",int(len(strucC.dihC))
        print "imp_cnt",imp_cnt
        print "int(len(strucC.dihC))",int(len(strucC.dihC))
        print "int(len( parmC.dtypC))",int(len( parmC.dtypC))
        print "int( len(parmC.imptypC))",int( len(parmC.imptypC))
        sys.exit(" dih cnt debug ")
                
    # Calculate totals
    n_atoms = len( strucC.ptclC  )
    n_bonds = int(len(strucC.bondC))
    n_angles = int(len(strucC.angleC))
    n_dihedrals = int(len(strucC.dihC))
    n_impropers =  int(len(strucC.impC))
    
    n_atypes = int(len( parmC.ljtypC )) #+ 1
    n_btypes = int(len( parmC.btypC )) #+ 1
    n_angtypes = int(len( parmC.atypC )) #+ 1
    n_dtypes = int(len( parmC.dtypC)) 
    imp_cnt = int( len(parmC.imptypC))
    if(imp_cnt > 0 ):
        n_imptypes = imp_cnt
    else:
        n_imptypes = 1
        
    # Calculate box size
    latvec = strucC.getLatVec()
    
    bmn_x = latvec[0][0]/-2.0
    bmx_x = latvec[0][0]/2.0
    bmn_y = latvec[1][1]/-2.0
    bmx_y = latvec[1][1]/2.0
    bmn_z = latvec[2][2]/-2.0
    bmx_z = latvec[2][2]/2.0

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
    for lj_p, ljObj_p  in parmC.ljtypC:
        F.write( "%10d %16.8f   # %5s \n" % ( lj_p , ljObj_p.get_mass() , ljObj_p.get_ptype1()  ) )
    F.write('\n')
    F.write(' Pair Coeffs \n')
    F.write('\n')
    for lj_p, ljObj_p  in parmC.ljtypC:
        F.write( "%10d %12.6f %12.6f  \n" % (lj_p, ljObj_p.get_epsilon(),  ljObj_p.get_sigma()  ) )
    F.write('\n')
    # Write Bond Coeffs
    if( len(parmC.btypC) > 0 ):
        F.write(' Bond Coeffs \n')
        F.write('\n')
        for btyp_p, btypObj_p  in parmC.btypC:    
            if( btypObj_p.get_type() == "harmonic"):
                F.write( "%10d %12.6f %12.6f # %5s %5s  \n" % (btyp_p,btypObj_p.get_kb(),btypObj_p.get_r0(), btypObj_p.get_ptype1(), btypObj_p.get_ptype2() ) )
        F.write('\n')

    # Write Angle Coeffs
    if( len(parmC.atypC) > 0 ):
        F.write(' Angle Coeffs \n')
        F.write('\n')
        for atyp_p, atypObj_p  in parmC.atypC:    
            if( atypObj_p.get_type() == "harmonic"):
                F.write( "%10d %12.6f %12.6f # %5s %5s  %5s   \n" % (atyp_p,atypObj_p.get_kb(),atypObj_p.get_theta0(), atypObj_p.get_ptype1(),atypObj_p.get_ptype2(),atypObj_p.get_ptype3() ) )
        F.write('\n')

    # Write Dihedral Coeffs
    if( len(parmC.dtypC) > 0 ):
        F.write(' Dihedral Coeffs \n')
        F.write('\n')
        imp_cnt = 0 
        for dtyp_p, dtypObj_p  in parmC.dtypC:    
            if( dtypObj_p.get_type() == "multiharmonic"):
                # K = K[1+d cons(n theta)]
                d = dtypObj_p.get_theat_s()
                K = dtypObj_p.get_kb()
                n = dtypObj_p.get_mult()
                w = 0.0 # Weight 
                F.write( "%10d %12.6f %d  %d %12.6f # %5s %5s  %5s  %5s   \n" % (dtyp_p-imp_cnt,K,n,d,w, dtypObj_p.get_ptype1(),dtypObj_p.get_ptype2(),dtypObj_p.get_ptype3(),dtypObj_p.get_ptype4()  ) )
            elif( dtypObj_p.get_type() == "rb" or  dtypObj_p.get_type() == "opls"  ):
                # Get opls parameters
                klist = dtypObj_p.get_oplsklist()
                F.write( "%10d  %12.6f  %12.6f  %12.6f  %12.6f # %5s %5s  %5s %5s \n" % (dtyp_p,klist[0],klist[1],klist[2],klist[3], dtypObj_p.get_ptype1(),dtypObj_p.get_ptype2(),dtypObj_p.get_ptype3(),dtypObj_p.get_ptype4()  ) )
            elif( dtypObj_p.get_type() == "improper"):
                error_line = " Will print improper later "
                imp_cnt += 1 
            else:
                error_line = " Unknow dihedral type %s "%(dtypObj_p.get_type() )
                sys.exit(error_line)


        F.write('\n')
    
    
    F.write(' Improper Coeffs \n')
    F.write('\n')
    if( len(parmC.imptypC) > 0 ):
        for imptyp_p, imptypObj_p  in parmC.imptypC:    
            if( imptypObj_p.get_type() == "improper"):
                e0,ke = imptypObj_p.getimp()
                F.write( "%10d %12.6f %12.6f # %5s %5s  %5s  %5s   \n" % (imptyp_p,ke,e0, imptypObj_p.get_ptype1(),imptypObj_p.get_ptype2(),imptypObj_p.get_ptype3(),imptypObj_p.get_ptype4()  ) )
    else:
        F.write( "    1 0.0 0.0  \n")

        
    F.write('\n')
    # Write Particles
    if( len(strucC.ptclC) > 0 ):

        F.write(' Atoms \n')
        F.write('\n')
        TOTAL_CHARGE = 0.0
        for pid_o, ptclObj_o  in strucC.ptclC:
            fftype_i = ptclObj_o.tagsDict["fftype"]
            chain_i = ptclObj_o.tagsDict["chain"]
            charge_i = ptclObj_o.charge
            #type_indx = int(ptclObj_o.type)
            lmptype_i = ptclObj_o.tagsDict["lmptype"]
            r_i = ptclObj_o.position            
            F.write( "%9d %9d %8d %12.8f %12.6f %12.6f %12.6f # %5s \n" % (pid_o,chain_i,lmptype_i,charge_i,r_i[0],r_i[1],r_i[2] ,fftype_i)  )
            TOTAL_CHARGE = TOTAL_CHARGE + float( charge_i )

        F.write('\n')
    # Write Bonds
    if( len(strucC.bondC) > 0 ):
        F.write(' Bonds \n')
        F.write('\n')
        b_cnt = 0 
        for b_o, bondObj_o  in strucC.bondC:
            #
            b_cnt += 1
            AT_i =  strucC.ptclC[ bondObj_o.pgid1 ].tagsDict["fftype"]
            AT_j =  strucC.ptclC[ bondObj_o.pgid2 ].tagsDict["fftype"]
            b_ind = int(bondObj_o.get_lmpindx())
            #
            F.write(  '%9d %8d %9d %9d # %5s %5s \n' % (b_cnt,b_ind,bondObj_o.pgid1,bondObj_o.pgid2, AT_i, AT_j ) )
        F.write('\n')
    
    # Write Angles
    if( len(strucC.angleC) > 0 ):

        F.write(' Angles \n')
        F.write('\n')
        a_cnt = 0 
        for a_o, angleObj_o in strucC.angleC:
            a_cnt += 1
            a_k = angleObj_o.pgid1
            a_i = angleObj_o.pgid2
            a_j = angleObj_o.pgid3
            a_ind = int(angleObj_o.get_lmpindx())
            F.write(  '%9d %8d %9d %9d %9d \n' % (a_cnt,a_ind,a_k,a_i,a_j) )
        F.write(  '\n' )
    
    # Write Dihedrals
    if( len(strucC.dihC) > 0 ):

        F.write(' Dihedrals \n')
        F.write('\n')
        dih_cnt = 0 
        for d_o,dihObj_o in strucC.dihC:
            dih_cnt += 1

            d_k = dihObj_o.pgid1
            d_i = dihObj_o.pgid2
            d_j = dihObj_o.pgid3
            d_l = dihObj_o.pgid4
            d_ind = int(dihObj_o.get_lmpindx())

            AT_k = strucC.ptclC[ d_k ].tagsDict["fftype"]
            AT_i = strucC.ptclC[ d_i ].tagsDict["fftype"]
            AT_i = strucC.ptclC[ d_j ].tagsDict["fftype"]
            AT_l = strucC.ptclC[ d_l ].tagsDict["fftype"]

            #error_line = " dih index is not in bounds "
            #error_line += " for atoms %d  %d  %d  %d  "%(d_k,d_i,d_j,d_l)
            #error_line += " for type atoms %s %s %s %s "%(AT_k,AT_i,AT_j,AT_l)
            #sys.exit(error_line)
            F.write(  '%9d %8d %9d %9d %9d %9d # %s %s %s %s \n' % (dih_cnt,d_ind,d_k,d_i,d_j,d_l,AT_k,AT_i,AT_j,AT_l) )
            
        F.write( '\n' )

    # Write Impropers
    if( len(strucC.impC) > 0 ):
        F.write(' Impropers \n')
        F.write('\n')
        imp_cnt = 0 
        for d_o,impObj_o in strucC.impC:
            if( impObj_o.get_type() == "improper"   ):
                imp_cnt += 1 
                d_k = impObj_o.pgid1
                d_i = impObj_o.pgid2
                d_j = impObj_o.pgid3
                d_l = impObj_o.pgid4
                imp_ind = int(impObj_o.get_lmpindx())
                F.write(  '%9d %8d %9d %9d %9d %9d \n' % (imp_cnt,imp_ind,d_k,d_i,d_j,d_l) )
            else:
                error_line = " non imporper type in improper container %s "%( impObj_o.get_type())
                sys.exit(error_line)
            
        F.write( '\n' )            

    F.close()


def update_data_prop(strucC):
    """
     Add in particle properties not included 
     """

    # Load periodic table 
    # elements = periodictable()
    
    for pid, pt_i  in strucC.ptclC:
        add_dict = pt_i.tagsDict
        add_dict["residue"] = 1
        add_dict["resname"] = "MOLRES"
        add_dict["gtype"] = pt_i.tagsDict["symbol"]
        pt_i.setTagsDict(add_dict)

        
    return (strucC)


def write_xmol(strucC,data_file):

    """
    Write a structure  to an xmol file

    Args:
        data_file    (str) xmol file name
    Reutrns
        null
    """
    # Open xmol file 
    F = open(data_file,"w")
    # Loop over structures
    comment = " generated by com2xyz.py "
    F.write(" %d \n" % len( strucC.ptclC ) )
    F.write(" %s \n"%comment)
    for pid, ptclObj  in strucC.ptclC:
        r_i = ptclObj.position
        symbol = str(ptclObj.tagsDict['symbol'])
        F.write( " %5s %16.8f %16.8f %16.8f \n"  % (symbol ,float(r_i[0]), float(r_i[1]),float(r_i[2]) ) )   
    F.close()

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
    
def main():
    """
    Read in data file and replicate it 
    """

    #        
    # Read options 
    #
    options, args = get_options()


    # Initialize mpi
    p = mpiBase.getMPIObject()

    all_mol_struc = StructureContainer()
    replicate_struc = StructureContainer()
    replicate_prama = ParameterContainer()



    if( len(options.mol1_data) > 0 and options.mol1_n > 0 ):
        mol1_struc = StructureContainer()
        mol1_struc , replicate_prama = read_data( mol1_struc , replicate_prama , options.mol1_data )
        update_data_prop(mol1_struc)
        print str(mol1_struc)
        all_mol_struc += mol1_struc

    if( len(options.mol2_data) > 0 and options.mol2_n > 0 ):
        mol2_struc = StructureContainer()
        mol2_struc , replicate_prama = read_data( mol2_struc , replicate_prama , options.mol2_data )
        update_data_prop(mol2_struc)
        print str(mol1_struc)
        all_mol_struc += mol2_struc


    if( len(options.mol3_data) > 0 and options.mol3_n > 0 ):
        mol3_struc = StructureContainer()
        mol3_struc , replicate_prama = read_data( mol3_struc , replicate_prama , options.mol3_data )
        update_data_prop(mol3_struc)
        print str(mol1_struc)
        all_mol_struc += mol3_struc

    
    # Don't need since not using norm_dihparam
    norm_dihparam = False 
    cov_nblist = []
    cov_nbindx = []
    #replicate_prama_clean,all_mol_struc  = set_param(all_mol_struc,replicate_prama,norm_dihparam,cov_nblist, cov_nbindx)


    if( len(options.mol1_data) > 0 and options.mol1_n > 0 ):
        replicate_struc = add_struc(options.verbose,p,replicate_struc,mol1_struc,options.mol1_n,options.atomic_cut,options.lc_expand,options.max_sys,options.max_mol_place,options.calc_overlap,options.seed)
    if( len(options.mol2_data) > 0 and options.mol2_n > 0 ):
        replicate_struc = add_struc(options.verbose,p,replicate_struc,mol2_struc,options.mol2_n,options.atomic_cut,options.lc_expand,options.max_sys,options.max_mol_place,options.calc_overlap,options.seed)
    if( len(options.mol3_data) > 0 and options.mol3_n > 0 ):
        replicate_struc = add_struc_ongrid(options.verbose,p,replicate_struc,mol3_struc,options.mol3_n,options.atomic_cut,options.lc_expand,options.max_sys,options.max_mol_place,options.calc_overlap)

    write_data(replicate_struc,replicate_prama,options.out_data)
    write_xmol(replicate_struc,options.out_xyz)

if __name__=="__main__":
    main()
