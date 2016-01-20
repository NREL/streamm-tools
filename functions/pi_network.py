#! /usr/bin/env python
"""
Radial distribution  code

 length - Angstroms
 mass   - AMU
 volume - Angstroms^3
"""

# Dr. Travis Kemper
# NREL
# Initial Date 6/30/2014
# travis.kemper@nrel.gov

import numpy as np
import datetime, os 
import copy

import json, math , sys

from itertools import izip

#MDanalysis
try:
    from MDAnalysis import *
    from MDAnalysis.core.distances import * ##distance_array
    #import MDAnalysis.core.units            # for bulk water density
except:
    import sys
    print "MDAnalysis module not build/configured correctly"
    sys.exit(0)

# Streamm toolkit modules 
from structureContainer import StructureContainer
from buildingblocks import Buildingblock
from particles     import Particle, ParticleContainer
from bonds         import Bond,     BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer
from impropers     import Improper, ImproperContainer

from parameters    import ParameterContainer
from parameters    import ljtype,   LJtypesContainer
from parameters    import bondtype, BondtypesContainer
from parameters    import angletype,AngletypesContainer
from parameters    import dihtype,  DihtypesContainer
from parameters    import imptype,  ImptypesContainer

from periodictable import periodictable

import units
import mpiBase


def get_options():
    """
    Set options
    """
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-o","--output_id", dest="output_id", default="pipi",type="string",help=" prefix for output files  ")
    #
    # Input files
    #
    parser.add_option("--in_cply", dest="in_cply", type="string", default="", help="Input cply file")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_dcd", dest="in_dcd", type="string", default="", help="Input trajectory file in compressed dcd format ")
    parser.add_option("--in_xtc", dest="in_xtc", type="string", default="", help="Input xtc file with atoms listed as atom type numbers")  
    #
    # Groups
    #
    parser.add_option("--group_chains", dest="group_chains", default=True,action="store_true", help="Group by molecules ")
    parser.add_option("--group_residues", dest="group_residues", default=False,action="store_true", help="Group by molecules ")    
    #
    # Frames
    #
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=False,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")
    #
    #
    # Filters
    #
    parser.add_option("--symbol", dest="symbol", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--label", dest="label", type="string", default="", help="select atoms of group by label  ")    
    parser.add_option("--fftype", dest="fftype", type="string", default="", help="select atoms of group by force field type  ")    
    parser.add_option("--chain", dest="chain", type="string", default="", help="select atoms of group by chain/molecule number  ")    
    parser.add_option("--resname", dest="resname", type="string", default="", help="select atoms of group by residue name  ")    
    parser.add_option("--residue", dest="residue", type="string", default="", help="select atoms of group by resudue number  ")    
    parser.add_option("--ring", dest="ring", type="string", default="", help="select atoms of group by particlesn a ring   ")    

    parser.add_option("--r_cut", dest="r_cut",  default=1.2,type=float, help=" e-e distance cut off ")
    parser.add_option("--r_e", dest="r_e",  default=0.6,type=float, help="distance of center of electron sphere above a conjugated carbon  ")
    parser.add_option("--grp_time", dest="grp_time",  default=False,action="store_true", help="print time series data for each group")
    
    (options, args) = parser.parse_args()
        
    return options, args
   
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

def addtagDic(dic_i,tag,tag_str,setint=False):
    """
    Take a string from input split it into values and add it to a dictionary list
    """
    if( len( tag_str ) ):
        dic_i[tag] = []
        for id_s in tag_str.split():
            if( setint ):
                dic_i[tag].append(int(id_s))
            else:
                dic_i[tag].append(id_s)
            
    return dic_i

def create_search(search_dic,f_symb,f_label,f_fftype,f_residue,f_resname,f_chain,f_ring):
    """
    Create a dictionary to pass to particle search
    """

    search_dic = addtagDic(search_dic,"symbol",f_symb)
    search_dic = addtagDic(search_dic,"label",f_label)
    search_dic = addtagDic(search_dic,"fftype",f_fftype)
    search_dic = addtagDic(search_dic,"residue",f_residue,setint=True)
    search_dic = addtagDic(search_dic,"resname",f_resname)
    search_dic = addtagDic(search_dic,"chain",f_chain,setint=True)
    search_dic = addtagDic(search_dic,"ring",f_ring,setint=True)
    
    return search_dic


def read_lmpdata( strucC , parmC , data_file, coordupdate=False):
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
        
    debug = False
    verbose = True

    set_chain_numbers = True

    if( not set_chain_numbers ):
        print " Warning not reading in chain numbers!!! "

    # Load periodic table 
    pt = periodictable()
    

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
            el = pt.getelementWithMass(m_i)
            
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

            pt_i = strucC.ptclC[ind+1]
            pt_i.position = r_i

            add_dict = pt_i.tagsDict
            # Not in cply file
            add_dict["lmptype"] = lmptype_i
            if( not coordupdate ):
                # Do not over write data from cply file
                pt_i.charge = q_i
                pt_i.mass = m_i
                add_dict["chain"] = chain_i
                add_dict["symbol"] = el.symbol
                add_dict["number"] = el.number
                add_dict["mass"] = el.mass
                add_dict["cov_radii"] = el.cov_radii
                add_dict["vdw_radii"] = el.vdw_radii

                add_dict["fftype"] = fftype_i
            #
                add_dict["ffmass"] = ATYPE_MASS[indx]
                add_dict["qgroup"] = chain_i
                add_dict["residue"] = chain_i
                add_dict["resname"] = "MOLR"
                add_dict["label"] = add_dict["symbol"]
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

def read_gro(strucC,in_gro,coordupdate=False,set_chaintoresidue=False,verbose=False,debug = False):
    """
    Read gromacs structure file

    Arguments:
        struc_o (StructureContainer)
        in_gro  (str) GROMACS .gro file

    Returns:
        struc_o (StructureContainer)
    """
    # atomicpy functions

    

    try:
        with open(in_gro,'r') as F:
            Lines = F.readlines()
            F.close()
    except IOError:
        print " Specified .gro file ",in_gro," does not exisit "
        sys.exit("Invalid file ")

    # Check to see if a previous read has occured
    pt_update = False 
    n_pt = int( Lines[1])
    if( len(strucC.ptclC) > 0 ):
        pt_update = True
        # Check of conistent number of atoms
        pt_cnt = len(strucC.ptclC)
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
            residue_i = int(line[0:5].strip())
            resname_i = line[5:10].strip()            
            g = line[10:15].strip()
            particle_i = int(line[15:20].strip())
            x = units.convert_nm_angstroms( float( line[20:28] ))
            y = units.convert_nm_angstroms( float(line[28:36]))
            z = units.convert_nm_angstroms( float(line[36:44]))
            #r_i =  numpy.array(  [float(x)*10,float(y)*10,float(z)*10] )
            r_i =    [x,y,z]
            if(debug):
                print " particle ",ptcl_cnt,g,r_i
            if( pt_update ):
                pt_i = strucC.ptclC[ptcl_cnt]
            else:
                pt_i = Particle(  )
            pt_i.position = r_i
            if( not coordupdate ):
                add_dict = pt_i.tagsDict
                add_dict["residue"] = int(residue_i)
                add_dict["resname"] = resname_i
                add_dict["label"] = str(g)
                # set as defualts 
                add_dict["fftype"] = "??"
                add_dict["qgroup"] = 1
                if( set_chaintoresidue ):
                    add_dict["chain"] = int(residue_i)
                else:
                    add_dict["chain"] = 1 
                pt_i.setTagsDict(add_dict)

            if( not pt_update ):
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
        gro_lines += "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"  % (atom_indx,pt_i.tagsDict["resname"][:5],pt_i.tagsDict["label"][:5],atom_indx,r_i_nm[0],r_i_nm[1],r_i_nm[2] )
        if( atom_indx > 99999 ):
            atom_indx = 1
            
    gro_lines += " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (units.convert_angstroms_nm(latvec[0][0]),units.convert_angstroms_nm(latvec[1][1]),units.convert_angstroms_nm(latvec[2][2]),units.convert_angstroms_nm(latvec[0][1]),units.convert_angstroms_nm(latvec[0][2]),units.convert_angstroms_nm(latvec[1][0]),units.convert_angstroms_nm(latvec[1][2]),units.convert_angstroms_nm(latvec[2][0]),units.convert_angstroms_nm(latvec[2][1])) 

    F = open( data_file, 'w' )
    F.write(gro_lines)
    F.close()
                    
def sigma_m(N,ave,ave_sq):
    """
    Calculate the standard deviation of the mean for a confidence
    interval of 90%. Will return zero for single valued data sets 
    """
    import numpy
    
    # Some website that probably does not exist
    #   http://mathworld.wolfram.com/Studentst-Distribution.html
    #   http://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
    #
    
    v = N - 1 #  Degrees of freedom
    
    # Set appropriate Students t prefactor 
    #if( v > 100 ):
    #	Coefint_pre = 1.64487
    #el
    
    if( v > 30 ):
	Coefint_pre = 1.66023
    elif( v > 10 ):
	Coefint_pre = 1.69726
    elif( v > 5 ):
	Coefint_pre = 1.81246
    elif( v == 5 ):
	Coefint_pre = 2.01505
    elif( v == 4 ):
	Coefint_pre = 2.13185
    elif( v == 3 ):
	Coefint_pre = 2.35336
    elif( v == 2 ):
	Coefint_pre = 2.91999
    elif( v == 1 ):
	Coefint_pre = 6.31375  
	
    if( N > 1 ):
	v_sqrt = numpy.sqrt(  N - 1 )
	sigma = numpy.sqrt(  ( ave_sq ) - (ave)**2 ) # Standard deviation
	sigma_mean = Coefint_pre*sigma/v_sqrt
    else:
	sigma = 0.0  # Set zero for unknow error 
	sigma_mean = 0.0  # Set zero for unknow error 
    
    return sigma,sigma_mean

def asphericity(Rnm_eg):
    """
    Calculate the asphericity from the eiganvalues of the radius of gyration tensor 
    """
    num = (Rnm_eg[0] - Rnm_eg[2])**2 + (Rnm_eg[1] - Rnm_eg[2])**2 + (Rnm_eg[0] - Rnm_eg[1])**2
    dem = 2*(Rnm_eg[0] + Rnm_eg[1] + Rnm_eg[2])**2
    Asphere = num/dem
    return Asphere



def molpbcs(strucC, cov_nblist, cov_nbindx ,verbose = False, debug = False ):
    """
    Apply PBC's to create whole molecules

    Assumes molecule is shorter than box length

    """
    
    debug2 = False
    latticevec = strucC.getLatVec()

    n_dim = 3

    if( debug ):
        print "latticevec",latticevec
        F = open("shift.rec","w")
        
    for mol_i in range( 1,strucC.ptclC.n_molecules()+1) :

        if( verbose ):
            print "  Checking molecule %d of %d "%(mol_i,strucC.ptclC.n_molecules())

        searchD = {'chain':mol_i}
        mol_list = strucC.ptclC.getParticlesWithTags(searchD)

        if( len( mol_list) > 0 ):
            ptclObj = strucC.ptclC[1]
            r_o = ptclObj.position
            pid_o  = 1

            part_shifted = [False]*len(strucC.ptclC)

            r_mol_mass = np.zeros(n_dim)
            shift = np.zeros(n_dim)
            total_mass = 0.0 

            # shift all atoms to be conected 
            
            for pid_i in sorted(mol_list):
                ptclObj_i = strucC.ptclC[pid_i]
                a_mass_i = ptclObj_i.mass
                r_i = ptclObj_i.position

                r_io = np.array(r_o) - np.array(r_i)

                # sum center of mass
                total_mass += a_mass_i


                
                shifted = False 
                for dim in range(n_dim):
                    shift_dim = round( r_io[dim]/  latticevec[dim][dim] )
                    r_i[dim] = r_i[dim]  + latticevec[dim][dim] * shift_dim
                    if( shift_dim != 0 ):
                        shifted = True 
                        
                    r_mol_mass[dim] = r_mol_mass[dim]  + a_mass_i*r_i[dim] 
                

                if( debug and shifted ):
                    shift_line = "\n  %d - %d   r = %f %f %f dr = %f %f %f "%(pid_i,pid_o,ptclObj_i.position[0],ptclObj_i.position[1],ptclObj_i.position[2],r_io[0],r_io[1],r_io[2])
                    F.write(shift_line)
                
                ptclObj_i.position = r_i
                r_o = r_i
                pid_o = pid_i

            # Shift molecular center of mass into box 
            for dim in range(n_dim):
                cent_mass_i = r_mol_mass[dim] /total_mass
                shift[dim] = latticevec[dim][dim] * round( cent_mass_i /  latticevec[dim][dim] )

            
            for pid_i in sorted(mol_list):
                ptclObj_i = strucC.ptclC[pid_i]
                for dim in range(n_dim):
                    ptclObj_i.position[dim] = ptclObj_i.position[dim] - shift[dim] 
                
                

    if( debug ):
        F.close()


        
def main():
    """
    
    """
    prop_dim = 3
    
    prtCl_time = False
    debug = False
    debug2  = False 
    #
    # Formated ouput varables
    #
    sperator_line = "\n---------------------------------------------------------------------"

    # Initialize mpi
    p = mpiBase.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()
 
    options, args = get_options()

    if(  options.group_residues ):
        options.group_chains  = False
        
    if( rank == 0 ):
        # record initial time 
        t_i = datetime.datetime.now()
        # Open log files 
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w")

        dump_file = options.output_id + ".dump"
        dump_out = open(dump_file,"w")
        dump_line = "# g_i,g_j,res_i,res_j,chain_i,chain_j,p_i,p_j,r_mag  "
        dump_out.write(dump_line)



        vmd_file = options.output_id + ".vmd"
        vmd_out = open(vmd_file,"w")
        vmd_line = " g_i g_j  residue_i residue_j chain_i chain_j n_ij  cmass_i cmass_j , n_ij  "
        vmd_out.write(vmd_line)


        group_file = options.output_id + "_group.xyz"
        group_out = open(group_file,"w")
        group_out.close()
        
        #group_line = "# g_i,r_cmass[0],r_cmass[1],r_cmass[2],g_radius "
        #group_out.write(group_line)

        time_series = open(str("%s.time"%(options.output_id)),"w")
        time_series.write("# Averages per frame of all considered groups ")
        time_series.write("\n # frame,r_sq_gy_t,Rgyeigval_sq_gy_t[0],Rgyeigval_sq_gy_t[1],Rgyeigval_sq_gy_t[2],asphere_t,dlsq_t")

        # Clear xmol file 
        pi_xmol = open("{}.xyz".format(options.output_id),"w")
        pi_xmol.close()
        
    p.barrier()
    #
    #  Initialize blank system 
    #
    struc_o = StructureContainer()
    param_o = ParameterContainer()

    if( rank == 0 ):

        log_line = sperator_line
        log_line += "\n  Reading in input files "
        log_line += sperator_line

        log_out.write(log_line)
        print log_line

    cply_read = False
    if( len(options.in_cply) > 0 ): 
        if( rank == 0 and options.verbose ):
            print "  Reading cply reference {} ".format(options.in_cply)
        struc_o.read_cply(options.in_cply)
        cply_read = True
    if(  len(options.in_data) ):
        if( rank == 0 and options.verbose ): print "  Reading  LAMMPS data file {} ".format(options.in_data)
        read_lmpdata( struc_o , param_o , options.in_data, coordupdate=cply_read)
    if(  len(options.in_gro) ):
        if( rank == 0 and options.verbose ): print "  Reading  GROMACS gro file {} ".format(options.in_gro)
        read_gro(struc_o,options.in_gro, coordupdate=cply_read)
    else:
        options.in_gro = "{}.gro".format(options.output_id)
        if( rank == 0 ):
            if(  options.verbose ):
                print "  No .gro file read in will print a {} for md analysis read in ".format(options.in_gro)
            write_gro(struc_o,options.in_gro)
            
    p.barrier()
    if( len( struc_o.ptclC) == 0 ):
        error_line = " No input "
        sys.exit(error_line)

    #  Build nieghbor list 
    rerun_bonded = False
    
    if( len(struc_o.bondC) <= 0 or rerun_bonded ):
        struc_o.bondC.clear()
        struc_o.ptclC.guess_radii()
        struc_o.build_bonded_nblist(max_nn=12.0,radii_buffer=1.25)
        struc_o.nblist_bonds()
    else:
        struc_o.bondC_nblist()

    # If no rings selected add all non-zero values to search list
    if( len(options.ring) == 0 ):
        ring_list = []
        options.ring = ''
        for pid_i, ptcCl in struc_o.ptclC:
            ring_i = ptcCl.tagsDict["ring"]
            if( ring_i not in ring_list and ring_i > 0 ):
                ring_list.append(ring_i)
                options.ring += " {}".format(ring_i)

        if( rank == 0 ):
            log_line = " Ring option changed to {} ".format(options.ring)
            log_out.write(log_line)
            if( options.verbose ):
                print log_line
            
    # Print system properties
    if( rank == 0 ):
        print struc_o
        log_out.write(str(struc_o))
        
    #   
    # Filter particles
    #
    ptclC_o = struc_o.ptclC   
    search_o = dict()
    search_o = create_search(search_o,options.symbol,options.label,options.fftype,options.residue,options.resname,options.chain,options.ring)

    if( rank == 0 ):
        if( options.verbose ): print " Filter input by ",search_o

    list_f = ptclC_o.getParticlesWithTags(search_o)
    sum_f = len(list_f)
    
    if( rank == 0 and options.verbose ):
        print "  %d particles found in search "%(sum_f)    
    # 
    # Read in trajectory 
    #
    traj_read = True 
    if( len(options.in_xtc ) and  len(options.in_gro) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s and %s "%(options.in_xtc,options.in_gro)
        universe =  Universe(options.in_gro, options.in_xtc)
    elif( len(options.in_dcd ) and  len(options.in_gro) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s and %s "%(options.in_gro,options.in_dcd)
        universe =  Universe( options.in_gro , options.in_dcd)
    elif(len(options.in_gro) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s "%(options.in_gro)
        universe =  Universe(options.in_gro)
    else:
        traj_read = False 
        error_line = "No Trajectory read in  "
        #print error_line
        sys.exit(error_line)

    p.barrier()

    if( traj_read ):
        n_frames = len(universe.trajectory)
        if( rank == 0 ):

            log_line = " Trajector read in with %d frames "%(n_frames)
            log_out.write(log_line)
            if( options.verbose ):
                print log_line
            
    # Set all resnum to -1 so they wont be selected by default 
    for pid_i, ptcCl in struc_o.ptclC:
        universe.atoms[pid_i-1].resnum = -1
        
    # Find max chain and residue numbers
    chain_max = -1
    residue_max = -1
    for pid_i, ptcCl in struc_o.ptclC:
        if( ptcCl.tagsDict["chain"] > chain_max ): chain_max =  ptcCl.tagsDict["chain"]
        if( ptcCl.tagsDict["residue"] > residue_max ): residue_max =  ptcCl.tagsDict["residue"]

    chain_mult =  float( len( str( abs( round(chain_max,0) )))*10.0/len( str( abs( round(chain_max,0) ))) )*10.0

    # 
    # Rename residue numbers as group numbers
    #
    debug_groups= True 
    if( rank == 0 ):
        log_line = "\n Max chain # {}".format(chain_max)
        log_line += "\n Max residue # {}".format(residue_max)
        if( options.group_chains ):
            log_line += "\n  Grouping selection by chain  "
        elif( options.group_residues ):
             log_line += "\n  Grouping selection by chain*{} + residue  ".format(chain_mult)         
        log_out.write(log_line)
        if( options.verbose ):
            print log_line

    #
    # Group particles 
    #
    g_ind_list = []
    group_pid_list = []
    g_list = []
    g_i = -1 
    for pid_i, ptcCl in struc_o.ptclC:
        tagsDict=ptcCl.tagsDict
        tagsDict["group"] = -1
        ptcCl.setTagsDict(tagsDict)
        if( pid_i in list_f ):
            # Set unique index group_ind
            if( options.group_chains ):
                group_ind = int(ptcCl.tagsDict["chain"])  # Set resnum to chain number for inter/intra rdf's
            elif( options.group_residues ):
                group_ind =  int(ptcCl.tagsDict["chain"]*chain_mult + ptcCl.tagsDict["residue"])  # Set resnum to chain number for inter/intra rdf's
                
            if( group_ind not in g_ind_list ):
                g_i += 1
                g_ind_list.append(group_ind)
                group_pid_list.append([])
                g_list.append(g_i)
            else:
                g_i = g_ind_list.index(group_ind)
                
            ptcCl.tagsDict["group"] = g_i
            group_pid_list[g_i].append(pid_i)
            if( traj_read ):
                universe.atoms[pid_i-1].mass = ptcCl.tagsDict["mass"] 
                universe.atoms[pid_i-1].resnum = g_i  # Set resnum to chain number for inter/intra rdf's

    n_groups = len(g_list)
    # Print system properties
    if( rank == 0 ):
        log_line = "  Groups found {} ".format(n_groups)
        log_out.write(log_line)
        if( options.verbose ):
            print log_line
                        
                    
    debug = False 
    if( debug ):
        for pid_i, ptcCl in struc_o.ptclC:
            print pid_i,ptcCl.tagsDict["group"],ptcCl.tagsDict["chain"] ,ptcCl.tagsDict["residue"] 
        sys.exit(" grpu asdf")

    # Get paticle and bond structures
    ptclC_o = struc_o.ptclC
    bondC_o  = struc_o.bondC
    if( traj_read ):
        resname_sel = str("resname * ")
        un =  universe.selectAtoms(resname_sel)
        coord = un.coordinates()
        # dist = numpy.zeros((len(struc_o.ptclC),len(struc_o.ptclC)), dtype=numpy.float64) 

    #
    # If multi-core split the number of groups onto each core
    #
    split_list = True
    debug_split = False 
    if( split_list ):
        g_list_nodes  = p.splitListOnProcs(g_list)
        
        p.barrier()

        log_line = "Processor %d has %d groups "%(rank,len(g_list_nodes))
        if( options.verbose ):
            print log_line
                        
        #sys.exit(" split testing ")
        if( debug_split ):            
            for g_i in g_list_nodes:
                log_line = "Processor %d has group %s "%(rank,str(g_i))
                if( options.verbose ):
                    print log_line
                        
    else:
        g_list_nodes = g_list



    if( rank == 0 ):
        # record initial time 
        t_i = datetime.datetime.now()
        # Open log files 
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w")
        dat_file = options.output_id + ".dat"
        dat_out = open(dat_file,"w")
        dat_line = "#  Date "
        dat_line += "\n#  in_cply {} ".format(options.in_cply)
        dat_line += "\n#  in_data {} ".format(options.in_data)
        dat_line += "\n#  group_chains {} ".format(options.group_chains)
        dat_line += "\n#  group_residues {} ".format(options.group_residues)
        dat_line += "\n#  r_cut {} A".format(options.r_cut)
        dat_line += "\n#  r_e {} A".format(options.r_e)
        dat_line += "\n# g_i g_j  residue_i residue_j chain_i chain_j n_ij  cmass_i cmass_j , n_ij  "
        dat_out.write(dat_line)
        
    p.barrier()

    # Open output for time seriers
    for g_i in g_list_nodes:
        if( options.grp_time ):
            # Print time sereies for each molecule
            time_m = "%s_%d.time"%(options.output_id,g_i)
            time_out = open(time_m,"w")
            time_out.write("# frame,r_gy_sq,Rgyeigval_sq[0],Rgyeigval_sq[1],Rgyeigval_sq[2] " )
            time_out.close()

    # Save current directory as home
    home_dir = os.getcwd()

    calc_frames = 0
    volume_frames = []
    density_frames = []

    # LV = struc_o.getLatVec()
    # Loop over frames
    # Allocate distance matrix 
    # dist = numpy.zeros((n_i,n_j), dtype=numpy.float64)

    
    finial_frame = options.frame_f
    if( options.readall_f or  options.frame_f > len(universe.trajectory) ):
        finial_frame = len(universe.trajectory)

    p.barrier()

    if( traj_read ):        
        for ts in universe.trajectory:
            if( options.frame_o <= ts.frame ):
                if( ts.frame <= finial_frame  ):
                    if( ts.frame%options.frame_step == 0 ):
                        calc_frames += 1 
                        if( rank == 0 ):
                            t_f_i = datetime.datetime.now()
                            log_line = "\n Processing frame {}/{} time {}".format(ts.frame,finial_frame,t_f_i)
                            log_out.write(log_line)
                            print log_line
                        
                        volume_frames.append( ts.volume  )
                        box = ts.dimensions                    
                        latvec_i = struc_o.LCtoLV( box )
                        density_frames.append(struc_o.getDensity() )

                        coor_i = un.coordinates()
                        # Update particle positions
                        for pid_i, part_i in struc_o.ptclC:
                            part_i.position = coor_i[pid_i-1]
                        struc_pi =  copy.deepcopy(struc_o)
                        #
                        # Calculate center of mass of each group 
                        #                        
                        if( rank == 0 ):
                            log_line = "\n Calculating center of mass and radius for each group "
                            log_out.write(log_line)
                            if( options.verbose ):
                                print log_line
                                
                            
                        group_cmass_p = []
                        group_radius_p = [] 
                        for g_i in g_list_nodes:
                            total_mass = 0.0
                            r_cmass = np.array( [0.0,0.0,0.0] )
                            prop_dim = 3                            
                            for p_i in group_pid_list[g_i]:
                                part_i = struc_pi.ptclC[p_i]
                                a_mass = part_i.tagsDict["mass"]  
                                r_i = part_i.position 
                                # sum center of mass
                                total_mass += a_mass
                                for d in range(prop_dim):
                                    r_cmass[d] += a_mass*r_i[d]


                            for d in range(prop_dim):
                                r_cmass[d] = r_cmass[d]/total_mass
                            group_cmass_p.append(r_cmass)

                            g_radius = 0.0 # np.array( [0.0,0.0,0.0] )
                            
                            for p_i in group_pid_list[g_i]:
                                part_i = struc_pi.ptclC[p_i]
                                r_i =  np.array(part_i.position)
                                dr_cmass = r_i - r_cmass
                                mag_dr_cmass = np.linalg.norm(dr_cmass)
                                if( mag_dr_cmass > g_radius ):  g_radius = mag_dr_cmass
                                
                            group_radius_p.append(g_radius)


                        p.barrier()
			group_cmass  = p.gatherList(group_cmass_p)
			group_radius  = p.gatherList(group_radius_p)
                        p.barrier()
                        if( rank == 0 ):
                            if( rank ==0 ):
                                group_out = open(group_file,"a")
                                group_line = " {} ".format(len(g_list_nodes))
                                group_line += "\n# fake type r_cmass[0],r_cmass[1],r_cmass[2],g_i,g_radius  frame {} ".format(calc_frames)
                                group_out.write(group_line)
                                for g_i in g_list:
                                    r_cmass = group_cmass[g_i]
                                    g_radius = group_radius[g_i]
                                    group_line = "\n Li {}  {} {} {}   {}  ".format(r_cmass[0],r_cmass[1],r_cmass[2],g_i,g_radius)
                                    group_out.write(group_line)
                                group_out.close()
                                
                        p.barrier()
                        if( rank == 0 ):
                            log_line = "\n Adding electrons "
                            log_out.write(log_line)
                            if( options.verbose ):
                                print log_line
                        #
                        # Add electron positions
                        # 
                        add_pi_debug = False
                        # for g_i in g_list_nodes:
                        #     for p_i in group_pid_list[g_i]:
                        for p_i, part_i in struc_o.ptclC:
                            g_i = part_i.tagsDict["group"]
                            if(  g_i > -1 ):
                                part_i = struc_pi.ptclC[p_i]                        
                                struc_pi.add_pi(p_i,options.r_e)                                
                                group_pid_list[g_i].append(len(struc_pi.ptclC)-1)
                                group_pid_list[g_i].append(len(struc_pi.ptclC))
                                if( add_pi_debug ):
                                    print " Add pi e to particle {} in group {} ".format(p_i,g_i)
                                    print "  index {} or index {} or index {}".format(p_i-1,len(struc_pi.ptclC)-2,len(struc_pi.ptclC)-1)

                        p.barrier()
                        if( rank == 0 ):
                            log_line = "\n  Updating pi xmol file  "
                            log_out.write(log_line)
                            if( options.verbose ):
                                print log_line
                            comment = " From pi atoms from frame {} ".format(ts.frame)
                            append = True 
                            struc_pi.ptclC.write_xmol("{}_e.xyz".format(options.output_id),comment,append)
                            
                        if( add_pi_debug ):
                            for g_i in g_list_nodes:
                                for p_i in group_pid_list[g_i]:
                                    part_i = struc_pi.ptclC[p_i]
                                    print part_i.type,part_i.tagsDict["fftype"],part_i.tagsDict["label"],part_i.tagsDict["residue"],part_i.tagsDict["chain"], part_i.tagsDict["group"]
                            sys.exit(" add_pi_debug = True ")


                        #distance_array(coor_i,coor_i, box, result=dist)

                        # pioverlap_nblist = []
                        # pioverlap_nbindx = []

                        #for pid_i, ptclObj_i  in self.ptclC:
                        #    for pid_j, ptclObj_j  in self.ptclC:
                        #        if( pid_j != pid_i ):


                        debug_loop = False
                        if( debug_loop ):

                            for g_i in g_list_nodes[0:len(g_list_nodes)-1]:
                                for g_j in g_list_nodes[g_i+1:len(g_list_nodes)]:
                                    print g_i,g_j
                            sys.exit(" loop list debug ")
                            
                        if( rank == 0 ):
                            log_line = "\n Calculating connections "
                            log_out.write(log_line)
                            if( options.verbose ):
                                print log_line

                        debug_rij = False
                        debug_mpi = True
                        debug_mpi = True
                        dat_line_list_p = []
                        dump_line_list_p = []
                        vmd_line_list_p = []
                        e_bonded_list_p = []
                        group_nblist = []
                        p.barrier()
                        for g_i in g_list_nodes[0:len(g_list_nodes)-1]:
                            group_nblist.append([])
                            if( debug_mpi):
                                print rank," checking group ",g_i
                                
                            for g_j in g_list_nodes[g_i+1:len(g_list_nodes)]:
                                if( g_j !=  g_i ):
                                    cmass_i = group_cmass[g_i]
                                    radius_i = group_radius[g_i]
                                    cmass_j = group_cmass[g_j]
                                    radius_j = group_radius[g_j]
                                    
                                    nb_list = [] 
                                    for p_i in group_pid_list[g_i]:
                                        part_i = struc_pi.ptclC[p_i]
                                        if( part_i.type == "e" ):
                                            r_i = np.array( part_i.position )
                                            res_i = part_i.tagsDict["residue"]
                                            chain_i = part_i.tagsDict["chain"]
                                            cmass_i = group_cmass[g_i]                                    
                                            for p_j in group_pid_list[g_j]:
                                                part_j = struc_pi.ptclC[p_j]
                                                if( part_j.type == "e" ):
                                                    r_j = np.array( part_j.position )
                                                    res_j = part_j.tagsDict["residue"]
                                                    chain_j = part_j.tagsDict["chain"]
                                                    cmass_j = group_cmass[g_j]
                                                    r_mag = np.linalg.norm( struc_pi.delta_r_c(r_i,r_j) )
                                                    if( r_mag <= options.r_cut ):
                                                        e_bonded_list_p.append(p_i)
                                                        e_bonded_list_p.append(p_j)                                                        
                                                        nb_list.append(g_j)
                                                        dump_line_list_p.append("\n  {} {}  {} {}  {} {}  {} {} {}   ".format(g_i,g_j,res_i,res_j,chain_i,chain_j,p_i,p_j,r_mag))
                                                        vmd_line_list_p.append("\n  index {} or index {}  ".format(p_i-1,p_j-1))
                                                        if( debug_rij ):
                                                            print rank 
                                                            print "    {}-{} rij {}".format(p_i,p_j,r_mag)
                                                            print "    label", part_i.tagsDict["label"],part_j.tagsDict["label"]
                                                            print "    chain", part_i.tagsDict["chain"],part_j.tagsDict["chain"]
                                                            print "    residue",part_i.tagsDict["residue"],part_j.tagsDict["residue"]
                                                            print "    group", part_i.tagsDict["group"],part_j.tagsDict["group"]
                                                            print "    R_i ",r_i
                                                            print "    R_j ",r_j
                                    n_ij =    len(nb_list)
                                    # group_nblist.append([]
                                    # group_nblist[g_i].append(nb_list)
                                    # p.barrier()                                    
                                    if(  n_ij > 0 ):
                                        # dat_line = " g_i g_j  residue_i residue_j chain_i chain_j n_ij  cmass_i cmass_j , n_ij  "
                                        dat_line_list_p.append( "\n  {} {}  {} {}  {} {}  {}  ".format(g_i,g_j,res_i,res_j,chain_i,chain_j,n_ij))

                            #print g_i,nb_list
                            #for g_j in nb_list:

                        p.barrier() # Barrier for MPI_COMM_WORLD
                        
                        dat_line_list = p.gatherList(dat_line_list_p)
                        dump_line_list = p.gatherList(dump_line_list_p)
                        vmd_line_list = p.gatherList(vmd_line_list_p)
                        e_bonded_list = p.gatherList(e_bonded_list_p)
                        p.barrier()                                    


                        debug_mpi = True
                        if( debug_mpi):
                            print "rank,dat_line_list_p",rank,dat_line_list_p
                        
                        p.barrier()                                    
                        if( rank == 0 ):

                            if( debug_mpi):
                                print "rank,dat_line_list",rank,dat_line_list
                                
                            t_f_f = datetime.datetime.now()
                            log_line = "  time {} {} {} ".format(t_f_i,t_f_f,t_f_f-t_f_i)
                            print log_line
                            log_line = "\n Printing data files "
                            log_out.write(log_line)
                            print log_line
                            
                            for dat_line in dat_line_list:
                                dat_out.write(dat_line)
                            for dump_line in dump_line_list:
                                dump_out.write(dump_line)
                            for vmd_line in vmd_line_list:
                                vmd_out.write(vmd_line)

                        p.barrier()
                        if( rank == 0 ):
                            log_line = "\n  Updating pi xmol file  "
                            log_out.write(log_line)
                            if( options.verbose ):
                                print log_line
                            
                            # Update electrons to be bonded
                            for p_i in e_bonded_list:
                                part_i = struc_pi.ptclC[p_i]
                                part_i.type = "Xe"
                                
                            #
                            # Print new xyz file with bonded electrons marked as b
                            #
                            comment = " From pi atoms from frame {} ".format(ts.frame)
                            append = True 
                            struc_pi.ptclC.write_xmol("{}_ebond.xyz".format(options.output_id),comment,append)
                            
                            
            p.barrier() # Barrier for MPI_COMM_WORLD
            
    if( rank == 0 ):
        time_series.close()
        comment = " From pi atoms "
        append = False 
        struc_pi.ptclC.write_xmol("{}.xyz".format(options.output_id),comment,append)

    if( options.verbose and rank == 0 ):
	print "      Finding averages "
    #
    # Find averages
    #
    debug = 0
        
    if( rank == 0 ):
        log_out.close()
	dat_out.close()
	dump_out.close()
	vmd_out.close()
        group_out.close()

if __name__=="__main__":
    main()
   
