#! /usr/bin/env python
"""
Find distribution of dihedral angles

k-i-j-l

"""

# Dr. Travis Kemper
# Initial Date July 2014
# travis.kemper@nrel.gov

from structureContainer import StructureContainer
from parameters import ParameterContainer
import particles

import mpiBase
import datetime, sys
import numpy as np

import  units
from periodictable import periodictable

# Scott's new classes
import particles
from particles import Particle
from particles import ParticleContainer
from bonds     import BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer
from impropers     import Improper, ImproperContainer

from structureContainer import StructureContainer
from parameters    import ParameterContainer
from parameters    import ljtype,   LJtypesContainer
from parameters    import bondtype, BondtypesContainer
from parameters    import angletype,AngletypesContainer
from parameters    import dihtype,  DihtypesContainer
from parameters    import imptype,  ImptypesContainer

#MDanalysis
try:
    from MDAnalysis import *
    from MDAnalysis.core.distances import * ##distance_array
    #import MDAnalysis.core.units            # for bulk water density
except:
    import sys
    print "MDAnalysis module not build/configured correctly"
    sys.exit(0)

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    #
    # Input files
    #
    parser.add_option("--in_cply", dest="in_cply", type="string", default="", help="Input cply file")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_dcd", dest="in_dcd", type="string", default="", help="Input trajectory file in compressed dcd format ")
    parser.add_option("--in_xtc", dest="in_xtc", type="string", default="", help="Input xtc file with atoms listed as atom type numbers")  
    #
    #  Options
    #
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.05, help=" Bin size in degrees ")
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=False,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")
    parser.add_option("-o","--output_id", dest="output_id", default="trans",type="string",help=" prefix for output files  ")
    #
    # Searchable properties
    # k
    parser.add_option("--symbol_k", dest="symbol_k", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--label_k", dest="label_k", type="string", default="", help="select atoms of group i by  label id ")    
    parser.add_option("--fftype_k", dest="fftype_k", type="string", default="", help="select atoms of group i by force field type  ")
    parser.add_option("--residue_k", dest="residue_k", type="string", default="", help="select atoms of group i by resudue number  ")    
    parser.add_option("--resname_k", dest="resname_k", type="string", default="", help="select atoms of group i by residue name  ")    
    parser.add_option("--chains_k", dest="chains_k", type="string", default="", help="select atoms of group i by chain number  ")    
    parser.add_option("--rings_k", dest="rings_k", type="string", default="", help="select atoms of group i by ring number  ")   
    # I
    parser.add_option("--symbol_i", dest="symbol_i", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--label_i", dest="label_i", type="string", default="", help="select atoms of group i by  label id ")    
    parser.add_option("--fftype_i", dest="fftype_i", type="string", default="", help="select atoms of group i by force field type  ")
    parser.add_option("--residue_i", dest="residue_i", type="string", default="", help="select atoms of group i by resudue number  ")    
    parser.add_option("--resname_i", dest="resname_i", type="string", default="", help="select atoms of group i by residue name  ")    
    parser.add_option("--chains_i", dest="chains_i", type="string", default="", help="select atoms of group i by chain number  ")    
    parser.add_option("--rings_i", dest="rings_i", type="string", default="", help="select atoms of group i by ring number  ")    
    # J
    parser.add_option("--symbol_j", dest="symbol_j", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--label_j", dest="label_j", type="string", default="", help="select atoms of group j by  label id ")    
    parser.add_option("--fftype_j", dest="fftype_j", type="string", default="", help="select atoms of group j by force field type  ")
    parser.add_option("--residue_j", dest="residue_j", type="string", default="", help="select atoms of group j by resudue number  ")    
    parser.add_option("--resname_j", dest="resname_j", type="string", default="", help="select atoms of group j by residue name  ")    
    parser.add_option("--chains_j", dest="chains_j", type="string", default="", help="select atoms of group j by chain number  ")    
    parser.add_option("--rings_j", dest="rings_j", type="string", default="", help="select atoms of group j by ring number  ")    
    # l
    parser.add_option("--symbol_l", dest="symbol_l", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--label_l", dest="label_l", type="string", default="", help="select atoms of group i by  label id ")    
    parser.add_option("--fftype_l", dest="fftype_l", type="string", default="", help="select atoms of group i by force field type  ")
    parser.add_option("--residue_l", dest="residue_l", type="string", default="", help="select atoms of group i by resudue number  ")    
    parser.add_option("--resname_l", dest="resname_l", type="string", default="", help="select atoms of group i by residue name  ")    
    parser.add_option("--chains_l", dest="chains_l", type="string", default="", help="select atoms of group i by chain number  ")    
    parser.add_option("--rings_l", dest="rings_l", type="string", default="", help="select atoms of group i by ring number  ")   
    #
    (options, args) = parser.parse_args()
        
    return options, args
   

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
        #gro_lines += "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"  % (pt_i.tagsDict["residue"],pt_i.tagsDict["resname"][:5],pt_i.tagsDict["fftype"][:5],atom_indx,r_i_nm[0],r_i_nm[1],r_i_nm[2] )
        if( atom_indx > 99999 ):
            atom_indx = 1
            
    gro_lines += " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (units.convert_angstroms_nm(latvec[0][0]),units.convert_angstroms_nm(latvec[1][1]),units.convert_angstroms_nm(latvec[2][2]),units.convert_angstroms_nm(latvec[0][1]),units.convert_angstroms_nm(latvec[0][2]),units.convert_angstroms_nm(latvec[1][0]),units.convert_angstroms_nm(latvec[1][2]),units.convert_angstroms_nm(latvec[2][0]),units.convert_angstroms_nm(latvec[2][1])) 

    F = open( data_file, 'w' )
    F.write(gro_lines)
    F.close()
                    
def calc_bins(val_binsize,val_min,val_max):
    
    val_floor = int(val_min/val_binsize)*val_binsize - 2.0*val_binsize
    val_ceil = int(val_max/val_binsize)*val_binsize + 2.0*val_binsize
    val_range = val_ceil - val_floor
    val_bins = int(val_range/val_binsize) + 1

    return val_floor,val_ceil,val_bins

def dih_hist(angle_list,struc_o,bin_size,hist_cnt):
    """
    Loop over dihedrals and calculate histogram 
    """
    
    debug = False
    
    #
    dih_cnt = 0

    for indx_kij in angle_list:
        dih_cnt += 1 
        a_k = indx_kij[0]
        a_i = indx_kij[1]
        a_j = indx_kij[2]
        a_l = indx_kij[3]

        angle_i = struc_o.getDihedral(a_k,a_i,a_j,a_l)
        abs_angle_i = np.absolute(angle_i)

        #if( abs_angle_i > 180.0 - bin_size/2.0 ):
        #    if( debug ):
        #        print " rescaling angle ",abs_angle_i," to ", abs_angle_i  -180.0 
        #    abs_angle_i += -180.0
        if( abs_angle_i <= 180.0 ):
            # Acount for round off errors 
            abs_angle_i += -0.0001
            
        bin_index = int(  abs_angle_i/bin_size ) 
        hist_cnt[bin_index] += 1

        
        if( debug ):
            hist_val = bin_size*float(bin_index) + bin_size/2.0
            if( bin_index == 0 ): hist_val = 0.0 
            if( hist_val == 179.0 ): hist_val = 180.0 
            print " index %d or index %d or index %d or index %d  %f bin_index %d hist_cnt %d hist_val %f "%(a_k-1,a_i-1,a_j-1,a_l-1,angle_i,bin_index, hist_cnt[bin_index],hist_val)

    return hist_cnt


def build_covnablist(strucC,debug = False ):
    """
    Build covalent neighbor list from elements and positions 
    """

    n_ptcl = len( strucC.ptclC )
    maxnnab = n_ptcl*12            # Assume max nieghbors as fcc 

    cov_buffer = 1.25
    
    radi_cov =  []
    cov_nblist = np.empty( maxnnab,  dtype=int )
    cov_nbindx = np.empty( maxnnab,  dtype=int )

    NNAB = 0


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

def get_dihatoms(struc_o,list_k,list_i,list_j,list_l,verbose = True,debug = False ):
    """
    Find sets of dihedrals in system

    k-i-j-l

    Arguments
        list_k (list) of atom indexes in of the first bonded atom in the dihedral 
        list_i (list) of atom indexes in of the second bonded atom in the dihedral 
        list_j (list) of atom indexes in of the third bonded atom in the dihedral 
        list_l (list) of atom indexes in of the fourth bonded atom in the dihedral
    Return
        angle_list (list) of four aotms in each dihedral 
    """

    import datetime
    # import topology

    
    
    #
    if(debug):
        print " list_k ",list_k
        print " list_i ",list_i
        print " list_j ",list_j
        print " list_l ",list_l
    #
    # Create neighbor list form bonds
    #
    if( len(struc_o.bondC) > 0 ):
        if( verbose ):
            log_line = "  Building neighbor list from {} bonds".format(len(struc_o.bondC))
            print log_line
        struc_o.bondC_nblist()
    else:
        if( verbose ): print " No bond found creating bonded niehgborlist based on radii "
        struc_o.build_bonded_nblist(max_nn=12.0,radii_buffer=1.25)
        struc_o.nblist_bonds()
    
    cov_nblist = struc_o.bonded_nblist
    cov_nbindx = struc_o.bonded_nbindx
    #
    # Find  atom groups k-i-j
    #
    angle_list = []
    #
    sum_angles = 0
    #
    # Find atom indices  of group i and j
    #
    #for p_k, ptcl_k  in struc_o.ptclC(list_k):
    for atom_k  in list_k:
        
        N_k_o = cov_nbindx[atom_k]
        N_k_f = cov_nbindx[atom_k+1]

        if(debug): print  "checking k ",atom_k,struc_o.ptclC[atom_k].type," with ",N_k_f - N_k_o," nbs"

        for indx_i in range( N_k_o,N_k_f):
            atom_i = cov_nblist[indx_i]
            add_i = False

            if( debug ): print " checking i ",atom_i,struc_o.ptclC[atom_i].type

            for  p_i in list_i:
                if( atom_i == p_i ): #and atom_i > atom_k ):
                    add_i = True
            if( add_i ): #atom_i in list_i ):


                N_i_o = cov_nbindx[atom_i]
                N_i_f = cov_nbindx[atom_i+1] 

                for indx_j in range( N_i_o,N_i_f):
                    atom_j = cov_nblist[indx_j]
                    add_j = False

                    if( debug ): print " checking j ",atom_j,struc_o.ptclC[atom_j].type

                    for  p_j  in list_j:
                        if( atom_j == p_j ): #and atom_j > atom_i ):
                            add_j = True
                    if( add_j ): # atom_j  in list_j  ):

                        N_j_o = cov_nbindx[atom_j]
                        N_j_f = cov_nbindx[atom_j+1] 

                        for indx_l in range( N_j_o,N_j_f):
                            atom_l = cov_nblist[indx_l]
                            add_l = False

                            if( debug ):
                                print " checking l ",atom_l,struc_o.ptclC[atom_l].type                                    
                            for  p_l  in list_l:
                                if( atom_l == p_l ): #and atom_l > atom_j ):
                                    add_l = True
                            if( add_l ): #atom_l  in list_l ):
                                # Check to make sure not in list already
                                add_dih = True
                                for indx_kij in angle_list:
                                    a_k = indx_kij[0]                
                                    a_i = indx_kij[1]
                                    a_j = indx_kij[2]
                                    a_l = indx_kij[3]
                                    if( atom_k == a_k and atom_i == a_i and atom_j == a_j and atom_l == a_l ):
                                        add_dih = False 
                                    if( atom_k == a_l and atom_i == a_j and atom_j == a_i and atom_l == a_k ):
                                        add_dih = False 
                                if(add_dih ):
                                    angle_list.append( [atom_k,atom_i,atom_j,atom_l] )

                                    if(debug):

                                        t_f = datetime.datetime.now()
                                        dt_sec  = t_f.second - t_i.second
                                        dt_min  = t_f.minute - t_i.minute
                                        if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0         
                                        print " found angle ",atom_k,atom_i,atom_j,atom_l
                                        print "  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "


    if ( debug ):
        dih_cnt = 0 
        for indx_kij in angle_list:
            dih_cnt += 1 
            a_k = indx_kij[0]                
            a_i = indx_kij[1]
            a_j = indx_kij[2]
            a_l = indx_kij[3]
            print dih_cnt," found angle ",a_k,a_i,a_j,a_l

        sys.exit('get_dihatoms debug')

    return angle_list


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



def norm_r_ij(r_i,r_j,latticevec):
    """
    Normailze difference between two vectors
    using cubic periodic boundry conditions 
    
    Args:
         r_i (vector numpy[3] ) position of particle i 
         r_j (vector numpy[3] ) position of particle j
         latticevec (numpy[3][3] ) 3 lattice vectors 
    Returns:
         sq_dr (float) square of vector from r_i to r_j
         
    """

    debug = False
    
    delta_ij = delta_r_c(r_i,r_j,latticevec)

    if( debug):
        print "delta_ij ",delta_ij
        print " mag ",np.linalg.norm(delta_ij)
    
    return (delta_ij)/np.linalg.norm(delta_ij)



def getVolume_c(latvec):
    """
    Calculate volume of Orthorhombic unit cell 
    Method:
    none cubic Volume = ( v_i x v_j ) . v_k / cubic Volume =v_i  v_j  v_k
    """
    
    vol = latvec[0][0]*latvec[1][1]*latvec[2][2]
    
    return vol


def getDensity(total_mass_i,latvec_i):

    """
    Calculate density of system in AMU/A^3 and convert to g/cm^3
    NOTE: mass units contained in PtclConatiner
    """

    volume_i = getVolume_c(latvec_i)
    density_i = units.convert_AMUA3_gcm3(total_mass_i/volume_i) 

    return density_i
        
def dihdist():
    """
    Read in files and create new files

    Arguments

    Return
    
    """

    debug = False 
    
    #
    # Formated ouput varables
    #
    sperator_line = "\n---------------------------------------------------------------------"
    

    # Initialize mpi
    p = mpiBase.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    #        
    # Read options 
    #
    options, args = get_options()
    #
    #  Initialize blank system 
    # 
    struc_o = StructureContainer(verbose=False)
    param_o = ParameterContainer()
    #
    # Open output files 
    #
    if( rank == 0 ):

        log_file = options.output_id + ".log"
        log_out = open(log_file,"w") 

        dat_file = options.output_id + ".dat"
        dat_out = open(dat_file,"w") 
        dat_out.write("#   frame, cos_kijl, a_k,a_i,a_j,a_l,chain_i,residue_i,residue_j, cos_kij, cos_ijl ")

    p.barrier()
    #
    # Read in system data 
    #
    if( rank == 0 ):
        log_line = sperator_line
        log_line += "\n  Reading in input files "
        log_line += sperator_line
        log_out.write(log_line)
        print log_line
    #
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
    if( len( struc_o.ptclC) == 0 ):
        error_line = " No input "
        sys.exit(error_line)
    
    # Get paticle and bond structures
    ptclC_o = struc_o.ptclC
    bondC_o  = struc_o.bondC
    
    p.barrier()

    # Print system properties
    if( rank == 0 ):
        print struc_o
        log_out.write(str(struc_o))

    p.barrier()
    
    hist_min = -1.0
    
    # Print system properties
    if( rank == 0 ):
        print struc_o
        log_out.write(str(struc_o))




    #   
    # Create sub lists for dihdral groups k-i-j-l
    #
    if( rank == 0 and options.verbose  ):
        print "  searching for list k i j l "
        
    search_k = dict()
    search_k = create_search(search_k,options.symbol_k,options.label_k,options.fftype_k,options.residue_k,options.resname_k,options.chains_k,options.rings_k)
    if( rank == 0 ):
        log_line = "\n Searching group _k {} ".format(search_k)
        log_out.write(log_line)
        if( options.verbose ): print log_line
    list_k = struc_o.ptclC.getParticlesWithTags(search_k)
    sum_k = len(list_k)
    if( rank == 0 ):
        log_line = "\n  group _k has {} particles ".format(sum_k)
        log_out.write(log_line)
        if( options.verbose ): print log_line

    

    search_i = dict()
    search_i = create_search(search_i,options.symbol_i,options.label_i,options.fftype_i,options.residue_i,options.resname_i,options.chains_i,options.rings_i)
    if( rank == 0 ):
        log_line = "\n Searching group i {} ".format(search_i)
        log_out.write(log_line)
        if( options.verbose ): print log_line
    list_i = struc_o.ptclC.getParticlesWithTags(search_i)
    sum_i = len(list_i)
    if( rank == 0 ):
        log_line = "\n  group i has {} particles ".format(sum_i)
        log_out.write(log_line)
        if( options.verbose ): print log_line


    search_j = dict()
    search_j = create_search(search_j,options.symbol_j,options.label_j,options.fftype_j,options.residue_j,options.resname_j,options.chains_j,options.rings_j)
    if( rank == 0 ):
        log_line = "\n Searching group i {} ".format(search_j)
        log_out.write(log_line)
        if( options.verbose ): print log_line
    list_j = struc_o.ptclC.getParticlesWithTags(search_j)
    sum_j = len(list_j)
    if( rank == 0 ):
        log_line = "\n  group i has {} particles ".format(sum_j)
        log_out.write(log_line)
        if( options.verbose ): print log_line



    search_l = dict()
    search_l = create_search(search_l,options.symbol_l,options.label_l,options.fftype_l,options.residue_l,options.resname_l,options.chains_l,options.rings_l)
    if( rank == 0 ):
        log_line = "\n Searching group i {} ".format(search_l)
        log_out.write(log_line)
        if( options.verbose ): print log_line
    list_l = struc_o.ptclC.getParticlesWithTags(search_l)
    sum_l = len(list_l)
    if( rank == 0 ):
        log_line = "\n  group i has {} particles ".format(sum_l)
        log_out.write(log_line)
        if( options.verbose ): print log_line



    if( rank == 0 and options.verbose  ):
        print "  searching finished "
        print "  Creating list of atomic indies of dihedrals "

    # Create list of dihderal angles 
    angle_list = get_dihatoms(struc_o,list_k,list_i,list_j,list_l)
    
    if( len(angle_list) == 0 ):
        error_line = " No dihdral angles found for selected atoms "
        sys.exit(error_line)
        
    if( rank == 0 ):
        log_line = " %d angles found "%(len(angle_list))
        log_out.write(log_line+"\n")
        if( options.verbose ):
            print log_line
    

    w_all_angles = False 
    if( rank == 0 and w_all_angles ):
        log_line = "Indecies of dihedral atoms to be histogramed"
        log_out.write(log_line+"\n")
        if( options.verbose ):
            print log_line
        for angle_i in angle_list:
            log_line = str(angle_i)
            log_out.write(log_line+"\n")
            if( options.verbose ):
                print log_line
    #
    # If multi-core split the number of angles onto each core
    #
    split_list = True
    if( split_list ):
        myChunk_i  = p.splitListOnProcs(angle_list)

        #
        #for angle_i in myChunk_i:
        log_line = "Processor %d has %s of %d "%(rank,len(myChunk_i),len(angle_list))
        if( rank == 0 ):
            log_out.write(log_line+"\n")
        if( options.verbose ):
            print log_line
                        
    else:
        myChunk_i = angle_list
    
    if( rank == 0 and options.verbose  ):
        print "  list creation finished "

    # Save lattice vector 
    latvec_i = struc_o.latvec 
    total_mass_i = struc_o.getTotMass() 

    chain_n = []
    residue_n = []
    resname_n = []
    for pid_i, ptclObj_i  in struc_o.ptclC:
        chain_n.append(ptclObj_i.tagsDict["chain"] )
        residue_n.append(ptclObj_i.tagsDict["residue"] )
        resname_n.append(ptclObj_i.tagsDict["resname"] )
                       
    
    
    # Delete StructureContainer to free memory
    del struc_o
    del param_o
    # 
    # Read in trajectory 
    #
    if( len(options.in_dcd) and len(options.in_gro) ):
        if(  rank == 0 and options.verbose ): print "  Reading  %s and %s "%(options.in_gro,options.in_dcd)
        universe =  Universe( options.in_gro , options.in_dcd)
    elif( len(options.in_xtc) and len(options.in_gro) ):
        if( rank == 0 and options.verbose ): print "  Reading  %s and %s "%(options.in_gro,options.in_xtc)
        universe =  Universe( options.in_gro , options.in_xtc)
    elif( len(options.in_gro) ):
        if(  rank == 0 and options.verbose ): print "  Reading  %s"%(options.in_gro)
        universe =  Universe( options.in_gro )
    else:
        error_line = " No input "
        sys.exit(error_line)
        
    p.barrier()

    if( rank == 0 and options.verbose ):
        print " Trajector read in with %d frames "%(len(universe.trajectory))

    p.barrier()
        

    # Select all atoms from the universe 
    uni_i = universe.selectAtoms(" resname * ")

    # Relabel 
    for pid_i in range(len(uni_i)):
        universe.atoms[pid_i].segid = str(chain_n[pid_i])
        universe.atoms[pid_i].resid = str(residue_n[pid_i])

        if( debug):

            print "resid", universe.atoms[pid_i].resid
            print "residue", universe.atoms[pid_i].residue
            print "segment",universe.atoms[pid_i].segment
            print "segid", universe.atoms[pid_i].segid

            #sys.exit(" albel debug")
    
    p.barrier()
    
    finial_frame = options.frame_f
    if( options.readall_f ):
        finial_frame = len(universe.trajectory)

    hist_cos = False 
    if( hist_cos ):

        #
        # Intialize lsits and counts
        #

        cos_max = 1.0
        cos_min = -1.0
        cos_floor,cos_ceil,cos_bins = calc_bins(options.bin_size,cos_min,cos_max)

        hist_cnt = np.zeros(cos_bins+1)
        hist_cnt_node = np.zeros(cos_bins+1)

        if( rank == 0 ):
            
            hist_file = options.output_id + ".hist"
            hist_out = open(hist_file,"w")
            hist_line  = "#   Histogram  \n"
            hist_line += "#   N bins %d  \n"%(cos_bins)
            hist_out.write(hist_line)
            
    if( rank == 0 ):

        # log_line += "\n  %d "%()
        log_line  = "" #"\n  Date %s "%(str(t_i))
        log_line += "\n  Number of processors  %d "%(size)

        log_line += sperator_line
        log_line += "\n  Particles of group l: %d "%(sum_l)
        log_line += "\n  Particles of group i: %d "%(sum_i)
        log_line += "\n  Particles of group j: %d "%(sum_j)
        log_line += "\n  Particles of group k: %d "%(sum_k)

        log_line += sperator_line
        log_line += "\n Frames "
        log_line += "\n     Initial frame  %d "%(options.frame_o)
        log_line += "\n     Step frame  %d  "%(options.frame_step)
        log_line += "\n     Final frame  %d  "%(finial_frame)
        log_line += sperator_line
        if( hist_cos ):

            log_line += "\n Histogram "
            log_line += "\n     Min  %f "%(cos_floor)
            log_line += "\n     Max  %f  "%(cos_ceil)
            log_line += "\n     Size  %f  "%(options.bin_size)
            log_line += "\n     N bins  %f  "%(cos_bins)

        log_out.write(log_line)
        if( options.verbose ):
            print log_line
            
    p.barrier()
    #
    # Initialize line count
    #
    frame_cnt = 0
    calc_frames = 0
    total_dih_cnt = 0 
    volume_frames = []
    density_frames = []
    cos_list = []
    cos_list_node = []
    cis_cnt = 0 
    trans_cnt = 0 
    conj_cnt = 0
    conj_cut = 0.9

    log_line += "#       Conj (<%f) %d %f \n"%( np.rad2deg( np.arccos(float(conj_cut))),conj_cnt,100.0*float(conj_cnt)/float(100) )
    if( options.verbose ):
        print log_line
    
    # 
    calc_frames = 0 
    for ts in universe.trajectory:
        if( options.frame_o <= ts.frame ):
            if( ts.frame <= finial_frame  ):
                if( ts.frame%options.frame_step == 0 ):
                        calc_frames += 1 
                        volume_frames.append( ts.volume  )
                        box = ts.dimensions
                        coor_i = uni_i.coordinates()
                        # Update latvec
                        for dim in range(3):
                            latvec_i[dim][dim] = box[dim]
                        density_frames.append( getDensity(total_mass_i,latvec_i) )
                        #
                        dih_cnt = 0
                        for indx_kij in myChunk_i:
                            dih_cnt += 1
                            total_dih_cnt += 1 
                            a_k = indx_kij[0]
                            a_i = indx_kij[1]
                            a_j = indx_kij[2]
                            a_l = indx_kij[3]

                            a_k_mda = a_k -1 
                            a_i_mda = a_i -1 
                            a_j_mda = a_j -1 
                            a_l_mda = a_l -1

                            chain_k = universe.atoms[a_k_mda].segid
                            chain_i = universe.atoms[a_i_mda].segid
                            chain_j = universe.atoms[a_j_mda].segid
                            chain_l = universe.atoms[a_l_mda].segid

                            residue_k = universe.atoms[a_k_mda].resid
                            residue_i = universe.atoms[a_i_mda].resid
                            residue_j = universe.atoms[a_j_mda].resid
                            residue_l = universe.atoms[a_l_mda].resid

                            r_k = np.array( coor_i[a_k_mda] )
                            r_i = np.array(coor_i[a_i_mda] )
                            r_j = np.array( coor_i[a_j_mda] )
                            r_l = np.array(coor_i[a_l_mda]  )

                            v1 = norm_r_ij(r_k, r_i,latvec_i)
                            v2 = norm_r_ij(r_i, r_j,latvec_i)
                            v3 = norm_r_ij(r_j, r_l,latvec_i)

                            v1v2 = np.cross(v1,v2)
                            v2v3 = np.cross(v2,v3)

                            debug_angle = False
                            if( debug_angle ):

                                print "v1",v1
                                print "v2",v2
                                print "v3",v3
                            
                            r_i_norm = v1v2/np.linalg.norm(v1v2)
                            r_j_norm = v2v3/np.linalg.norm(v2v3) 
                            cos_kijl = np.dot(r_i_norm,r_j_norm)

                            aux_angles = True 
                            if( aux_angles ):                                
                                cos_kij = np.dot(v1,v2)
                                cos_ijl = np.dot(v2,v3)


                            if( cos_kijl < 0.0 ):
                                trans_cnt += 1
                            else:
                                cis_cnt += 1
                            if(  np.absolute( cos_kijl ) <  conj_cut  ):
                                conj_cnt += 1
                            
                            if( hist_cos ):
                                bin_index = int(  ( cos_kijl- hist_min )/options.bin_size  )
                                hist_cnt_node[bin_index] += 1

                            debug_cos = False  
                            if( debug_cos ):
                                if( cos_kijl >= 1.0 ):
                                    ang_deg = 0.0
                                elif(  cos_kijl <= -1.0 ):
                                    ang_deg = 180.0
                                else:    
                                    cos_ang = np.arccos(cos_kijl )
                                    ang_deg = np.rad2deg( cos_ang )
                                print  "cos_kijl ",cos_kijl,ang_deg
                                
                            calc_angle = False
                            if( calc_angle ):

                                angle_i = getAngle(v1v2,v2v3)

                                if( debug):
                                    print "v1 ",v1
                                    print "v2 ",v2
                                    print "v3 ",v3
                                    print "v1v2 ",v1v2
                                    print "v2v3 ",v2v3
                                    print "angle_i ",angle_i
                                debug2 = False
                                if( debug2 ):
                                    print "angle_i ",a_k,a_i,a_j,a_l,angle_i
                                    print "index ",a_k_mda," or index ",a_i_mda," or index ",a_j_mda," or index ",a_l_mda


                                #
                                # Find sign of angle 
                                #
                                v1v3 = np.cross(v1,v3)
                                sign_v = np.dot(v2,v1v3)

                                if( sign_v > 0.0  ):
                                    angle_i = -1.0*angle_i

                                #hist_cnt = dih_hist(angle_list,struc_o,options.bin_size,hist_cnt)
                                abs_angle_i = np.absolute(angle_i)

                                if( abs_angle_i <= 180.0 ):
                                    # Acount for round off errors 
                                    abs_angle_i += -0.0001

                                # bin_index = int(  abs_angle_i/options.bin_size )
                                # print bin_index,abs_angle_i,options.bin_size
                                # hist_cnt_node[bin_index] += 1
                                if( rank == 0 and options.verbose ):
                                    log_line = " %d %f   %f   %d %d %d %d"%(ts.frame, cos_kijl,abs_angle_i, a_k,a_i,a_j,a_l )
                                    print log_line
                                    
                            #dat_line = " %d %f   %d %d %d %d"%(ts.frame, cos_kijl, a_k,a_i,a_j,a_l )
                            dat_line = "\n %d %f   %d %d %d %d  %s  %s %s"%(ts.frame, cos_kijl, a_k,a_i,a_j,a_l,chain_i,residue_i,residue_j )
                            if( aux_angles ):
                                dat_line = "\n %d %f   %d %d %d %d  %s  %s %s  %f %f "%(ts.frame, cos_kijl, a_k,a_i,a_j,a_l,chain_i,residue_i,residue_j, cos_kij, cos_ijl )
                            dat_out.write(dat_line)
                            cos_list.append(cos_kijl)


                        if( rank == 0 ):
                            log_line =  "Frame %4d with volume %f " % (ts.frame, ts.volume)
                            log_out.write(log_line+"\n")
                            if( options.verbose ): 
                                print log_line



    if( hist_cos ):
        for bin_index in range( 0,cos_bins+1):
            # Sum hist_cnt of each bin on each processor 
            hist_cnt[bin_index] = p.allReduceSum(hist_cnt_node[bin_index])

            # print "hist_cnt[bin_index] ",bin_index,hist_cnt[bin_index] 
        
    p.barrier() # Barrier for MPI_COMM_WORLD
    
    dat_out.close()
    log_out.write("\n")

    if( rank == 0 ):
        box_vol_ave = np.average( volume_frames )
        box_den_ave = np.average( density_frames )

        d_vol =  volume_frames[-1] - volume_frames[0]
        d_den = density_frames[-1] - density_frames[0]

        log_line  = "#    Frames %d \n" %  (calc_frames)
        log_line += "#       Average Box Volume %f A^3 \n" % (box_vol_ave) 
        log_line += "#       Delta Box Volume %f A^3 \n" % (d_vol) 
        log_line += "#       Average Box density %f g/cm^3 \n" % (box_den_ave) 
        log_line += "#       Delta Box density %f  g/cm^3 \n" % (d_den)
        log_line += "#    Total_cnts %d  \n" % (total_dih_cnt)
        log_line += "#       Trans %d %f \n"%(trans_cnt,100.0*float(trans_cnt)/float(len(cos_list)))
        log_line += "#       Cis %d %f \n"%(cis_cnt,100.0*float(cis_cnt)/float(len(cos_list)))
        log_line += "#       Conj (<%f) %d %f \n"%( np.rad2deg( np.arccos(float(conj_cut))),conj_cnt,100.0*float(conj_cnt)/float(len(cos_list)) )
        log_line += "\n# "
        if( options.verbose ):
            print log_line
        log_out.write(log_line)
        
        
        if( hist_cos ):

            hist_sum = np.sum(hist_cnt)

            # Write output 
            #
            hist_line = "#    Frames %d " %  (calc_frames)
            hist_line += "\n#    Bin-size %f  " % (options.bin_size)
            hist_line += "\n#    Total_cnts %d  " % (hist_sum)
            hist_line += "\n#    Total_cnts %d  " % (hist_sum)
            hist_line += "\n#    Average Box Volume %f A^3 " % ( box_vol_ave) 
            hist_line += "\n#    Delta Box Volume %f A^3 " % ( volume_frames[-1] - volume_frames[0] ) 
            hist_line += "\n#    Average Box density %f g/cm^3" % ( box_den_ave) 
            hist_line += "\n#    Delta Box density %f  g/cm^3 " % ( density_frames[-1] - density_frames[0] ) 
            hist_line += "\n# "
            hist_line += "\n# bin index ; cnt    ; cnt/Total_cnts  "
            hist_out.write(hist_line)

            if( options.verbose ):
                print hist_line

            if( debug ):
                sum_check = 0 
            for bin_index in range( 0,cos_bins+1):
                val_cnt =  hist_cnt[bin_index]

                if( val_cnt > 0.0 ):

                    hist_val = options.bin_size*float(bin_index) + hist_min + options.bin_size/2.0
                    bin_index_check = int(  ( hist_val- hist_min )/options.bin_size  )
                    if( bin_index_check != bin_index ):
                        error_line = "Bin interptratation incorrect for valeue %f "%(hist_val)
                        error_line += "\n %f != %f "%(  bin_index_check , bin_index)
                        sys.exit(error_line)


                        #if( bin_index == 0 ): hist_val = 0.0
                        #if( hist_val == 179.0 ): hist_val = 180.0 

                    #print " hist_cnt[bin_index] ",bin_index, hist_cnt[bin_index]

                    cnt_fnorm =   float(val_cnt)/float(hist_sum) #float(calc_frames)


                    if( hist_val >= 1.0 ):
                        ang_deg = 0.0
                    elif(  hist_val <= -1.0 ):
                        ang_deg = 180.0
                    else:    
                        cos_ang = np.arccos(hist_val )
                        ang_deg = np.rad2deg( cos_ang )

                    hist_out.write("\n  %d %f %f %f %f " % (bin_index,hist_val,val_cnt,cnt_fnorm,ang_deg) )

                    if( debug ): sum_check += val_cnt

            if( debug ): print "sum_check ",sum_check
        
        
        log_out.close()
        if( hist_cos ): hist_out.close()
    
	
	    
if __name__=="__main__":
    dihdist()
   
	    
