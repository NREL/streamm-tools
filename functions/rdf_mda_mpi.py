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

from itertools import izip
import json, math , sys
import numpy as np
import datetime


#MDanalysis
try:
    from MDAnalysis import *
    from MDAnalysis.core.distances import * ##distance_array
    #import MDAnalysis.core.units            # for bulk water density
except:
    print "MDAnalysis module not build/configured correctly"
    sys.exit(0)

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

import mpiBase
import units
from periodictable import periodictable


def get_options():
    """
    Set options
    """
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-o","--output_id", dest="output_id", default="rdf",type="string",help=" prefix for output files  ")
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
    parser.add_option("--mol_inter",dest="mol_inter", default=False,action="store_true", help="Use only inter molecular rdf's ")
    parser.add_option("--mol_intra",dest="mol_intra", default=False,action="store_true", help="Use only intra molecular rdf's")
    parser.add_option("--exclude14",dest="exclude14", default=False,action="store_true", help="Do not include 1-4 bonded neighbors particles in rdf")
    parser.add_option("--truedensity",dest="truedensity", default=False,action="store_true", help="Use the true density of group j, this makes inter/intra molecular rdfs not components of the total rdf but true independent rdfs")
    
    # Bins
    parser.add_option("--r_cut", dest="r_cut", type=float, default=20.0, help=" Cut off radius in angstroms ")
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.10, help=" Bin size in angstroms")
    #
    # Searchable properties
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
    # Frames
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=True,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")
    parser.add_option("--set_chaintoresidue", dest="set_chaintoresidue", default=False,action="store_true", help=" set_chaintoresidue during .gro file read in with no cply file")
        
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
        if( atom_indx > 99999 ):
            atom_indx = 1
            
    gro_lines += " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (units.convert_angstroms_nm(latvec[0][0]),units.convert_angstroms_nm(latvec[1][1]),units.convert_angstroms_nm(latvec[2][2]),units.convert_angstroms_nm(latvec[0][1]),units.convert_angstroms_nm(latvec[0][2]),units.convert_angstroms_nm(latvec[1][0]),units.convert_angstroms_nm(latvec[1][2]),units.convert_angstroms_nm(latvec[2][0]),units.convert_angstroms_nm(latvec[2][1])) 

    F = open( data_file, 'w' )
    F.write(gro_lines)
    F.close()
                    
def main(debug=False):
    """
    Calculate radial distribution (RDF) based on specified atom types

    Input:
        - cply file ".cply" containing all the information for structure
        - frames to average RDF over
        - calculation specifications 
            - atomic groups to calculate RDF between
            - which frames to use 
    Output:
        - log file ".log" containing calculation information
            the amount of information in the log file can be increased using the verbose option ( -v )         
    """
    #
    # Formated ouput varables
    #
    sperator_line = "\n---------------------------------------------------------------------"
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    #
    options, args = get_options()
    #
    if( rank == 0 ):
        # Open log files 
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w") 
    p.barrier()

    # Check options
    if( options.mol_inter and options.mol_intra and rank == 0  ):
	print " Options --mol_inter and --mol_intra are mutually exclusive "
	sys.exit("Error is specified options ")
    #
    #  Initialize blank system 
    #
    struc_o = StructureContainer(verbose=False)
    param_o = ParameterContainer()
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

    debug_chains = False  
    if( debug_chains ):
        for pid, pt_i  in struc_o.ptclC:
            print pid,pt_i.tagsDict["chain"]
        sys.exit(" chain debug ")
            
    # Print system properties
    if( rank == 0 ):
        print struc_o
        log_out.write(str(struc_o))

    p.barrier()

    search_i = dict()
    search_i = create_search(search_i,options.symbol_i,options.label_i,options.fftype_i,options.residue_i,options.resname_i,options.chains_i,options.rings_i)
    if( rank == 0 ):
        log_line = "\n Searching group i {} ".format(search_i)
        log_out.write(log_line)
        if( options.verbose ): print log_line
    list_i = struc_o.ptclC.getParticlesWithTags(search_i)
    N_i = len(list_i)
    if( rank == 0 ):
        log_line = "\n  group i has {} particles ".format(N_i)
        log_out.write(log_line)
        if( options.verbose ): print log_line

    #
    # If multi-core split the number of particles in group $i$ in the molecule onto each core
    #
    list_i_p  = p.splitListOnProcs(list_i)
    N_i_p = len(list_i_p)
    if(options.verbose ):
        if( len( list_i_p) > 0 ):
            print " cpu {} has {} particle numbers {} - {}  ".format(rank,len(list_i_p),list_i_p[0],list_i_p[len(list_i_p)-1])

    search_j = dict()
    search_j = create_search(search_j,options.symbol_j,options.label_j,options.fftype_j,options.residue_j,options.resname_j,options.chains_j,options.rings_j)
    if( rank == 0 ):
        log_line = "\n Searching group j {} ".format(search_j)
        log_out.write(log_line)
        if( options.verbose ): print log_line
    list_j = struc_o.ptclC.getParticlesWithTags(search_j)
    N_j = len(list_j)
    if( rank == 0 ):
        log_line = "\n  group j has {} particles ".format(N_j)
        log_out.write(log_line)
        if( options.verbose ): print log_line


    probabilityperpair = 1.0 
    consider_ilessj = False
    
    #if( list_i == list_j):
    #    consider_ilessj = True
    #    probabilityperpair = 2 
    #    if( options.verbose and rank == 0 ):
    #        print " Groups i and j are the same so i <= j will only be considered "
    
    # Create list for groups i and j
    #   the coordinate list of uni_i_p and uni_j
    #   is passed back from md analysis with indexes 0-N_i
    #   so:
    mda_list_i_p = range(N_i_p)
    mda_list_j = range(N_j)
    # Create include matrix of pairs to include
    # include_ij = [[None]*n_j]*n_i
    # had to use numpy list since regular list was asigning values include_ij[:][p_j] instead of include_ij[p_i][p_j]
    # not sure why
    if( N_i_p > 0 ): include_ij =  numpy.zeros((N_i_p,N_j), dtype=numpy.int)    
    # Allocate distance matrix 
    if( N_i_p > 0 ):  dist = numpy.zeros((N_i_p,N_j), dtype=numpy.float64)
    #
    # Find appropriate pairs
    #
    if( rank == 0 ):
        log_line = "\n  Finding {} x {}  pairs ".format(N_i_p,N_j)
        log_out.write(log_line)
        if( options.verbose ): print log_line
        
    debug_norm = False  
    debug_add = False  
    totalpair_cnt_p = 0
    pair_cnt_p = 0
    Nj_i_list_p = np.zeros( N_i , dtype=numpy.int)    
    Nj_i_list = np.zeros( N_i , dtype=numpy.int)    
    for a_i in  mda_list_i_p:
        pid_i =  list_i_p[a_i]
        Nj_sum = 0
        Nj_i = 0.0 
        for a_j in  mda_list_j:
            pid_j =  list_j[a_j]
            add_ij = 1
            if( pid_i == pid_j ):
                add_ij = 0
            if( pid_i != pid_j ):
                Nj_sum += 1 
            if( options.mol_inter and struc_o.ptclC[pid_i].tagsDict["chain"] == struc_o.ptclC[pid_j].tagsDict["chain"] ):
                add_ij = 0
            elif( options.mol_intra and struc_o.ptclC[pid_i].tagsDict["chain"] != struc_o.ptclC[pid_j].tagsDict["chain"] ):
                add_ij = 0
            #elif( consider_ilessj and atom_j.number <= atom_i.number ):
            #    add_ij = 0                
            #elif( options.exbonded and  struc_o.bondC.ij_bonded(atom_i.number+1,atom_j.number+1) ):
            #    add_ij = 0
            #    if( debug_add):
            #        print " found bonded ",atom_i.number+1,atom_j.number+1,struc_o.ptclC[atom_i.number+1].tagsDict["label"],struc_o.ptclC[atom_j.number+1].tagsDict["label"]
            if( debug_add ):
                print " debug_add ", pid_i , pid_j, struc_o.ptclC[pid_i].tagsDict["chain"] ,struc_o.ptclC[pid_j].tagsDict["chain"],add_ij
            include_ij[a_i][a_j] = add_ij
            if( add_ij == 1 ):
                Nj_i += probabilityperpair
                pair_cnt_p +=1
        if( options.truedensity ):
            Nj_i_list_p[a_i] = Nj_i
        else:            
            Nj_i_list_p[a_i] = Nj_sum
            totalpair_cnt_p += Nj_sum
        if( debug_norm ):
            print " atom ",pid_i, " has ",Nj_i,"allowed nieghbors and ",Nj_sum," possible neighbors "
            sys.exit("debug_norm")
            
    if( debug_norm ):
        sys.exit("debug_norm")


    p.barrier()
    pair_cnt = p.allReduceSum(pair_cnt_p)
    totalpair_cnt = p.allReduceSum(totalpair_cnt_p)
    for a_i in range(N_i):
        Nj_i_list[a_i] = p.allReduceSum(Nj_i_list_p[a_i])

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
        
    # Relabel segments to correspond to lists i on processor p
    for pid_i in list_i_p:
        if( debug ):
            print " Adding from list i",pid_i,struc_o.ptclC[pid_i].tagsDict["label"],struc_o.ptclC[pid_i].position," from u ",universe.atoms[pid_i-1].name , universe.coord[pid_i-1]
        universe.atoms[pid_i-1].resname = "grpi"
    uni_i_p = universe.selectAtoms(" resname grpi ")
    p.barrier()
    
    # Relabel segments to correspond to lists j
    for pid_j in list_j:
        if( debug ):
            print " Adding from list j",pid_j,struc_o.ptclC[pid_j].tagsDict["label"]," from u ",universe.atoms[pid_j-1].name
        universe.atoms[pid_j-1].resname = "grpj"
    uni_j = universe.selectAtoms(" resname grpj ")

    finial_frame = options.frame_f
    if( options.readall_f ):
        finial_frame = len(universe.trajectory)
            
    if( rank == 0 ):
        # log_line += "\n  %d "%()d
        log_line = sperator_line
        log_line += "\n  Date %s "%(str(datetime.datetime.now()))
        log_line += "\n  Number of processors  %d "%(size)
        log_line += "\n  Particles of group i: %d  "%(N_i)
        log_line += "\n  Particles of group j: %d  "%(N_j)
        log_line += "\n  Total pairs found %d "%(totalpair_cnt)
        log_line += "\n  Allowed pairs found %d "%(pair_cnt)
        if( options.truedensity ):
            log_line += "\n  The  Allowed pairs will be used to calculate the density "
        else:
            log_line += "\n  The  Total pairs will be used to calculate the density "
            
        log_line += "\n  probabilityperpair %f "%(probabilityperpair)
        log_line += "\n  <Nj(i)>  %f "%(np.average(Nj_i_list) )
        log_line += sperator_line
        log_line += "\n Frames "
        log_line += "\n     Initial frame  %d "%(options.frame_o)
        log_line += "\n     Step frame  %d  "%(options.frame_step)
        log_line += "\n     Final frame  %d  "%(finial_frame)
        log_out.write(log_line)
        if( options.verbose ):
            print log_line


    # Calculate rdf relate values
    sq_r_cut = options.r_cut**2
    n_bins = int(options.r_cut/options.bin_size) + 1 
    rdf_cnt_p = np.zeros(n_bins)    
    rdf_nn_cnt_p = np.zeros(n_bins)    # Nearest neighbor count 
    rdf_frames = 0
    volume_sum = 0.0 
    for ts in universe.trajectory:
        if( options.frame_o <= ts.frame ):
            if( ts.frame <= finial_frame  ):
                if( ts.frame%options.frame_step == 0 ):
                    rdf_frames += 1 
                    if( rank == 0 and options.verbose ):
                        print "Calculation %d frame %d/%d " % (rdf_frames,ts.frame, finial_frame)
                    volume_sum += ts.volume      # correct unitcell volume
                    box = ts.dimensions
                    coor_i = uni_i_p.coordinates()
                    coor_j = uni_j.coordinates()
                    distance_array(coor_i,coor_j, box, result=dist)  
                    for a_i in mda_list_i_p:
                        Nj_i = Nj_i_list_p[a_i]
                        a_i_hasnieghbor = False
                        r_ij_nn = options.r_cut   # Nearest Neighbor distance  
                        for a_j in mda_list_j:
                            # if( debug ):print " r_ij ",a_i,a_j,dist[a_i,a_j],include_ij[a_i][a_j]
                            if( include_ij[a_i][a_j] == 1 and dist[a_i,a_j] <= options.r_cut):
                                bin_index = int( round( dist[a_i,a_j] / options.bin_size) )
                                rdf_cnt_p[bin_index] += probabilityperpair/Nj_i
                                if( debug ):
                                    print " r_ij ",a_i,a_j,dist[a_i,a_j],include_ij[a_i][a_j]
                                    print rank," adding ij at index ",bin_index, probabilityperpair/Nj_i
                                a_i_hasnieghbor = True
                                if( dist[a_i,a_j] < r_ij_nn ): r_ij_nn = dist[a_i,a_j]                            
                        if( a_i_hasnieghbor ):
                            bin_nn_index = int( round( r_ij_nn / options.bin_size) )
                            rdf_nn_cnt_p[bin_nn_index] += probabilityperpair



   
    p.barrier()			
    if( options.verbose and rank == 0 ):
	print "      Finding averages "
    #
    # Find averages
    #
    rdf_cnt = np.zeros(n_bins)   
    rdf_nn_cnt = np.zeros(n_bins)   
    for bin_index in range( n_bins):
	# Sum rdf_cnt of each bin on each processor 
	rdf_cnt[bin_index] = p.allReduceSum(rdf_cnt_p[bin_index])
	rdf_nn_cnt[bin_index] = p.allReduceSum(rdf_nn_cnt_p[bin_index])
    p.barrier() # Barrier for MPI_COMM_WORLD
    
    #
    # Calculate rdf results 
    #
    '''

    g_ij(r) = n_j(r_ij)/ (rho_j 4 pi r^2 dr )
    g_ij(r) = n_j(r_ij)/ (rho_j 4/3 pi( r_out^3 - r_in^3)

    rho_j = N_j / V_ave

    g(r) = n(r) /( rho dV )

    n(r) = 1/Ni sum_i^{N_i} sum_j^{N_j} \gamma( r - r_{ij})

    g(r) = 1/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) / rho_j(i) ]

    Regular density
    g(r) =  1/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) / (  sum_j^{N_j} \gamma(  pair ij )/<V> )  ]
    g(r) =  <V>/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) /   sum_j^{N_j} \gamma(  pair ij ) ]
    
    True density 
    g(r) =  1/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) / (  sum_j^{N_j} \gamma( allowed pair ij )/<V> )  ]
    g(r) =  <V>/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) /   sum_j^{N_j} \gamma( allowed pair ij ) ]


    Nj_i = sum_j^{N_j} \gamma(  pair ij )
    rdf_cnt_p = sum_f^{N_frames}   sum_j^{N_j} \gamma( r - r_{ij})
    sum_j^{N_j} \gamma( r - r_{ij})  = rdf_cnt_p/N_frames

    '''
    total_cnts = np.sum( rdf_cnt )
    total_nn_cnts = np.sum( rdf_nn_cnt )
    box_vol_ave = volume_sum/float(rdf_frames)
    #  rho_box_j = 
    rho_box_j = np.average(Nj_i_list)  / box_vol_ave
    vol_cut = 4.0*math.pi/3.0*options.r_cut**3
    #  Calculate density based on spheres
    #    This is handy if you have vacuum or different phases in your system
    #    However, by no including some interaction the density can be off
    #       such ass mol_inter,mol_intra, exbonded
    n_spheres = float(N_i)*float(rdf_frames)
    vol_sheres  = vol_cut*n_spheres # /float(probabilityperpair) 
    rho_sphere_j =  np.average(Nj_i_list)  / vol_sheres 
    cnt_sum_j = 0.0 
    nn_cnt_sum_j = 0.0
    r_ij = []
    g_r_box = [] #.append( cnt_r_frame/dr_vol/rho_box_j )
    g_r_sphere = [] #.append( cnt_r_frame/dr_vol/sphere_den_j )
    g_r_nn_box = [] #.append( nn_cnt_r_frame/dr_vol/rho_box_j )
    nb_r = [] #.append(nb_cnt)
    nb_sum = [] #.append( cnt_sum_j )
    nn_nb_r = [] #.append(nb_cnt)
    nn_nb_sum = [] #.append( nn_cnt_sum_j )

    for bin_index in range(n_bins):
        
        r_val = options.bin_size*float(bin_index)
        dr_sq = r_val*r_val
        r_in = r_val - options.bin_size*0.5
        r_out = r_val + options.bin_size*0.5
        cnt_r_frame = float( rdf_cnt[bin_index] ) /float(rdf_frames)
        nn_cnt_r_frame = float( rdf_nn_cnt[bin_index] ) /float(rdf_frames)
        dr_vol = 4.0*math.pi/3.0*( r_out**3 - r_in**3 )
        dr_vol_apx = 4.0*math.pi*(  r_val**2 )*options.bin_size

        # n(r)  = 1/N_i  sum_j^{N_j} \gamma( r - r_{ij}) 
        nb_cnt = cnt_r_frame/float( N_i )
        cnt_sum_j += nb_cnt

        nn_nb_cnt = nn_cnt_r_frame/float( N_i )
        nn_cnt_sum_j += nn_nb_cnt
        
        r_ij.append(r_val)
        # g(r) = <V> * n(r) / dV 
        g_r_box.append( box_vol_ave*nb_cnt/dr_vol )
        g_r_nn_box.append( box_vol_ave*nn_nb_cnt/dr_vol )
        g_r_sphere.append( vol_sheres*nb_cnt/dr_vol )
        nb_r.append(nb_cnt)
        nb_sum.append( cnt_sum_j )
        nn_nb_r.append(nn_nb_cnt)
        nn_nb_sum.append( nn_cnt_sum_j )

        if( debug):
            print "n_bins,cnt_r_frame,dr_vol,rho_box_j",bin_index,nb_cnt,cnt_sum_j,nn_nb_cnt,nn_cnt_sum_j,dr_vol,rho_box_j

    if( rank == 0 ):
        
        # Write output
        #
        dat_lines ="#   Input "
        dat_lines +="\n#     Searching group i {} ".format(search_i)
        dat_lines +="\n#      N_i %d " % ( N_i )
        dat_lines +="\n#     Searching group j {} ".format(search_j)
        dat_lines +="\n#      N_j %d " % (N_j )
        dat_lines +="\n#    Frames: "
        dat_lines +="\n#      Initial  %d " %  (options.frame_o)
        dat_lines +="\n#      Step %d " %  (options.frame_step)
        dat_lines +="\n#      Final  %d  " %  (finial_frame)
        dat_lines +="\n#      Count  %d  " %  (rdf_frames)
        dat_lines +="\n#    Bin-size %f  " % (options.bin_size)
        dat_lines +="\n#    Cut-off %f  " % (options.r_cut)
        dat_lines +="\n#    Pairs found %d "%(pair_cnt)
        dat_lines +="\n#    probabilityperpair %f "%(probabilityperpair)
        if( options.truedensity ):
            dat_lines +="\n#   The true density is used in that only allowed pairs ij are used for the number density of group j "
            dat_lines +="\n#    <Nj(i)>  %f "%(np.average(Nj_i_list) )
            dat_lines +="\n#    Total_cnts %d  " % (total_cnts)
            dat_lines +="\n#    Total nearest neighbor cnts %d  " % (total_nn_cnts)
            dat_lines +="\n#    Box averages "
            dat_lines +="\n#      Average Box Volume %f " % ( box_vol_ave)
        dat_lines +="\n#      Box density j %f N A^-3 " % (rho_box_j )
        #dat_lines +="\n#    Sphere averages "
        #dat_lines +="\n#      N spheres %d " % (n_spheres )
        #dat_lines +="\n#      Sphere volume  %f A^3 " % (vol_cut )
        #dat_lines +="\n#      Average Sphere density j  %f N A^3 " % (rho_sphere_j )
        dat_lines +="\n#    "
        dat_lines +="\n#    r (Angstroms) "
        dat_lines +="\n#    g_r_box ()  "
        dat_lines +="\n#    g_r_nn_box () "
        dat_lines +="\n#    nb_r (pair_ij count/frame/N_i) "
        dat_lines +="\n#    nb_sum ( Sum(r<r_ij) pair_ij count/frame/N_i) "
        dat_lines +="\n#    nn_nb_r ( nearest neighbor pair_ij count/frame/N_i) "
        dat_lines +="\n#    nn_nb_sum ( Sum(r<r_ij) nearest neighbor pair_ij count/frame/N_i) "
        dat_lines +="\n# r ; g_r_box, g_r_nn_box, nb_r ,nb_sum,nn_nb_r,nn_nb_sum "
        dat_lines +="\n#   1        2         3         4      5      6       7        "
        for b_i in range(n_bins):
            r_val = options.bin_size*float(bin_index)
            dat_lines += "\n  {} {} {} {} {} {} {} ".format(r_ij[b_i],g_r_box[b_i], g_r_nn_box[b_i],nb_r[b_i],nb_sum[b_i],nn_nb_r[b_i],nn_nb_sum[b_i])
        if( options.verbose ):
            print dat_lines
        # Write data file 
        dat_file = options.output_id + ".dat"
        dat_out = open(dat_file,"w") 
        dat_out.write(dat_lines)
        dat_out.close()            
            
        log_line="\n  Finished   "
        print log_line
        log_out.write(log_line)
	log_out.close()
	

if __name__=="__main__":
    main(debug=False)
   
