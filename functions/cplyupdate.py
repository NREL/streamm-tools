#! /usr/bin/env python

# Dependency modules 
import numpy as np 
import sys
import copy

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

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("--in_cply", dest="in_cply", type="string", default="", help="Input .cply reference file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input .data LAMMPS file ")
    parser.add_option("--mol_pbcs", dest="mol_pbcs",default=False,action="store_true", help="Apply pbcs to creat whole molecules")    
    parser.add_option("-o","--out_id", dest="out_id", type="string", default="out", help="Output id")
    #
    (options, args) = parser.parse_args()
        
    return options, args


def read_lmpdata( strucC , parmC , data_file):
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

            #tagsD = {"chain":chain_i,"symbol":el.symbol,"number":el.number,"mass":el.mass,"cov_radii":el.cov_radii,"vdw_radii":el.vdw_radii}
            #if( pt_overwrite ):
            pt_i = strucC.ptclC[ind+1]
            pt_i.position = r_i

            update_coord = True
            if( not update_coord ):

                pt_i.type = el.symbol
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
                add_dict["fftype"] = fftype_i 
                # 
                add_dict["ffmass"] = ATYPE_MASS[indx]
                add_dict["residue"] = chain_i
                add_dict["qgroup"] = chain_i
                add_dict["resname"] = "MOLR"
                add_dict["ring"] = 0
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
    Read in gaussian fchk file and create an xyz file 
    """

    #        
    # Read options 
    #
    options, args = get_options()
    #
    #  Initialize blank system 
    #

    bb_o = Buildingblock()
    bb_o.read_cply(options.in_cply)

    param_o = ParameterContainer()
    bb_o,param_o = read_lmpdata( bb_o , param_o , options.in_data)
    print bb_o
    print param_o

    if( options.mol_pbcs  ):
        bb_o.bondC_nblist()
        molpbcs(bb_o,bb_o.bonded_nblist, bb_o.bonded_nbindx,debug=False)

    comment = " From data file {} ".format( options.in_data)
    append = False
    out_xyz = "{}.xyz".format(options.out_id)
    bb_o.ptclC.write_xmol("{}".format(out_xyz),comment,append)

    if( options.verbose ): print " Writing cply file {}.cply".format(options.out_id)
    bb_o.write_cply("{}.cply".format(options.out_id),write_ff=True,write_bonds=True)
    
if __name__=="__main__":
    main()
