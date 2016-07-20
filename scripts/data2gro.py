#! /usr/bin/env python

# Dependency modules 
import numpy as np 
import sys
import copy

# Streamm toolkit modules 
from structureContainer import StructureContainer
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
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input .data LAMMPS file ")
    parser.add_option("--out_id", dest="out_id", type="string", default="out", help="Output id")
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
    

def add_gtype_prop(strucC):
    """
     Add in residue
     """

    for pid, pt_i  in strucC.ptclC:
        add_dict = pt_i.tagsDict
        add_dict["gtype"] = pt_i.tagsDict["symbol"]
        pt_i.setTagsDict(add_dict)                    
        
def add_residue_prop(strucC):
    """
     Add in residue
     """

    for pid, pt_i  in strucC.ptclC:
        add_dict = pt_i.tagsDict
        add_dict["residue"] = pt_i.tagsDict["chain"]
        add_dict["qgroup"] = pt_i.tagsDict["chain"]
        add_dict["resname"] = "MOLRES"
        pt_i.setTagsDict(add_dict)                    
        
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

    struc_o,param_o = read_lmpdata( struc_o , param_o , options.in_data)
    print struc_o
    print param_o
    add_residue_prop(struc_o)
    add_gtype_prop(struc_o)

    print " Writing xyz file {}.xyz".format(options.out_id)
    comment = " From data file {} ".format( options.in_data)
    append = False 
    struc_o.ptclC.write_xmol("{}.xyz".format(options.out_id),comment,append)
    print " Writing gromacs files {} ".format(options.out_id)
    write_gro(struc_o,"{}.gro".format(options.out_id))
    write_top(struc_o,"{}.top".format(options.out_id),"{}.itp".format(options.out_id))
    write_itp(param_o,"{}.itp".format(options.out_id))

    
    
if __name__=="__main__":
    main()
