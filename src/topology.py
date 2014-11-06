"""
Process topology information
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov

# Conversions
KJ_KCAL = 0.23901
GRO_SIG = 5.61230943
NM_ANG = 10.0

import sys
import numpy as np 
import datetime
import copy

from bonds         import Bond,     BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer

from parameters import ParameterContainer
#from particles import ParticleContainer
from structureContainer import StructureContainer
from parameters import imptype

#from periodictable import periodictabled

import atomtypes, pbcs, xmol 

def create_top(strucC,ff_charges): # Move out of class (or derived class)
    """
    Find topology information for force-field input files
    """



    # New options that need to be passed
    rerun_bonded = False  
    set_biaryl = True 
    set_ptma = False

    limdih = 0
    limitdih_n = 1
    
    verbose = True
    debug = True



    if( debug):
        print " Read in %d Bonds "%len(strucC.bondC)
        print " Read in %d Angles "%len(strucC.angleC)
        print " Read in %d Dihedrals "%len(strucC.dihC)

    # Check bonds
    if( len(strucC.bondC) <= 0 or rerun_bonded ):
        # Create nearest neighbor list
        cov_nblist, cov_nbindx = build_covnablist(strucC)
        if( rerun_bonded ):
            strucC.bondC.clear()
        nblist_bonds(strucC,cov_nblist, cov_nbindx)
        if( verbose ):
            print " Found bonds %d "%(len(strucC.bondC))
    else:
        # Build nieghbor list with provided bonds 
        cov_nblist, cov_nbindx = bonded_nblist(strucC) 
        if( verbose ):
            print " Read in bonds %d "%(len(strucC.bondC))
    # Check angles
    if( len(strucC.angleC) <= 0  or rerun_bonded ):
        if( rerun_bonded ):
            strucC.angleC.clear()
        nblist_angles(strucC,cov_nblist, cov_nbindx)
        if( verbose ):
            print " Found angles %d "%(len(strucC.angleC))
    # Check dihedrals
    if( len(strucC.dihC) <= 0  or rerun_bonded ):
        if( rerun_bonded ):
            strucC.dihC.clear()
        nblist_dih(strucC,cov_nblist, cov_nbindx)
    if( verbose ):
        print " Found dihedrals %d "%(len(strucC.dihC))

    # Find rings
    strucC , ring_nblist, ring_nbindex = find_rings(strucC,cov_nblist, cov_nbindx)

    debug = False 
    if( debug):
        for pid_o, ptclObj_o  in strucC.ptclC:
            print "create_top particles ",pid_o,ptclObj_o.tagsDict["symbol"],ptclObj_o.tagsDict["resname"],ptclObj_o.tagsDict["ring"]
            # Write xyz file
            ptclObj_o.tagsDict['symbol'] = ptclObj_o.tagsDict["ring"]
            comment = "ring check"
            append = False
            xmol_file = "rings.xyz"
            xmol.write(strucC.ptclC, xmol_file,comment,append)

        sys.exit("create_top 2 print particles in strucC.ptclC")


    # Asign atom types
    strucC = atomtypes.oplsaa( ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx )
    
    if(set_biaryl):
        strucC = atomtypes.biaryl_types( ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx )
        strucC = atomtypes.interring_types( ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx )

    if(set_ptma):
        strucC = atomtypes.set_pmmatypes( ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx )
        strucC = atomtypes.set_ptmatypes( ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx )

    #if( set_qgroup ):
    strucC = set_chargegroups( ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx )
    qgroup_maxvol(strucC)
    
    if(debug):
        print len(strucC.ptclC)
        for pid_i in range(1,len(strucC.ptclC)+1):
            print pid_i, strucC.ptclC[pid_i].tagsDict["number"],strucC.ptclC[pid_i].tagsDict["ring"],strucC.ptclC[pid_i].tagsDict["fftype"]

    return strucC,cov_nblist, cov_nbindx


def set_param(struc_o,param_all,norm_dihparam,cov_nblist, cov_nbindx):
    """
    Set force-field parameters
    """

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
    #
    paramC_p.set_nbfunc( param_all.get_nbfunc() )
    paramC_p.set_combmixrule( param_all.get_combmixrule() )
    paramC_p.set_genpairs( param_all.get_genpairs() )
    paramC_p.set_fudgeLJ( param_all.get_fudgeLJ() )
    paramC_p.set_fudgeQQ( param_all.get_fudgeQQ() )
    #
    # Count atom types
    #
    debug = 0
    for pid_o, ptclObj_o  in ptclC_o:    
        new_type = True
        lj_p = 0
        fftype_i = ptclObj_o.tagsDict["fftype"]
        #mass_i = ptclObj_o.mass
        for lj_p, ljObj_p  in ljtypC_p:    
            if( fftype_i == ljObj_p.ptype1 ):
                new_type = False
                ptclObj_o.type = str(lj_p)
                ptclObj_o.tagsDict["lmptype"] = lj_p
                if(debug): print ' maches previous ', lj_p,ljObj_p.ptype1
                
        if( new_type ):
            # Find type in ljtypC_all
            cnt_check = 0
            type_found = False 
            ptclObj_o.type = str(lj_p+1)
            ptclObj_o.tagsDict["lmptype"] = lj_p+1
            for lj_all, ljObj_all  in ljtypC_all:
                if( fftype_i == ljObj_all.ptype1):
                    cnt_check += 1
                    type_found = True 
                if( type_found ):
                    # Set mass of particle type
                    #   This is need for some force-field input files (LAMMPS)
                    #mass_i = 
                    ljObj_all.setpid( ptclObj_o.tagsDict["number"] )
                    #ljObj_all.setmass( mass_i )
                    
                    
                    ljtypC_p.put(ljObj_all)
                    type_found = False 
                    
                    
            if( cnt_check < 1 ):
                print " No LJ parameters were found for atom type %s ",fftype_i
                #raise TypeErrorb
            elif( cnt_check > 1 ):
                print " Multiple LJ parameters were found for atom type %s ",fftype_i
                #raise TypeError
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
        bond_len = np.linalg.norm(pbcs.delta_r_c(r_i,r_j,struc_o.get_latvec() ))
        
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
        r_ik = pbcs.delta_r_c(r_i,r_k,struc_o.get_latvec() )
        r_ij = pbcs.delta_r_c(r_i,r_j,struc_o.get_latvec() )
        angle_kij = pbcs.getAngle(r_ik,r_ij)
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
                    print " dihen_norm ",dihen_norm

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
                print  j
                #    el_j = ELN[j]
                #    ELCNT[j] = ELCNT[j] + 1

        sys.exit('debug build_covnablist')

    return (cov_nblist, cov_nbindx)


def build_nablist_bonds(strucC):
    """
    Build nieghbor list using bonded information 
    """
    import sys,numpy

    debug = 0
    NNAB  = 0
    
    maxnnab = len(strucC.bondC)*2 + 1
    #NBLIST = numpy.empty( maxnnab,  dtype=int )
    #NBINDEX = numpy.empty( maxnnab,  dtype=int )
    
    # python style nieghbor list
    nblist_py = [] #numpy.empty( maxnnab,  dtype=int )
    
    NBLIST = []
    NBINDEX = []
    NBLIST.append( 0 )
    
    for atom_i in range( len(strucC.ptclC)):
	nblist_py.append( [  ] )
        
    # Loop over bonds and record neighbors 
    for b in range(  len(BONDS) ) :
	bnd_i = BONDS[b][0]
	bnd_j = BONDS[b][1]
	
	nblist_py[bnd_i].append( bnd_j )
	nblist_py[bnd_j].append( bnd_i )
	
    # Translate 2D into 1D array
    #   mostly to match perviously writen fortran code
    
    for nbs_i in range( len(nblist_py)):
	nlist_i = nblist_py[nbs_i]
	NBINDEX.append( NNAB + 1 )
	for nb_i in  nlist_i:
	    NNAB +=  1
	    NBLIST.append( nb_i )

    NBINDEX.append( NNAB + 1 )
    
	
    debug = 0
    if ( debug ):
       print ' total nbs ',NNAB,' total bonds ',len(BONDS)
       for i in range( len(ELN)):
	    N_o = NBINDEX[ i  ]
	    N_f = NBINDEX[ i + 1 ] - 1
	    NNAB = N_f - N_o + 1
	    print ' atom ', i,' ',ELN[i],NNAB,N_o,nblist_py[i]
	    
	#sys.exit('debug')
    return (NBLIST,NBINDEX)



def bonded_nblist(strucC):
    """
    Create neighbor list of bonded particles
    """
        
    debug = False

    NNAB  = 0

    maxnnab = len(strucC.bondC)*2 + 1

    if(debug):
        print " maxnnab",maxnnab

    #NBLIST = numpy.empty( maxnnab,  dtype=int )
    #NBINDEX = numpy.empty( maxnnab,  dtype=int )

    # python style nieghbor list
    nblist_py = [] #numpy.empty( maxnnab,  dtype=int )

    NBLIST = []
    NBINDEX = []
    NBLIST.append( 0 )

    # First create an N diminsional list of index lists for each particle

    for p_i, prtclC in strucC.ptclC:
        nblist_py.append( [  ] )
    nblist_py.append( [  ] )

    # bassed on bonds add index of neighbros to particle index of nblist_py
    for b_i,bondObj in  strucC.bondC:
        bnd_i = bondObj.pgid1 
        bnd_j = bondObj.pgid2 

        nblist_py[bnd_i].append( bnd_j )
        nblist_py[bnd_j].append( bnd_i )

        if(debug): print " adding bond bonded_nblist",bnd_i,bnd_j

    # Translate 2D into 1D array
    #   mostly to match perviously writen fortran code
    for nblis_indx in range( len(nblist_py)):
        # loop over  each particle p_i and get list of neighbors nlist_i
        nlist_i = nblist_py[nblis_indx]
        NBINDEX.append( NNAB + 1 )
        # Loop over neighbor list of each particle nlist_i and get neighbor p_j
        p_i = nblis_indx #+ 1
        for p_j in  nlist_i:

            if(debug):
                print "nblis_indx",nblis_indx
                print " p_i ",p_i," p_j ",p_j

            #if( p_j > p_i):
            # remove redundent neighbors 
            NNAB +=  1
            # add to neighbor list 
            NBLIST.append( p_j )

    NBINDEX.append( NNAB + 1 )

    if ( debug ):
        print ' total nbs ',NNAB

        for p_i, prtclC in strucC.ptclC:
            N_i_o = NBINDEX[p_i]
            N_i_f = NBINDEX[p_i+1]
            print " atom ",p_i,prtclC.type, " has ",N_i_f - N_i_o

            for indx_j in range( N_i_o,N_i_f):
                atom_j = NBLIST[indx_j]
                print "      nb ",atom_j, strucC.ptclC[atom_j].type

        sys.exit('bonded_nblist debug')

    return (NBLIST,NBINDEX)


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


def id_ring3(pid_o,ptclObj_o,strucC,cov_nblist, cov_nbindx):
    import sys
    """
    Find atoms in conjugated rings 
    """
    
    max_paths = 1000
    
    debug = False
    
    R_SET = []
    RING_ATOMS = []
    BAD_PATH = []

    #a_i = pid_o -1 
    one = 1
    zero = 0
    atoms_in_ring = 0
    # relabel based on neighbors
    n_ptcl = len(strucC.ptclC)

    if(debug):
        print " for a system of %d partiles checking particle %d "%(n_ptcl,pid_o)
        
    for pid_i in range(n_ptcl+1): 
        R_SET.append(one)
        
    R_SET[pid_o] = 0 
    r_term = 1
    cnt = 0
    p_cnt = 0
    if(debug): print ' initializing ring ',pid_o," w NN ", calc_nnab(pid_o,cov_nbindx)

    last_i = pid_o
    NNAB_last = calc_nnab(pid_i,cov_nbindx)
    while ( r_term ):
        N_o = cov_nbindx[last_i]
        N_f = cov_nbindx[last_i+1] - 1

        if( debug): print " checking neighbors ",N_o,N_f

        for n_indx in range( N_o,N_f+1):
            j = cov_nblist[n_indx]
            cnt = cnt + 1
            # If non of the neighboring atoms of atom j are conjugated
            #    make the path as bad 
            if( cnt > NNAB_last+1 ):
                p_cnt += 1
                if(debug): print ' bad path found resetting at atom ',j,strucC.ptclC[j].tagsDict["number"]
                for ring_a in range(len(RING_ATOMS)):
                    j = RING_ATOMS[ring_a]
                    
                BAD_PATH.append(j)
                # Reset lists and counts
                RING_ATOMS = []
                for pid_k in range(n_ptcl):
                    R_SET[pid_k] = 1
                atoms_in_ring = 0
                R_SET[pid_o] = 0 
                cnt = 0
                last_i = pid_o
                
                if(debug): print '  resetting last_i to ',pid_o,ptclObj_o.tagsDict["number"]
                    
                for bad_i in range( len(BAD_PATH)):
                    j = BAD_PATH[bad_i]
                    R_SET[j] = 0
                    if(debug): print '  bad path atoms ',strucC.ptclC[j].tagsDict["number"]
                break

            # If path traces back to original atom ring has been found 
            if( atoms_in_ring  > 1 and j == pid_o ):
                r_term = 0
                if(debug): print ' ring found with ',atoms_in_ring 
                break

            number_j = strucC.ptclC[j].tagsDict["number"]
            NNAB_j = calc_nnab(j,cov_nbindx)
            ELCNT_j = calc_elcnt(j,strucC,cov_nblist,cov_nbindx)
            ring_type = False 
	    
            if(debug):
		print '           with ', number_j
		print '           with ', NNAB_j,
		print '           with ', R_SET[j]

            # Test if atom j is conjugated 
            if ( number_j == 6 and NNAB_j == 3 and R_SET[j] == 1 ):                      ring_type = True
            if ( number_j == 16 and NNAB_j == 2 and R_SET[j] == 1 and ELCNT_j[6] == 2 ): ring_type = True
            if ( number_j == 7 and NNAB_j >= 2 and R_SET[j] == 1 ):                      ring_type = True
            if( ring_type ):
                atoms_in_ring = atoms_in_ring + 1
                RING_ATOMS.append(j)
                R_SET[j] = 0
                last_i = j
                cnt = 0
                NNAB_last = NNAB_j 
                if(debug): print '   atom ',j,number_j,' added '
                break
            
        if (p_cnt > max_paths ):
            if(debug): print '     max paths considered ',max_paths
            r_term = 0
            
    return RING_ATOMS


def id_ring(pid_o,ptclObj_o,strucC,cov_nblist, cov_nbindx):
    import sys
    """
    Find atoms in conjugated rings 
    """
    
    debug = True
    if( pid_o ==7 ): debug = True
    
    R_SET = []
    RING_ATOMS = []
    BAD_PATH = []

    #a_i = pid_o -1 
    one = 1
    zero = 0
    atoms_in_ring = 0
    # relabel based on neighbors
    n_ptcl = len(strucC.ptclC)

    if(debug):
        print " for a system of %d partiles checking particle %d "%(n_ptcl,pid_o)
        
    for pid_i in range(n_ptcl+1): 
        R_SET.append(one)
        
    R_SET[pid_o] = 0 
    r_term = 1
    cnt = 0
    p_cnt = 0
    if(debug): print ' initializing ring ',pid_o," w NN ", calc_nnab(pid_o,cov_nbindx)

    last_i = pid_o
    NNAB_last = calc_nnab(pid_i,cov_nbindx)
    while ( r_term ):
        N_o = cov_nbindx[last_i]
        N_f = cov_nbindx[last_i+1] - 1

        if( debug): print " checking neighbors ",N_o,N_f

        for n_indx in range( N_o,N_f+1):
            j = cov_nblist[n_indx]
            cnt = cnt + 1
            # If non of the neighboring atoms of atom j are conjugated
            #    make the path as bad 
            if( cnt > NNAB_last+1 ):
                p_cnt = p_cnt + 1
                if(debug): print ' bad path found resetting at atom ',j,strucC.ptclC[j].tagsDict["number"]
                for ring_a in range(len(RING_ATOMS)):
                    j = RING_ATOMS[ring_a]
                    
                BAD_PATH.append(j)
                # Reset lists and counts
                RING_ATOMS = []
                for pid_k in range(n_ptcl):
                    R_SET[pid_k] = 1
                atoms_in_ring = 0
                R_SET[pid_o] = 0 
                cnt = 0
                last_i = pid_o
                
                if(debug): print '  resetting last_i to ',pid_o,ptclObj_o.tagsDict["number"]
                    
                for bad_i in range( len(BAD_PATH)):
                    j = BAD_PATH[bad_i]
                    R_SET[j] = 0
                    if(debug): print '  bad path atoms ',strucC.ptclC[j].tagsDict["number"]
                break

            # If path traces back to original atom ring has been found 
            if( atoms_in_ring  > 1 and j == pid_o ):
                r_term = 0
                if(debug): print ' ring found with ',atoms_in_ring 
                break

            number_j = strucC.ptclC[j].tagsDict["number"]
            NNAB_j = calc_nnab(j,cov_nbindx)
            ELCNT_j = calc_elcnt(j,strucC,cov_nblist,cov_nbindx)
            ring_type = False 
	    
            if(debug):
		print '           with ', number_j
		print '           with ', NNAB_j,
		print '           with ', R_SET[j]

            # Test if atom j is conjugated 
            if ( number_j == 6 and NNAB_j == 3 and R_SET[j] == 1 ):                      ring_type = True
            if ( number_j == 16 and NNAB_j == 2 and R_SET[j] == 1 and ELCNT_j[6] == 2 ): ring_type = True
            if ( number_j == 7 and NNAB_j >= 2 and R_SET[j] == 1 ):                      ring_type = True
            if( ring_type ):
                atoms_in_ring = atoms_in_ring + 1
                RING_ATOMS.append(j)
                R_SET[j] = 0
                last_i = j
                cnt = 0
                NNAB_last = NNAB_j 
                if(debug): print '   atom ',j,number_j,' added '
                break
            
        if (p_cnt > 100 ):
            if(debug): print '     max paths considered ',100
            r_term = 0
            
    return RING_ATOMS

def ring_type(pid_i,ptclObj_o,strucC,cov_nblist, cov_nbindx):
    """
    Criteria for a conjugated particle 
    """
    debug = False
    
    ring_type = False 
    number_j = strucC.ptclC[pid_i].tagsDict["number"]
    NNAB_j = calc_nnab(pid_i,cov_nbindx)
    ELCNT_j = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)

    if(debug):
        print '            atomic # %d  NN  %d '%(number_j, NNAB_j)

    # Test if atom j is conjugated 
    if ( number_j == 6 and NNAB_j == 3 ):                      ring_type = True
    if ( number_j == 16 and NNAB_j == 2 and ELCNT_j[6] == 2 ): ring_type = True
    if ( number_j == 7 and NNAB_j >= 2  ):                      ring_type = True

    return ring_type


def id_ring2(pid_o,ptclObj_o,strucC,cov_nblist, cov_nbindx):
    import sys
    """
    Find atoms in conjugated rings 
    """

    max_ring_particles = 30
    max_paths = 100
    
    debug = True
    #if( pid_o ==7 ): debug = True
    
    RING_ATOMS = []
    BAD_PATH = []

    #a_i = pid_o -1 
    one = 1
    zero = 0
    atoms_in_ring = 0
    # relabel based on neighbors
    n_ptcl = len(strucC.ptclC)

    if(debug):
        print " for a system of %d partiles checking particle %d "%(n_ptcl,pid_o)
        
    for pid_i in range(n_ptcl+1): 
        BAD_PATH.append([])
        
    
    r_term = 1
    cnt = 0
    path_cnt = 0
    if(debug): print ' initializing ring with particle ',pid_o 

    pid_i = pid_o
    last_j = 0 
    NNAB_last = calc_nnab(pid_i,cov_nbindx)
    while ( r_term ):
        N_o = cov_nbindx[pid_i]
        N_f = cov_nbindx[pid_i+1] - 1

        if( debug): print " checking pid_i %d w %d neighbors "%(pid_i,calc_nnab(pid_i,cov_nbindx))
        #           j
        #        /
        #  pid_i - j
        #        \
        #           j
        #
        path_found = False 
        for n_indx in range( N_o,N_f+1):
            j = cov_nblist[n_indx]
            if(debug):
                print " checking NN %d "%(j)

            # If path traces back to original atom ring has been found 
            if( atoms_in_ring  > 1 and j == pid_o ):
                r_term = 0
                path_found = True
                if(debug): print ' ring found with ',atoms_in_ring 
                break
            j_good = True
            if( debug ): print " check if j is part of a bad path ",BAD_PATH[pid_i] 
            if( j in BAD_PATH[pid_i] ):
                if(debug): print " j marked as a bad path "
                j_good = False
            # Prevent looping back
            if( j in RING_ATOMS or j == pid_o ):
                if(debug): print " j already in path "
                j_good = False
            if(  j_good ):
                if( ring_type(j,ptclObj_o,strucC,cov_nblist, cov_nbindx) ):
                    atoms_in_ring += 1
                    RING_ATOMS.append(j)
                    last_j = pid_i  
                    pid_i = j           # set j to 
                    if(debug):
                        print "   atom %d added atom # %d "%(j,strucC.ptclC[j].tagsDict["number"])
                        print "      last_j ",last_j
                        print "      pid_i ",pid_i

                    path_found = True 
                    break
            else:
                if(debug): print " skipping previously found to prevent looping back ",last_j

        if( not path_found ):
            if(debug):
                print " path not found pid_i %d last_j %d "%(pid_i,last_j)
                
                
            if( j == pid_o ):
                if(debug): print ' path has been reduced back to original atom %d '%(j)
                RING_ATOMS = []
                r_term = 0
                
            BAD_PATH[last_j].append( pid_i )
            if( debug): print "setting bad path %d for last_j %d "%(pid_i,last_j)
            
            if( len(RING_ATOMS) > 1 ):
                if(debug): print " pre pop RING_ATOMS ",RING_ATOMS
                RING_ATOMS.pop( len(RING_ATOMS) -1 )
                if(debug): print " post pop RING_ATOMS ",RING_ATOMS

                pid_i = last_j
                last_j = RING_ATOMS[-1]
            else:
                pid_i = pid_o
                
            
            if( debug):
                print " Bad path found no conjugated atoms " # attached to  %d stepping back to %d "%(BAD_PATH[-1],pid_i)
                print "RING_ATOMS ",RING_ATOMS
                print "BAD_PATH ",BAD_PATH
            
        
        if ( len(RING_ATOMS) > max_ring_particles ):
            if(debug):
                print '   maximum possible particles in ring %d exceeded '%max_ring_particles
                print "   pid_i %d last_j %d "%(pid_i,last_j)
            # Set last atom as a bad path
            #BAD_PATH[pid_i]
            
            # Reset lists and counts
            RING_ATOMS = []
            atoms_in_ring = 0
            last_j = 0
            pid_i = pid_o

            path_cnt += 1
            if( path_cnt > max_paths ):
                r_term = 0
            
    return RING_ATOMS


def find_rings(strucC,cov_nblist, cov_nbindx):
    """
    Find conjugate rings
    """
    import sys

    debug = False
    
    RINGLIST = []
    RINGINDEX = []
    RING_NUMB = []
    
    IN_RING = []

    n_ptcl = len(strucC.ptclC)
    one = 0 
    zero = 0 
    ring_cnt = 0
    a_cnt = 0
    # relabel based on neighbors
    for pid_i, ptclObj_i  in strucC.ptclC:
        IN_RING.append(one)
        RINGLIST.append(zero)
        RING_NUMB.append(zero)
        RINGINDEX.append(zero)

    RING_NUMB.append(zero)
    RINGLIST.append(zero)
    RINGINDEX.append(zero)

        
    # relabel based on neighbors
    for pid_i, ptclObj_i  in strucC.ptclC:
        RING_ATOMS = []
        # If particle i is conjugated 
        if( ring_type(pid_i,ptclObj_i,strucC,cov_nblist, cov_nbindx) ):
            RING_ATOMS = id_ring3(pid_i,ptclObj_i,strucC,cov_nblist, cov_nbindx)
        if( debug ):
            print pid_i,RING_ATOMS
            
        if ( len(RING_ATOMS) > 1 ):
            # If member of ring alread part of another ring add new ring to existing ring
            r_numb = 0 
            for ring_a in range(len(RING_ATOMS)):
                j = RING_ATOMS[ring_a]
                if( RING_NUMB[j] != 0 ):
                    r_numb = RING_NUMB[j]
            if( r_numb == 0 ):
                ring_cnt = ring_cnt + 1
                r_numb = ring_cnt 
                #print ' new ring ',ring_cnt
                 
            for ring_a in range(len(RING_ATOMS)):
                j = RING_ATOMS[ring_a]
                if( RING_NUMB[j] == 0 ): 
                    a_cnt = a_cnt + 1
                    RING_NUMB[j] = r_numb

    a_cnt = 0
    for r_numb in range(1,ring_cnt+1):
        RINGINDEX[r_numb] = a_cnt + 1
        for i in range(n_ptcl):
            if( RING_NUMB[i] == r_numb ):
                a_cnt = a_cnt + 1
                RINGLIST[a_cnt] = i
                
    RINGINDEX[ring_cnt+1] = a_cnt + 1


    # Set attached H to have same ring #
    attach_H = False
    if( attach_H ):

        for r_numb in range(1,ring_cnt+1):
            N_o = RINGINDEX[r_numb]
            N_f = RINGINDEX[r_numb+1] - 1
            print ' Ring ', r_numb
            for r_indx in range(N_o,N_f+1):
                i = RINGLIST[r_indx]
                N_o = cov_nbindx[i]
                N_f = cov_nbindx[i+1] - 1
                for n_indx in range( N_o,N_f+1):
                    j = cov_nblist[n_indx]
                    if( strucC.ptclC[j].tagsDict["number"] == 1 ):
                        RING_NUMB[j] = r_numb
		

    # update particles 
    for pid_i, ptclObj_i  in strucC.ptclC:
        ptclObj_i.tagsDict["ring"] =  RING_NUMB[pid_i]
        print ' atom  ',pid_i,ptclObj_i.tagsDict["number"],ptclObj_i.tagsDict["ring"]
	    		
	    

        #debug = True
    if(debug):
        print ring_cnt , " rings found "
        for r_numb in range(1,ring_cnt+1):
            N_o = RINGINDEX[r_numb]
            N_f = RINGINDEX[r_numb+1] - 1
            print ' Ring ', r_numb
            for r_indx in range(N_o,N_f+1):
                i = RINGLIST[r_indx]
                print ' type ',strucC.ptclC[i].tagsDict["number"],' or ',i,RING_NUMB[i]
                #sys.exit('find_rings')
                
        for pid_i, ptclObj_i  in strucC.ptclC:    
            print ' atom  ',pid_i,ptclObj_i.tagsDict["number"],ptclObj_i.tagsDict["ring"]
	    
	sys.exit('top.find_rings')

    return ( strucC, RINGLIST, RINGINDEX  )

def get_nrings(strucC):
    """
    Return the number of rings in a structure
    """
    nrings = 0 
    for pid_i, ptclObj_i  in strucC.ptclC:
        if( ptclObj_i.tagsDict["ring"] > nrings ):
            nrings = ptclObj_i.tagsDict["ring"]
    return nrings

def set_chargegroups(ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
    """
    Set charge groups 
    """
    debug = False
    verbose = True 
    
    if( verbose ):
	print "   Setting charge groups "
    
    
    # initialize 
    n_groups = 0


    #CG_SET = []
    #for pid_i, ptclObj_i  in strucC.ptclC:
    #        CG_SET.append(True) 

    for pid_i, ptclObj_i  in strucC.ptclC:
        ptclObj_i.tagsDict["qgroup"] = 0
        
    # set rings as charge groups
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group 
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            if( ptclObj_i.tagsDict["ring"] > 0 ):
                ptclObj_i.tagsDict["qgroup"] = ptclObj_i.tagsDict["ring"]
                print " setting qgroup  ring ",pid_i,ptclObj_i.tagsDict["symbol"],ptclObj_i.tagsDict["qgroup"],ptclObj_i.tagsDict["ring"]
                if(  ptclObj_i.tagsDict["ring"] >  n_groups ):
                    n_groups = n_groups + 1
                
    #if( debug): sys.exit('ring charge groups')
            
    # set methyls and alkyls 
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            if( ptclObj_i.tagsDict["number"] == 6 and ELCNT[1] > 0 ):  # CH_n ( Methyls / alkyl ) 
                n_groups = n_groups + 1
                #CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["number"] == 1 ):
                        ptclObj_j.tagsDict["qgroup"] = n_groups
                        #CG_SET[j] = 0
			
            if(  ptclObj_i.tagsDict["fftype"]  == 'C' ):  # 
                n_groups = n_groups + 1
                #CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["number"] == 8 ):
                        ptclObj_j.tagsDict["qgroup"] = n_groups

                    #if ( ELN[j] == 8 ):
                    #    CHARN[j] = n_groups
                    #    CG_SET[j] = 0

            if(  ptclObj_i.tagsDict["fftype"] == 'ON' ):  # nitroxide
                n_groups = n_groups + 1
                #CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] == 0 ):
                        ptclObj_j.tagsDict["qgroup"] = n_groups

                    #j = cov_nblist[j_indx]
                    #if ( CG_SET[j]  ):
		    #CHARN[j] = n_groups
                    #	CG_SET[j] = 0

    # check unset heavy atom groups
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            if( ptclObj_i.tagsDict["number"] == 6 or ptclObj_i.tagsDict["number"] == 8 ):  #
		q_g = 0 
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] != 0 ):
                        q_g = ptclObj_j.tagsDict["qgroup"]
                        break

		if( q_g  == 0 ):
			
		    n_groups = n_groups + 1
		    #CG_SET[i] = 0
		    ptclObj_i.tagsDict["qgroup"] = n_groups
		else:
		    #CG_SET[i] = 0
		    ptclObj_i.tagsDict["qgroup"] = q_g
		    
		#print " unset carbon set to ",ptclObj_i.tagsDict["qgroup"]
			

            if( ptclObj_i.tagsDict["fftype"] == 'CA' ):  # conjugate
                n_groups = n_groups + 1
                # CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
 
                    

    # set light atom groups 
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            
            if ( ptclObj_i.tagsDict["number"] == 1 or  ptclObj_i.tagsDict["number"] == 9 ):
                j_set = 0
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] != 0 ):
                        n_groups = ptclObj_j.tagsDict["qgroup"]
                        j_set = 1
                if( j_set ):
                    ptclObj_i.tagsDict["qgroup"] = n_groups
                    #CG_SET[i] = 0
                else :
                    print ptclObj_i.tagsDict["symbol"] ,' connected to unset atom '
                    print ptclObj_i.tagsDict["symbol"] , " - ",ptclObj_j.tagsDict["symbol"] 
                    print ptclObj_i.tagsDict["fftype"] , " - ",ptclObj_j.tagsDict["fftype"] 
                    #for i in range(NA):
                    #    print i,ASYMB[i],ATYPE[i]
                    sys.exit('charge group error')
                    
            if (ptclObj_i.tagsDict["fftype"] == 'LP' ):
                j_set = 0
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] != 0 ):
                        n_groups = ptclObj_j.tagsDict["qgroup"]
                        j_set = 1
                if( j_set ):
                    ptclObj_i.tagsDict["qgroup"] = n_groups

                else :
                    print ' lone pair connected to unset atom '
                    print ptclObj_i.tagsDict["symbol"] , " - ",ptclObj_j.tagsDict["symbol"] 
                    print ptclObj_i.tagsDict["fftype"] , " - ",ptclObj_j.tagsDict["fftype"] 
                    sys.exit('charge group error')
                    

		
    debug = 0
    if( debug ):
        for pid_i, ptclObj_i  in strucC.ptclC:
            # If unset set chrge group
             print i,ASYMB[i],CG_SET[i],ptclObj_i.tagsDict["qgroup"],  RING_NUMB[i]
	     
	     
    # Check for unset atoms 
    for pid_i, ptclObj_i  in strucC.ptclC:
        if( ptclObj_i.tagsDict["qgroup"] == 0 ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            
            print pid_i,ptclObj_i.tagsDict["symbol"],ptclObj_i.tagsDict["fftype"],NNAB,ELCNT[6],ELCNT[1],' not set '
	    sys.exit('charge group error')
            #else :
            # print ASYMB[i],ATYPE[i],NNAB,ELCNT[6],ELCNT[1],GTYPE[i],' set '

    return strucC

def qgroup_maxvol(strucC):
    """
    Find the maximum volume of the charge groups
    """

    qgrp_max = 0 
    for pid_i, ptclObj_i  in strucC.ptclC:
        if( ptclObj_i.tagsDict["qgroup"] > qgrp_max ):
            qgrp_max = ptclObj_i.tagsDict["qgroup"]


    print "    Charge groups ",qgrp_max
    a_max = 0
    a_min = 1e16

    dr_x_max = -1e16
    dV_max = -1e16

    for q_g in range( 1,qgrp_max+1 ):
        atom_cnt = 0

        r_x_min = 1e16 
        r_y_min = 1e16
        r_z_min = 1e16
        r_x_max = -1e16 
        r_y_max = -1e16
        r_z_max = -1e16

        for pid_i, ptclObj_i  in strucC.ptclC:
            if( ptclObj_i.tagsDict["qgroup"] == q_g ):
                atom_cnt += 1

                #print R[i]

                if( ptclObj_i.position[0] < r_x_min  ): r_x_min =  ptclObj_i.position[0]
                if( ptclObj_i.position[1] < r_y_min  ): r_y_min =  ptclObj_i.position[1] 
                if( ptclObj_i.position[2] < r_z_min  ): r_z_min =  ptclObj_i.position[2] 
                if( ptclObj_i.position[0] > r_x_max  ): r_x_max =  ptclObj_i.position[0] 
                if( ptclObj_i.position[1] > r_y_max  ): r_y_max =  ptclObj_i.position[1] 
                if( ptclObj_i.position[2] > r_z_max  ): r_z_max =  ptclObj_i.position[2] 


        if( atom_cnt == 0 ): print "      Group ",q_g," has 0 atoms "

        if( atom_cnt > a_max ): a_max = atom_cnt
        if( atom_cnt < a_min ): a_min = atom_cnt

        dr_x = r_x_max - r_x_min
        dr_y = r_y_max - r_y_min
        dr_z = r_z_max - r_z_min

        if( dr_x > dr_x_max ): dr_x_max = dr_x 
        if( dr_y > dr_x_max ): dr_x_max = dr_y
        if( dr_z > dr_x_max ): dr_x_max = dr_z 

        # Cubic PBC's
        #dr_x_pbc,dr_y_pbc,dr_z_pbc = prop.pbc_r_c(dr_x,dr_y,dr_z,LV)	    

        dV = dr_x*dr_y*dr_z
        if( dV > dV_max ): dV_max = dV

        #print "      Group ",q_g," has ",atom_cnt," atoms in volume ",dV," A^3 "

    print "    Max atoms in group ",a_max
    print "    Min atoms in group ",a_min
    print "    Max dr and dV ",dr_x_max, dV_max


def find_conections(strucC,cov_nblist, cov_nbindx, ring_nblist, ring_nbindex ):
    """
    Find inter conjugated ring  atoms 
    """
    import sys, numpy 
           
    RING_CONNECT = []
    debug = 0
    NA = len(strucC.ptclC)
    ADDED =  numpy.zeros(NA+1, dtype=int )    
    # Find Ring linkers
    for pid_i, ptclObj_i  in strucC.ptclC:        
        if ( ptclObj_i.tagsDict["number"] == 6 ) :
            NNAB = calc_nnab(pid_i,cov_nbindx)
            if( NNAB == 3 ):  # if sp2
                ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
                N_o = cov_nbindx[pid_i]
                N_f = cov_nbindx[pid_i+1] - 1
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if ( ADDED[pid_j] == 0 ):
                        NNAB_j = calc_nnab(pid_j,cov_nbindx)
                        if( NNAB_j == 3 ):  # if sp2
                            if(   ptclObj_i.tagsDict["ring"] !=  ptclObj_j.tagsDict["ring"] ): # ring linker
                                ADDED[pid_i] = 1
                                ADDED[pid_j] = 1
                                if(debug):
                                    print i,' and ',j,' link ',ELN[i],ELN[j]
                                    print ' ', RING_NUMB[i] , RING_NUMB[j]
                                RING_CONNECT.append( [pid_i,pid_j])

    return RING_CONNECT


def set_cply_tags(verbose, strucC,cov_nblist, cov_nbindx):
    #
    #  Set tags for new cply file 
    #     Use ctype tag and bonding enviroment
    #
    
    debug = True
    
    if( verbose): print "    setting cply tags "
    for pid_i, ptclObj_i  in strucC.ptclC:

        ptclObj_i.tagsDict["cplytag"] = ""
    for pid_i, ptclObj_i  in strucC.ptclC:

        if( debug ):
            print "CHECKING CTYPE ",pid_i,ptclObj_i.tagsDict["linkid"],ptclObj_i.tagsDict["cplytag"]
	
	# Set terminal attached carbons 
	if( ptclObj_i.tagsDict["linkid"] == "T" and ptclObj_i.tagsDict["number"] == 6 ):
	    #
	    term_con_cnt = 0 # Terminal connections count
	    #
	    # Find terminal hydrogen
	    #		    
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
	    for j_indx in range( N_o,N_f+1):
                pid_j = cov_nblist[j_indx]
                ptclObj_j = strucC.ptclC[pid_j]
                if( ptclObj_j.tagsDict["linkid"] == "X" and ptclObj_j.tagsDict["number"] == 1 ):
		    # Check to see if it has been bonded to another unit 
		    if( ptclObj_j.tagsDict["residue"] ==  ptclObj_j.tagsDict["residue"] ):
			term_con_cnt += 1
			ptclObj_j.tagsDict["cplytag"] = "termcap_H(" + str(pid_j) + ")_on_C("+str(pid_i)+")"
		      
	    if( term_con_cnt == 1 ):
		ptclObj_i.tagsDict["cplytag"] = "term_C(" + str(pid_i) + ")"
	    if( term_con_cnt > 1 ):
		print " Number of terminal atoms attached to atom ",pid_i," greater than 1 "
		sys.exit(" Error in terminal connections ")
	    
	# Set functional attached carbons 
	if( ptclObj_i.tagsDict["linkid"] == "F" and ptclObj_i.tagsDict["number"] == 6 ):
	    #
	    term_con_cnt = 0 # Terminal connections count
	    # 
	    # Find function hydrogen
	    #
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
	    for j_indx in range( N_o,N_f+1):
                pid_j = cov_nblist[j_indx]
                ptclObj_j = strucC.ptclC[pid_j]
                if( ptclObj_j.tagsDict["linkid"] == "R" and ptclObj_j.tagsDict["number"] == 1 ):
		    if( ptclObj_j.tagsDict["residue"] ==  ptclObj_j.tagsDict["residue"] ):
			term_con_cnt += 1
			ptclObj_j.tagsDict["cplytag"] = "funccap_H(" + str(pid_j) + ")_on_C("+str(pid_i)+")"

	    if( term_con_cnt == 1 ):
		ptclObj_i.tagsDict["cplytag"] = "func_C(" + str(pid_i) + ")"
		
	    if( term_con_cnt > 1 ):
		print " Number of functional atoms attached to atom ",pid_i," not equal to 1 "
		sys.exit(" Error in functional connections ")
		
    return strucC	

def zero_unitq(verbose,strucC,cov_nblist, cov_nbindx,zero_term,zero_func):
    import numpy
    import sys 
    
    #
    # Sum exsessive charges into carbon atoms 
    #
    for pid_i, ptclObj_i  in strucC.ptclC:                
	# Find each terminal hydrogen
	if(  ptclObj_i.tagsDict["linkid"]  == "X"  and zero_term  ):
	    term = pid_i 
	    # Check to be sure hydrogen
	    if( ptclObj_i.tagsDict["number"]  != 1 ):
		print " Non hydrogen used as terminal group atom ",pid_i  
		sys.exit(" Code unable to process multi atom (nonhyrdogen) terminal group ")
	    if( verbose ):
		print " Terminal atom found ",pid_i," ",ptclObj_i.tagsDict["fftype"]
	    #
	    # Loop over nieghbors to find attached atom
	    #
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
	    for j_indx in range( N_o,N_f+1):                
                pid_j = cov_nblist[j_indx]
                ptclObj_j = strucC.ptclC[pid_j]
                
		term_con_cnt = 0 # Terminal connections count
		if( ptclObj_j.tagsDict["linkid"].strip() == "T" ):
		    term_con_cnt += 1
		    if( verbose ):
			print " Terminal connection found ",pid_j," ",ptclObj_j.tagsDict["fftype"]
		    term_con = pid_j
	    # Check to be sure multiple atoms not found
	    if( term_con_cnt < 1 ):
		print " No terminal connections found "
		sys.exit(" Error in terminal connections ")
	    if( term_con_cnt > 1 ):
		print " Multiple terminal connections found "
		sys.exit(" Error in terminal connections ")
		
	    # Sum charges into base monomer unit
            strucC.ptclC[term_con].charge = strucC.ptclC[term_con].charge  + strucC.ptclC[term].charge 
	    strucC.ptclC[term].charge   = 0.0
	    
	# Find each functional hydrogen
	if( ptclObj_i.tagsDict["linkid"] == "R"  and zero_func ):
	    term = pid_i 
	    # Check to be sure hydrogen
	    if(  ptclObj_i.tagsDict["number"] != 1 ):
		print " Non hydrogen used as functional group "
		sys.exit(" Code unable to process multi atom (nonhyrdogen) functional group ")
	    if( verbose ):
		print " Functional atom found ",pid_i," ",ptclObj_i.tagsDict["fftype"]
	    #
	    # Loop over nieghbors to find attached atom
	    #
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
	    for j_indx in range( N_o,N_f+1):
                pid_j = cov_nblist[j_indx]
                ptclObj_j = strucC.ptclC[pid_j]
		term_con_cnt = 0 # Terminal connections count
		if( ptclObj_j.tagsDict["linkid"].strip() == "F" ):
		    term_con_cnt += 1
		    if( verbose ):
			print " functional connection found ",pid_j," ",ptclObj_j.tagsDict["fftype"]
		    term_con = pid_j
	    # Check to be sure multiple atoms not found
	    if( term_con_cnt  < 1 ):
		print " No functional connections found "
		sys.exit(" Error in functional connections ")
	    # Check to be sure multiple atoms not found
	    if( term_con_cnt > 1 ):
		print " Multiple functional connections found "
		sys.exit(" Error in functional connections ")
	    # Sum charges into base monomer unit
            strucC.ptclC[term_con].charge = strucC.ptclC[term_con].charge  + strucC.ptclC[term].charge 
	    strucC.ptclC[term].charge   = 0.0

    return strucC



'''

    def create_top(self,ff_charges): # Move out of class (or derived class)
        """
        Find topology information for force-field input files 
        """

        # Version 1 will be dependent on Atomicpy
        import elements , lammps ,gromacs , atom_types, top 

        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        GTYPE = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( [ float( ptclObj.position[0] ),float( ptclObj.position[1] ),float( ptclObj.position[2] )]  )  )
            AMASS.append( float( ptclObj.mass)  )
            CHARGES.append( float( ptclObj.charge ) )
            MOLNUMB.append( int( ptclObj.tagsDict["chain"] ) )
            RESID.append( int( ptclObj.tagsDict["residue"] ) )
            RESN.append(  ptclObj.tagsDict["resname"]  )
            GTYPE.append( ptclObj.tagsDict["gtype"] )
            
            
        # Direct copy of top.print_ff_files
            

        # Read in ff file
        #itp_file = 'ff.itp'
        ##print "   Read in parameters from ",itp_file
        #FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(itp_file)
                 
        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        verbose = True 

        # Find atomic number based on atomic symbol 
        ELN = elements.asymb_eln(ASYMB)
        AMASS = elements.eln_amass(ELN)

        #   Initialize  topology values
        #GTYPE = top.initialize_gtype( ELN )  # need to do this in frag read in 
        CHARN = top.initialize_charn( ELN )        

        NA = len(ELN)
        
        find_bonds = False 
        if( find_bonds ):

            #   Build covalent nieghbor list for bonded information 
            NBLIST, NBINDEX = top.build_covnablist(ELN,R)
            BONDS = top.nblist_bonds(NA,NBLIST, NBINDEX)

            print_bonds = True
            if( print_bonds ):
                for b_indx in range(len(BONDS)):
                    print " bond ",BONDS[b_indx][0]+1,BONDS[b_indx][1]+1

                sys.exit(" print bonds from proximity ")
        else:
            BONDS = []
            for b_i,bondObj in  self.bondC:
                pt_i = bondObj.pgid1
                pt_j = bondObj.pgid2
                BONDS.append( [pt_i-1,pt_j-1])
                print "create_top  bond ",pt_i,pt_j

            NBLIST,NBINDEX = groups.build_nablist_bonds(ELN,BONDS)
            #sys.exit(" passing bonds checking ")
            

        print_bonds = True
        if( print_bonds ):
            for b_indx in range(len(BONDS)):
                print " bond ",BONDS[b_indx][0]+1,BONDS[b_indx][1]+1
            
        
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
        #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)

        #
        # Set charge groups
        #
        CG_SET = []
        one = 1
        for i in range( len(ELN) ):
            CG_SET.append(one)
        #

        if( verbose ):
            print "      Finding atom types  "
            if( ff_charges ):
                print "      Using ff charges "

        d_mass = 0
        d_charge = 0

        find_rings = False 
        NA = len(ELN)
        if( find_rings ):
            RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
        else:
            zero = 0.0

            RINGLIST = []
            RINGINDEX = []
            RING_NUMB = []

            # relabel based on neighbors
            for i in range(NA):
                RINGLIST.append(zero)
                RING_NUMB.append(zero)
                RINGINDEX.append(zero)

            RINGLIST.append(zero)
            RINGINDEX.append(zero)

        # Asign oplsaa atom types
        ATYPE, CHARGES = atom_types.oplsaa( ff_charges,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB)

        #ATYPE , CHARGES = atom_types.biaryl_types( ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
        ATYPE ,RESID, CHARGES = atom_types.set_pmmatypes(ff_charges, ELN, ATYPE,GTYPE,RESID,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB )

        
        ATYPE,RESID,CHARGES,CG_SET,CHARN = atom_types.set_ptmatypes( ff_charges, ELN,ASYMB, ATYPE,GTYPE,RESID,CHARGES,AMASS,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB,CG_SET,CHARN )

        debug = False 
        if(debug):
            # relabel based on neighbors
            NA = len(ELN)
            for i in range(NA):
                print  i+1,ATYPE[i] ,RESID[i], CHARGES[i]
            sys.exit(" atom chagre check 1 ")

        #Refind inter ring types
        ATYPE , CHARGES  = atom_types.interring_types(ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )

        # Update structure information
        pt_cnt = 0
        for pid, ptclObj  in self.ptclC:
             ptclObj.tagsDict["fftype"] = ATYPE[pt_cnt]
             ptclObj.mass = AMASS[pt_cnt]
             ptclObj.charge = CHARGES[pt_cnt]
             ptclObj.tagsDict["ring"] = RING_NUMB[pt_cnt]
             pt_cnt += 1 

        # Add bonds to system
        #for i in range( len(BONDS) ):
        #    #
        #    a_i = BONDS[i][0] 
        #    a_j = BONDS[i][1]
        #    b_i = Bond( a_i+1, a_j+1 )
        #    self.bondC.put(b_i)

            #  Sudo code
            #   Read in parameter file
            #      gromacs itp file
            #      tinker parameter file
            #      lammps parameters in data file
            #   Find bonds
            #      from gaussian output optimization
            #      distance cut-off
            #         system_i = system_i.bonds()
            #   Find Rings
            #   Guess atom types
            #      amber
            #      oplsaa
            #         oligomer = oligomer.guess_oplsaatypes()

'''
"""
def find_types(struc_o,param_all):
    Find all the force-field parameters associated with a structure 

    ptclC_o =  struc_o.ptclC
    bondC_o  = struc_o.bondC
    angleC_o  = struc_o.angleC
    dihC_o  = struc_o.dihC
    
    #param_i = ParameterContainer()
    #
    # Create new parameter container to hold each unique type
    #   which exsists in the structure container struc_i
    #   this is necessary for parameter outputs to no have
    #   redundent values 
    #
    paramC_o = ParameterContainer()
    
    ljtypC_p =  paramC_p.ljtypC
    btypC_p =  paramC_p.btypC
    atypC_p =  paramC_p.atypC
    dtypC_p =  paramC_p.dtypC
    
    #
    # Count atom types
    #
    debug = 1
    for pid_o, ptclObj_o  in ptclC_o:    
        new_type = 1
        lpid_p = 0
        for lpid_p, ptclObj_p  in ptclC_p:    
            if( ptclObj_o.tagsDict["fftype"] == ptclObj_p.tagsDict["fftype"] ):
                new_type = 0
                ptclObj_o.type = str( lpid_p )
                if(debug): print ' maches previous ',ptclObj_p.tagsDict["fftype"]
        if( new_type ):
            ptclC_p.put(ptclObj_o)
            ptclObj_o.type = str( lpid_p + 1 )
            if(debug): print ' new type found ',ptclObj_o.tagsDict["fftype"]

    debug = 0
    if(debug):
        print " types found %d "%(len(ptclC_p))
        for lpid_p, ptclObj_p  in ptclC_p:
            print ptclObj_p.tagsDict["fftype"],ptclObj_p.type
        print "  All particles should have new type labeled as interger stored as a string "
        for pid_o, ptclObj_o  in ptclC_o:
            print ptclObj_o.tagsDict["fftype"],ptclObj_o.type
        sys.exit('find_types')



    #
    # Count bonds types
    #
    debug = 1
    if( debug ):
        print " len( bondC_o)",len( bondC_o)
    if( len( bondC_o) > 0 ):
        for b_o,bondObj_o in bondC_o:
            new_type = 1
            b_p = 0
            AT_i =  ptclC_o[ bondObj_o.pgid1 ].tagsDict["fftype"]
            AT_j =  ptclC_o[ bondObj_o.pgid2 ].tagsDict["fftype"]
            for b_p,bondObj_p in bondC_p:
                p_i = ptclC_o[ bondObj_p.pgid1 ].tagsDict["fftype"]
                p_j = ptclC_o[ bondObj_p.pgid2 ].tagsDict["fftype"]
                if ( AT_i == p_i  and  AT_j == p_j ):
                    new_type = 0
                    bondObj_o.type = b_p
                    break
                if ( AT_j == p_i  and  AT_i == p_j ):
                    new_type = 0
                    bondObj_o.type = b_p
                    break
                
            # If new 
            if( new_type ):
                bondC_p.put(bondObj_o)
                bondObj_o.type = str( b_p + 1 )
                #if(debug): print ' New bond type %s - %s  type -> %d '%(AT_i ,  AT_j, b_p + 1 )
                
            if(debug):
                print 'Bond type %s - %s  type -> %s '%(AT_i ,  AT_j,bondObj_o.type  )

    if(debug): sys.exit('BTYPE_IND')
    #
    # Count angles types
    #
    if( len(angleC_o) > 0):
        # 
        for a_o,angleObj_o in angleC_o:
            new_type = 1
            AT_k =  ptclC_o[ angleObj_o.pgid1 ].tagsDict["fftype"]
            AT_i =  ptclC_o[ angleObj_o.pgid2 ].tagsDict["fftype"]
            AT_j =  ptclC_o[ angleObj_o.pgid3 ].tagsDict["fftype"]
            a_p = 0
            for a_p,angleObj_p in angleC_p:
                p_k = ptclC_o[ angleObj_p.pgid1 ].tagsDict["fftype"]
                p_i = ptclC_o[ angleObj_p.pgid2 ].tagsDict["fftype"]
                p_j = ptclC_o[ angleObj_p.pgid3 ].tagsDict["fftype"]
                if ( AT_k == p_k and  AT_i == p_i and  AT_j == p_j ):
                    new_type = 0
                    angleC_p.put(angleObj_o)
                    break
    
                if ( AT_j == p_k and  AT_i == p_i and  AT_k == p_j ):
                    new_type = 0
                    angleC_p.put(angleObj_o)
                    break
    
            # If new 
            if( new_type ):
                angleC_p.put(angleObj_o)
                bondObj_o.type = str( a_p + 1 )
                if(debug):
                    print  ' new angle type ',AT_i ,  AT_j , AT_k

"""
