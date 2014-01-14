# group atoms

# Dr. Travis Kemper
# NREL
# Initial Date 9/2013
# travis.kemper@nrel.gov

# length - angstroms
# mass   - AMU
# volume - angstroms^3




def build_nablist(options,ELN,R,LV):
    # General neighbor list based on a cut off
    
    import sys,  numpy
    import elements, prop
    import datetime

    debug = 0
    p_time = 0

    #radi_cov =  elements.covalent_radi()

    NNAB = -1
    
    r_cut = 25.0
    sq_r_cut = r_cut*r_cut
        
    na = len( ELN )
    
    # Calculate
    
    den_buffer = 1.25 
    vol_cut = 4*numpy.pi/3*(r_cut**3)
    nbs_atom = vol_cut*options.n_den*den_buffer
    
    if( debug ):
	print " neighbors per atom ",vol_cut,nbs_atom
	#sys.exit("debug")
    
    maxnnab = na*nbs_atom

    cov_buffer = 1.25 
    
    radi_cov =  []
    NBLIST = numpy.empty( maxnnab,  dtype=int )
    NBINDEX = numpy.empty( maxnnab,  dtype=int )

    if( debug):
	print " R cut ",r_cut,sq_r_cut

    if( p_time ): t_i = datetime.datetime.now()

    for atom_i in range(na-1):
        NBINDEX[atom_i] = NNAB + 1
        #el_i = ELN[atom_i]
        r_i =  R[atom_i] 
	#rc_i = radi_cov[el_i]*cov_buffer
        for atom_j in range( atom_i+1, na ):
            if( atom_i != atom_j ):
                #el_j = ELN[atom_j]
                r_j =  R[atom_j] 
		sq_r_ij = prop.sq_drij_c(options,r_i,r_j,LV)
		#mag_dr =  numpy.linalg.norm(r_ij)
                #r_ij = delta_r(r_i,r_j)
                #r_cov = rc_i + radi_cov[el_j]*cov_buffer
		#sq_r_cov = sq_r_cov*sq_r_cov
                if( sq_r_ij <= sq_r_cut ):
                    NNAB += 1
                    NBLIST[NNAB] =  atom_j
                    if( debug ):
			print atom_i,atom_j , NBINDEX[atom_i] , NNAB 
			print r_i
			print r_j
			print sq_r_ij
                        #print ' atom i/j ', atom_i,atom_j,el_i,el_j
                        #print ' cov radi ', radi_cov[el_i] , radi_cov[el_j]
                        #print '   r_ij ',r_ij
                        #print '   r_cov ',r_cov
                    
    if( p_time ):
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	print "  build_nablist dt ",dt_min,dt_sec

    if(debug): sys.exit('debug')
    
    # Account for final atom, which has no connections 
    NBINDEX[atom_i+1] =  NNAB + 1
    NBINDEX[atom_i+2] =  NNAB + 1


    debug = 0
    if ( debug ):
        for i in range(len(ELN)):
            N_o = NBINDEX[ i  ]
            N_f = NBINDEX[ i + 1 ] - 1
            NNAB = N_f - N_o + 1
            print ' atom ', i,' ',ELN[i],NNAB,N_o,N_f    
            #
            # Find number of elements
            #
            #ELCNT = numpy.zeros(120)
            for indx in range( N_o,N_f+1):
                j = NBLIST[indx]
                print ELN[j],j
                #    el_j = ELN[j]
                #    ELCNT[j] = ELCNT[j] + 1

        sys.exit('debug')


    return (NBLIST, NBINDEX)


def index_nablist(options,ELN,R,LV):
    # neighbor list for certian atoms
    
    import sys,  numpy
    import elements, prop
    import datetime

    debug = 0
    p_time = 0

    #radi_cov =  elements.covalent_radi()

    NNAB = -1
    
    r_cut = 25.0
    sq_r_cut = r_cut*r_cut
        
    na = len( ELN )
    
    # Calculate
    
    den_buffer = 1.25 
    vol_cut = 4*numpy.pi/3*(r_cut**3)
    nbs_atom = vol_cut*options.n_den*den_buffer
    
    if( debug ):
	print " neighbors per atom ",vol_cut,nbs_atom
	#sys.exit("debug")
    
    maxnnab = na*nbs_atom

    cov_buffer = 1.25 
    
    radi_cov =  []
    NBLIST = numpy.empty( maxnnab,  dtype=int )
    NBINDEX = numpy.empty( maxnnab,  dtype=int )

    if( debug):
	print " R cut ",r_cut,sq_r_cut

    if( p_time ): t_i = datetime.datetime.now()

    for atom_i in range(na-1):
        NBINDEX[atom_i] = NNAB + 1
        #el_i = ELN[atom_i]
        r_i =  R[atom_i] 
	#rc_i = radi_cov[el_i]*cov_buffer
        for atom_j in range( atom_i+1, na ):
            if( atom_i != atom_j ):
                #el_j = ELN[atom_j]
                r_j =  R[atom_j] 
		sq_r_ij = prop.sq_drij_c(options,r_i,r_j,LV)
		#mag_dr =  numpy.linalg.norm(r_ij)
                #r_ij = delta_r(r_i,r_j)
                #r_cov = rc_i + radi_cov[el_j]*cov_buffer
		#sq_r_cov = sq_r_cov*sq_r_cov
                if( sq_r_ij <= sq_r_cut ):
                    NNAB += 1
                    NBLIST[NNAB] =  atom_j
                    if( debug ):
			print atom_i,atom_j , NBINDEX[atom_i] , NNAB 
			print r_i
			print r_j
			print sq_r_ij
                        #print ' atom i/j ', atom_i,atom_j,el_i,el_j
                        #print ' cov radi ', radi_cov[el_i] , radi_cov[el_j]
                        #print '   r_ij ',r_ij
                        #print '   r_cov ',r_cov
                    
    if( p_time ):
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	print "  build_nablist dt ",dt_min,dt_sec

    if(debug): sys.exit('debug')
    
    # Account for final atom, which has no connections 
    NBINDEX[atom_i+1] =  NNAB + 1
    NBINDEX[atom_i+2] =  NNAB + 1


    debug = 0
    if ( debug ):
        for i in range(len(ELN)):
            N_o = NBINDEX[ i  ]
            N_f = NBINDEX[ i + 1 ] - 1
            NNAB = N_f - N_o + 1
            print ' atom ', i,' ',ELN[i],NNAB,N_o,N_f    
            #
            # Find number of elements
            #
            #ELCNT = numpy.zeros(120)
            for indx in range( N_o,N_f+1):
                j = NBLIST[indx]
                print ELN[j],j
                #    el_j = ELN[j]
                #    ELCNT[j] = ELCNT[j] + 1

        sys.exit('debug')


    return (NBLIST, NBINDEX)


def build_nablist_bonds(ELN,BONDS):
    import sys,numpy

    debug = 0
    NNAB  = 0
    
    maxnnab = len(BONDS)*2 + 1
    #NBLIST = numpy.empty( maxnnab,  dtype=int )
    #NBINDEX = numpy.empty( maxnnab,  dtype=int )
    
    # python style nieghbor list
    nblist_py = [] #numpy.empty( maxnnab,  dtype=int )
    
    NBLIST = []
    NBINDEX = []
    NBLIST.append( 0 )
    
    for atom_i in range( len(ELN)):
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


def tempo( ATYPE,ELN,NBLIST,NBINDEX ,options ):
    import sys,numpy, top 
    
    debug = 0

    NA = len(ATYPE)
        
    a_gcnt = -1
    g_cnt = 0 
    group_index = []
    group_list = []

    a_added = [1]*NA
    group_numb = [0]*NA
    
    # Find TEMPO rings 
    for  atom_i in range( NA ):
        if(  ATYPE[atom_i] == 'NN' ):            
            N_o = NBINDEX[atom_i]
            N_f = NBINDEX[atom_i+1] - 1
            # nirtoxide published values for tempo 
            nitroxy_carb = 0
            for indx in range( N_o,N_f+1):
                atom_j = NBLIST[indx]
                if ( ATYPE[atom_j] == 'ON' ):
                    nitroxy_carb = 1
		    g_cnt += 1 
            if( nitroxy_carb ):
                if( a_added[atom_i] ):
                    group_index.append( a_gcnt + 1 )
                    a_added[atom_i] = 0
		    group_numb[atom_i] = g_cnt
                    group_list.append( atom_i )
                    a_gcnt = a_gcnt + 1                
                for indx in range( N_o,N_f+1):
                    atom_j = NBLIST[indx]
                    if( a_added[atom_j] ):
                        a_added[atom_j] = 0
			group_numb[atom_j] = g_cnt
                        group_list.append( atom_j )
                        a_gcnt = a_gcnt + 1   
                    Nj_o = NBINDEX[atom_j]
                    Nj_f = NBINDEX[atom_j+1  ] - 1
                    for indx_j in range(Nj_o,Nj_f+1):
                        atom_k = NBLIST[indx_j]  
                        if( a_added[atom_k] ):
                            a_added[atom_k] = 0
			    group_numb[atom_k] = g_cnt
                            group_list.append( atom_k )
                            a_gcnt = a_gcnt + 1    
                        Nk_o = NBINDEX[atom_k]
                        Nk_f = NBINDEX[atom_k+1  ] - 1
                        for indx_k in range(Nk_o,Nk_f+1):
                            atom_l = NBLIST[indx_k]  
                            if( a_added[atom_l] ):
                                a_added[atom_l] = 0
				group_numb[atom_l] = g_cnt
                                group_list.append( atom_l )
                                a_gcnt = a_gcnt + 1

    group_index.append( a_gcnt + 1 )
    
    if( options.verbose ):
        print "      TEMPO groups found ",len(group_index) - 1, g_cnt 
        
    if( debug):        
        for g_indx in range( len(group_index) - 1):
            N_o = group_index[g_indx]
            N_f = group_index[g_indx + 1 ] - 1
            NA_g = N_f - N_o + 1
            if( options.verbose): print "        Group ",g_indx, "with ",NA_g,"atoms "
            # Find atomic formula 
            ELCNT = top.calc_elcnt(g_indx,ELN,group_list,group_index) 
            for el_indx in range( len( ELCNT) ):
                if( ELCNT[el_indx] and debug ):  print "  ",el_indx,ELCNT[el_indx]
                
            for a_indx in range( N_o,N_f):
                i = group_list[a_indx]
                #print ATYPE[i]
                #print ASYMB[i],R[i][0],R[i][1],R[i][2])
	sys.exit(" debug list")
        
    return ( group_index,group_list,group_numb,g_cnt)


def ff_opt_groups( group_index_i,group_list_i, ELN,ASYMB,R,ATYPE,GTYPE,CHARGES,CHARN,AMASS,RESID,RESN,BONDS,ANGLES,DIH ,LV,NBLIST,NBINDEX,options):
    import os, sys
    import  gaussian, top, gromacs ,lammps , file_io, prop, atom_types, gaussian
    import numpy as np
    import datetime

    # Set debug 
    debug = 0
    p_time = 0
    
    # Read in parameter file
    
    FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(options.ff_file)

    # Intialize optimize as positions as the input positions 
    R_opt = R
 
    # Loop over all groups 
    for g_indx in range( len(group_index_i) - 1):
	
        # Get group index range
        N_o = group_index_i[g_indx]
        N_f = group_index_i[g_indx + 1 ] - 1
        NA_g = N_f - N_o + 1
        if( options.verbose ): print "        Printing gaussian input for group ",g_indx," with atoms ",NA_g
	
        F_id = "g_" + str(g_indx) 

        # Record id
        group_id = F_id
	
	
        # Place group into local array
	N_i = []
        ELECTRONS_i = 0 
        for a_indx in range( N_o,N_f+1):
            i = group_list_i[a_indx]
	    if( options.rm_VS  ):
		if( int(ELN[i]) > 0 ):
		    ELECTRONS_i += ELN[i]
		    N_i.append( i )
		    
		    print i,ELN[i]
	    else:
		ELECTRONS_i += ELN[i]
		N_i.append( i )
        
	# Extract group topology information
	ELN_j,ASYMB_j,R_j,ATYPE_j,GTYPE_j,CHARGES_j,RESID_j,RESN_j,CHARN_j,AMASS_j,BONDS_j,ANGLES_j,DIH_j,IMPS_j,REF_j,GHOST_j = top.pass_i(N_i,ELN,ASYMB,R,ATYPE,GTYPE,CHARGES,CHARN,AMASS,RESID,RESN,BONDS,ANGLES,DIH,options)
	
	if( options.zero_q ):    
	    for atom_i in range( len(CHARGES_j)):
		CHARGES_j[atom_i] = 0.0 
	    

	#print_oniom(  id_name, ASYMB,R,ATYPE,CHARGES,ELECTRONS_i,Q,M,FIX,ONMTAG,options)
	if( options.software == "gromacs"):

	    # Create input files
	    g_id = "min_"+str(g_indx)
	    g_top = g_id + '.top'
	    g_gro = g_id + '.gro'
	    g_mdp = g_id + '.mdp'
	    g_log = g_id + '.log'
	    rm_bck = "rm \"#\"*"
	    s_suf = options.gromacs_sufix + '_d'

	    if( file_io.file_exists(g_top) ): os.remove(g_top)
	    if( file_io.file_exists(g_gro) ): os.remove(g_gro)
	    if( file_io.file_exists(g_mdp) ): os.remove(g_mdp)
	    if( file_io.file_exists(g_log) ): os.remove(g_log)
	    os.system(rm_bck)

	    gromacs.print_gro(g_gro,GTYPE_j,RESID_j,RESN_j,R_j,LV)

	    # Identify total number of atom types for lammps output 
	    ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ATYPE_j,AMASS_j,BONDS_j,ANGLES_j,DIH_j)

	    BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
	    ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
	    DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,options)	    
	    IMPTYPE_F  = top.imp_parameters()

	    const = []
	    angle = []

	    gromacs.print_top( g_top,ASYMB_j , ELN_j,ATYPE_j, GTYPE_j, CHARN_j , CHARGES_j, AMASS_j,RESN_j, RESID_j ,BONDS_j , ANGLES_j , DIH_j , IMPS_j
		    ,const,angle
		    ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
		    ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)

	    gromacs.print_min_nopbcs(g_mdp)

	    #sys.exit(" debug gro ")

	    # Run minimization 
	    ff_energy_f = gromacs.run_gromacs(g_gro,g_top,g_mdp,s_suf,options )

	    # Get optimize positions 
	    R_k = gromacs.get_coord(g_id,options)
	
	elif( options.software == "gaussian" ):
	    
	    # Create input files
	    g_id = "min_"+str(g_indx)
	    g_com = g_id + '.com'
	    g_chk = g_id + '.chk'
	    g_fchk = g_id + '.fchk'
	    g_log = g_id + '.log'

	    if( file_io.file_exists(g_com) ): os.remove(g_com)
	    if( file_io.file_exists(g_chk) ): os.remove(g_chk)
	    if( file_io.file_exists(g_fchk) ): os.remove(g_fchk)
	    if( file_io.file_exists(g_log) ): os.remove(g_log)
            
	
	    ELECTRONS_j = 0
	    for atom_j in range(len(ELN_j)):
		ELECTRONS_j += ELN_j[atom_j] 
	    
	    gaussian.print_com( g_id, ASYMB_j,R_j,ATYPE_j,CHARGES_j,ELECTRONS_j,options)
	    gaussian.run(options,g_id)
	    
	    fchk_file = g_id +"/"+g_fchk
	    
	    NA, ELN_k, R_k, TOTAL_ENERGY_k = gaussian.parse_fchk( fchk_file )
	    
	#
	#elif( options.ff_software == "lammps" ):
	#    
	#    rest_file = 'rest.in'
	#    data_file = "out.data"
	#    lammps.print_rest(rest_file,data_file,DIH_CONST_ANGLE,DIH_ATOMS[cent_indx],options)
	#    
	#    # Print lammps structure file 
	#    data_file = "out.data"
	#    lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
	#       BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
	#       ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
	#       DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
	#       RESN,ATYPE_IND,CHARGES,R , ATYPE,
	#       BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LAT_CONST)
	#    
	# Update positions of non ghost atoms 
	for atom_j in range(len(ELN_j)):
	    if( GHOST_j[atom_j] == 0 ):
		atom_i = REF_j[atom_j]
		R[atom_i] = R_k[atom_j]
			
    return R


def cent_mass( atom_index, atom_list, R, AMASS ,options):
    import sys, numpy 
    
    debug = 0
    
    prop_dim = 3
    
    # Intialize center of mass list 
    r_mass = numpy.zeros( [ len(atom_index),prop_dim])
    
    if( debug ):
	size 
    
    # Loop over all groups 
    for g_indx in range( len(atom_index) - 1 ):
        # Intialize center of mass
        total_mass = 0.0
	
        # Loop over atoms in group 
        N_o = atom_index[g_indx ]
        N_f = atom_index[g_indx + 1 ] - 1
	
        for a_indx in range( N_o,N_f+1):
            i = atom_list[a_indx]
            a_mass_l = AMASS[i]
            # sum center of mass
            total_mass = total_mass + a_mass_l
            for d in range(prop_dim):
                r_mass[g_indx][d] = r_mass[g_indx][d]  + a_mass_l*R[i][d]
                
        # Normalize
        for d in range(prop_dim):
            r_mass[g_indx][d] = r_mass[g_indx][d]/total_mass
	    
	if(debug):
	    print "    g_i , cent_mass ",g_indx,r_mass[g_indx]
        
    return  r_mass
    

def nablist(group_cnt,group_cent,LV,options):
    import numpy 
    import prop
    
    debug = 0
    
    g_nb_ind = []
    g_nb_list = []
    g_r_list = []
    g_nb_cnt = -1
    
    g_cut = options.groups_com_cut
    g_cut_sq = g_cut*g_cut
#    
#    # Get fractional coordinates of each group's center of mass
#    f_cent = []
#    for g_i in range(group_cnt+1):
#	print g_i, " getting frac"
#	r_i = group_cent[g_i]
#	f_i = structure.r_frac(r_i,LV)
#	f_cent.append( numpy.array(f_i) )

    for g_i in range(group_cnt):
	g_nb_ind.append(g_nb_cnt + 1 )
	#f_i = f_cent[g_i]
	r_i = group_cent[g_i]
		    
        for g_j in range( g_i+1, group_cnt):
	    
	    r_j = group_cent[g_j]
	    
	    
	    sq_r_ij = prop.sq_drij_c(options,r_i,r_j,LV)
	    
	    #f_j = f_cent[g_j]
	    #f_ij = f_j - f_i
	    #df_pbc = structure.f_pbc(f_ij,options)
	    #dr_pbc = structure.f_real(df_pbc,LV)
	    #mag_dr =  numpy.linalg.norm(dr_pbc)
	    #if(  mag_dr < g_cut ):
		
            if( sq_r_ij <= g_cut_sq ):
		mag_dr = numpy.sqrt( sq_r_ij)
		g_nb_cnt += 1
		g_nb_list.append(g_j)
		g_r_list.append(mag_dr)
		
		if( debug):
		    print "  group ",g_i,g_j ,mag_dr

    g_nb_ind.append(g_nb_cnt + 1 )
    
    return (g_nb_ind,g_nb_list,g_r_list )
    