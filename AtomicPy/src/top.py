#! /usr/bin/env python
# Process topology information 

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov

# Conversions
KJ_KCAL = 0.23901
GRO_SIG = 5.61230943
NM_ANG = 10.0

def  initialize_gtype(ELN):
    import sys,elements
    
    GTYPE = []
    #
    elsymbol = elements.set_elsymbol()
    for atom_i in range( len(ELN) ):
        el_i = ELN[atom_i]
        atomic_symbol = elsymbol[el_i]
        GTYPE.append( atomic_symbol + str(atom_i+1) )
         
    return GTYPE

def initialize_resid( ELN ):

    RESID = []
    residue_id = 'CHN'
    
    for atom_i in range( len(ELN) ):
        RESID.append( residue_id )
        
    return RESID

def initialize_resn( ELN ):

    RESN = []
    residue_numb = 1
    
    for atom_i in range( len(ELN) ):
        RESN.append( residue_numb )
        
    return RESN
    
def initialize_charn( ELN ):

    CHARN = []
    char_numb = 1
    
    for atom_i in range( len(ELN) ):
        CHARN.append( char_numb )

    return CHARN

          
def initialize_charges( ELN ):

    CHARGES = []
    q = 0.0 
    
    for atom_i in range( len(ELN) ):
        CHARGES.append( q )

    return CHARGES

def calc_nnab(i,NBLIST,NBINDEX):
    #
    # Find number of elements 
    #
    N_o = NBINDEX[ i  ]
    N_f = NBINDEX[  i+1  ] - 1 
    NNAB = N_f - N_o + 1
    return NNAB

def calc_elcnt(i,ELN,NBLIST,NBINDEX):
    import numpy
    #
    # Find number of elements 
    #
    ELCNT = numpy.zeros(120, dtype =int )
    N_o = NBINDEX[ i  ]
    N_f = NBINDEX[  i+1  ] - 1 
    for indx in range( N_o,N_f+1):
        j = NBLIST[indx]
        el_j = int( ELN[j] )
	if( el_j >= 0 ):
	    ELCNT[el_j] = ELCNT[el_j] + 1

    return ELCNT 

def build_covnablist(ELN,R):
    #
    # Build covalent neighbor list from elements and positions 
    #
    import sys, elements, numpy 
    import datetime

    debug = 0
    p_time = 0
    
    na = len( ELN )
    maxnnab = na*12

    cov_buffer = 1.25
    
    radi_cov =  []
    NBLIST = numpy.empty( maxnnab,  dtype=int )
    NBINDEX = numpy.empty( maxnnab,  dtype=int )

    radi_cov =  elements.covalent_radi()


    NNAB = 0

    if( p_time ): t_i = datetime.datetime.now()

    for atom_i in range(na):
        NBINDEX[atom_i] = NNAB + 1
        el_i = ELN[atom_i]
        r_i = numpy.array( R[atom_i] )
	rc_i = radi_cov[el_i]*cov_buffer
	
	NNAB_i = NNAB 
        for atom_j in range( na ):
            if( atom_i != atom_j ):
                el_j = ELN[atom_j]
                r_j = numpy.array( R[atom_j] )
		r_ij = r_j - r_i
		mag_dr =  numpy.linalg.norm(r_ij)
                #r_ij = delta_r(r_i,r_j)
                r_cov = rc_i + radi_cov[el_j]*cov_buffer
		

    
                if( mag_dr <= r_cov ):
                    NNAB = NNAB + 1
                    NBLIST[NNAB] =  atom_j
						
		    if( debug ):
			print '    atom i/j ', atom_i+1,atom_j+1,el_i,el_j
			print '       cov radi ', radi_cov[el_i] , radi_cov[el_j]
			print '       r_i ',r_i
			print '       r_j ',r_j
			print '       r_ij ',r_ij 
			print '       |r_ij| ',mag_dr 
			print '       r_cov ',r_cov
	
		    
		    
		    
	if( debug ):
	    print "  atom ",atom_i + 1 ,ELN[atom_i], " has ",NNAB - NNAB_i  ," bonded nieghbors "
		    
	
    if( p_time ):
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	print "  build_nablist dt ",dt_min,dt_sec

    if(debug): sys.exit('debug build_covnablist')
    
    # Account for final atom position
    NBINDEX[atom_i+1] =  NNAB + 1


    debug = 0
    if ( debug ):
        for i in range(len(ELN)):
            N_o = NBINDEX[ i  ]
            N_f = NBINDEX[ i + 1 ] - 1
            NNAB = N_f - N_o + 1
            print ' atom ', i,' ',ELN[i],NNAB,N_o    
            #
            # Find number of elements
            #
            #ELCNT = numpy.zeros(120)
            for indx in range( N_o,N_f+1):
                j = NBLIST[indx]
                print ELN[j],j
                #    el_j = ELN[j]
                #    ELCNT[j] = ELCNT[j] + 1

        sys.exit('debug build_covnablist')


    return (NBLIST, NBINDEX)


def nblist_bonds(NA,NBLIST, NBINDEX):
    import sys
    #
    # Generate bonds from neighbor list
    #
    BONDS = []*2
 
    for atom_i in range(NA):
        N_o = NBINDEX[atom_i ]
        N_f = NBINDEX[atom_i+ 1 ] - 1
        for indx_j in range( N_o,N_f+1):
            atom_j = NBLIST[indx_j]
            if( atom_j > atom_i):
                BONDS.append([atom_i,atom_j])
            
    return BONDS

def nblist_angles(NA,NBLIST, NBINDEX):
    import sys
    #
    # Generate angles from neighbor list 
    #
    
    ANGLES = []*3
 
    for atom_i in range(NA):
        N_o = NBINDEX[atom_i ]
        N_f = NBINDEX[atom_i+ 1 ] - 1
        NNAB = N_f - N_o + 1
        if( NNAB >= 2 ):
            for indx_j in range( N_o,N_f):
                atom_j = NBLIST[indx_j]
                for indx_k in range( indx_j+1,N_f+1):
                    atom_k = NBLIST[indx_k]
                    if ( atom_j != atom_k ):
                        ANGLES.append([atom_j,atom_i,atom_k])

    return ANGLES


def nblist_dih(NA,NBLIST, NBINDEX):
    import sys
    #
    # Generate dihedrals from neighbor list 
    #
    DIH = []*4
 
    limit_dih = 1
    limit_n = 1 
 
    for atom_i in range(NA):
        N_o = NBINDEX[atom_i ]
        N_f = NBINDEX[atom_i+ 1 ] - 1
        NNAB = N_f - N_o + 1
        for indx_j in range( N_o,N_f+1):
            atom_j = NBLIST[indx_j]
            if( atom_j >atom_i):  # for double counting
		dih_ij_cnt = 0 
                for indx_k in range( N_o,N_f+1 ):
                    atom_k = NBLIST[indx_k]
                    if ( atom_k != atom_j ):
                        No_j =  NBINDEX[ atom_j  ]
                        Nf_j = NBINDEX[atom_j+ 1 ] - 1
                        for indx_l in range( No_j,Nf_j+1):
                            atom_l = NBLIST[indx_l]
                            if ( atom_l != atom_i and atom_l != atom_k ):
				dih_ij_cnt += 1
				if( limit_dih ):
				    if( dih_ij_cnt <= limit_n ):
					DIH.append([atom_k,atom_i,atom_j,atom_l])
				else:
				    DIH.append([atom_k,atom_i,atom_j,atom_l])

    return DIH

def nblist_imp(NA,NBLIST, NBINDEX,ELN):
    import sys
    #
    #
    # Generate impropers from neighbor list 
    #
    IMPS = []

    return  IMPS
        

def id_ring(a_i,ELN,NBLIST,NBINDEX):
    import sys
        
    debug = 0
    
    R_SET = []
    RING_ATOMS = []
    BAD_PATH = []

    one = 1
    zero = 0
    atoms_in_ring = 0
    # relabel based on neighbors
    NA = len(ELN)
    for i in range(NA):
        R_SET.append(one)
        
    R_SET[a_i] = 0 
    r_term = 1
    cnt = 0
    p_cnt = 0
    if(debug): print ' initializing ring ',a_i,ELN[a_i]

    last_i = a_i
    NNAB_last = calc_nnab(a_i,NBLIST,NBINDEX)
    while ( r_term ):
        N_o = NBINDEX[last_i]
        N_f = NBINDEX[last_i+1] - 1

        for n_indx in range( N_o,N_f+1):
            j = NBLIST[n_indx]
            cnt = cnt + 1
            if( cnt > NNAB_last+1 ):
                p_cnt = p_cnt + 1
                if(debug): print ' bad path found resetting '
                for ring_a in range(len(RING_ATOMS)):
                    j = RING_ATOMS[ring_a]
                    
                    
                BAD_PATH.append(j)
                RING_ATOMS = []
                for i in range(NA):
                    R_SET[i] = one
                atoms_in_ring = 0
                R_SET[a_i] = 0 
                cnt = 0
                last_i = a_i
                
                if(debug): print '  resetting last_i to ',a_i,ELN[a_i]
                    
                for bad_i in range( len(BAD_PATH)):
                    j = BAD_PATH[bad_i]
                    R_SET[j] = 0
                    if(debug): print '  bad path atoms ',ELN[j]
                break
                    
            if( atoms_in_ring  > 1 and j == a_i ):
                r_term = 0
                if(debug): print ' ring found with ',atoms_in_ring 
                break
                        
            NNAB_j = calc_nnab(j,NBLIST,NBINDEX)
            ELCNT_j = calc_elcnt(j,ELN,NBLIST,NBINDEX)
            ring_type = 0
	    
	    debug = 0
            if(debug):
		print '           with ', ELN[j]
		print '           with ', NNAB_j,
		print '           with ', R_SET[j]
		
            if ( ELN[j] == 6 and NNAB_j == 3 and R_SET[j] == 1 ): ring_type = 1 # ATYPE[j] == 'CA' ):
            if ( ELN[j] == 16 and NNAB_j == 2 and R_SET[j] == 1 and ELCNT_j[6] == 2 ): ring_type = 1 # ATYPE[j] == 'CA' ):
            if ( ELN[j] == 7 and NNAB_j >= 2 and R_SET[j] == 1 ): ring_type = 1 # ATYPE[j] == 'CA' ):
            if( ring_type ):
                atoms_in_ring = atoms_in_ring + 1
                RING_ATOMS.append(j)
                R_SET[j] = 0
                last_i = j
                cnt = 0
                NNAB_last = NNAB_j 
                if(debug): print '   atom ',j,ELN[j],' added '
                break
            
        if (p_cnt > 100 ):
            if(debug): print '     max paths considered ',100
            r_term = 0
            
    return RING_ATOMS

def find_rings(ELN,NBLIST,NBINDEX):
    #
    # Add hydrogens to ring 
    import sys
        
    RINGLIST = []
    RINGINDEX = []
    RING_NUMB = []

    debug = 0
    
    print
    print ' finding rings '
    print

    IN_RING = []
    
    one = 0 
    zero = 0 
    ring_cnt = 0
    a_cnt = 0
    # relabel based on neighbors
    NA = len(ELN)
    for i in range(NA):
        IN_RING.append(one)
        RINGLIST.append(zero)
        RING_NUMB.append(zero)
        RINGINDEX.append(zero)
        
    RINGLIST.append(zero)
    RINGINDEX.append(zero)
        
        
    # relabel based on neighbors
    for i in range(NA):
        RING_ATOMS = id_ring(i,ELN,NBLIST,NBINDEX)
            
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
        for i in range(NA):
            if( RING_NUMB[i] == r_numb ):
                a_cnt = a_cnt + 1
                RINGLIST[a_cnt] = i
                
    RINGINDEX[ring_cnt+1] = a_cnt + 1


    # Set attached H to have same ring #

    for r_numb in range(1,ring_cnt+1):
	N_o = RINGINDEX[r_numb]
	N_f = RINGINDEX[r_numb+1] - 1
	print ' Ring ', r_numb
	for r_indx in range(N_o,N_f+1):
	    i = RINGLIST[r_indx]
	    N_o = NBINDEX[i]
	    N_f = NBINDEX[i+1] - 1
	    for n_indx in range( N_o,N_f+1):
		j = NBLIST[n_indx]
		if( ELN[j] == 1 ):
		    RING_NUMB[j] = r_numb
		
		
	    

    debug = 0
    if(debug):
        print ring_cnt , " rings found "
        for r_numb in range(1,ring_cnt+1):
            N_o = RINGINDEX[r_numb]
            N_f = RINGINDEX[r_numb+1] - 1
            print ' Ring ', r_numb
            for r_indx in range(N_o,N_f+1):
                i = RINGLIST[r_indx]
                print ' type ',ELN[i],' or ',i,RING_NUMB[i]
                #sys.exit('find_rings')
                
        for atom_i in range( len(ELN)):
            print ' atom  ',atom_i,ELN[atom_i],RING_NUMB[atom_i]
	    
	#sys.exit('top.find_rings')

    return ( RINGLIST, RINGINDEX , RING_NUMB )

def find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB):
    import sys, numpy 
           
    RING_CONNECT = []
    debug = 0
    NA = len(ELN)
    ADDED =  numpy.zeros(NA, dtype=int )    
    # Find Ring linkers
    for i in range(NA):
        
        if ( ELN[i] == 6 ) :
          NNAB = calc_nnab(i,NBLIST,NBINDEX)
          if( NNAB == 3 ):  # if sp2
              ELCNT = calc_elcnt(i,ELN,NBLIST,NBINDEX)
              N_o = NBINDEX[ i  ]
              N_f = NBINDEX[  i+1  ] - 1
              for indx in range( N_o,N_f+1):
                  j = NBLIST[indx]
                  if ( ADDED[j] == 0 ):
                      NNAB_j = calc_nnab(j,NBLIST,NBINDEX)
                      if( NNAB_j == 3 ):  # if sp2
                          if(  RING_NUMB[i] !=  RING_NUMB[j] ): # ring linker
                              ADDED[i] = 1
                              ADDED[j] = 1
                              if(debug):
                                  print i,' and ',j,' link ',ELN[i],ELN[j]
                                  print ' ', RING_NUMB[i] , RING_NUMB[j]
                              RING_CONNECT.append( [i,j])

    return RING_CONNECT


def ring_natoms(r_numb,RINGINDEX):

    RN_o = RINGINDEX[r_numb]
    RN_f = RINGINDEX[r_numb+1] - 1
    nring = RN_f - RN_o + 1
    
    return nring 



def special_types(ATYPE, ASYMB , ELN , atom_types):
    import sys
    import file_io

    debug = 1

    # Open types file
    if( file_io.file_exists(atom_types )):
        if( debug ): print "Reading in ",atom_types
        F = open(atom_types , 'r' )
        Lines = F.readlines()
        F.close()

    else:
        print " Specified .gro file ",atom_types," does not exisit "
        sys.exit("Invalid file ")

    # Allocate types
    s_type = []
    s_numb = []
    s_symb = []

    # Loop over lines of top file
    for line in Lines:
        col = line.split()
        if ( len(col) > 1 and col[0][:1] != ';' ):
            s_type.append( col[0].strip() )
            s_numb.append( col[1].strip())
            s_symb.append( col[2].strip())


    for a_indx in range( len(ATYPE) ):
        a_type = ATYPE[a_indx]
        for s_indx in range( len( s_type)):
            if( s_type[s_indx] == a_type.strip() ):
                ELN[a_indx] = s_numb[s_indx]
                ASYMB[a_indx] = s_symb[s_indx]
		
		if(debug): print "   editing atom ",a_indx,a_type," to atomic# ",ELN[a_indx] 

    return ( ASYMB,ELN )


def set_chargegroups(verbose,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LV):
    import sys, prop
    
    if( verbose ):
	print "   Setting charge groups "
    
    NA = len(ATYPE) 
    debug = 0
    # initialize 

    n_groups = 0
    # set rings as charge groups

    debug = 0
    for i in range(NA):
        if( CG_SET[i] ):
            if( RING_NUMB[i] > 0 ):
                CHARN[i] = RING_NUMB[i]
                CG_SET[i] = 0 
                if(debug): print i,ASYMB[i],   RING_NUMB[i]
                if(  RING_NUMB[i] >  n_groups ):
                    n_groups = n_groups + 1
                
    if( debug): sys.exit('ring charge groups')
            
    # set methyls and alkyls 
    for i in range(NA):
        if( CG_SET[i] ):
            NNAB = calc_nnab(i,NBLIST,NBINDEX)
            ELCNT = calc_elcnt(i,ELN,NBLIST,NBINDEX)
            N_o = NBINDEX[i]
            N_f = NBINDEX[i+1] - 1
            if( ELN[i] == 6 and ELCNT[1] > 0 ):  # CH_n ( Methyls / alkyl ) 
                n_groups = n_groups + 1
                CG_SET[i] = 0
                CHARN[i] = n_groups
                for j_indx in range( N_o,N_f+1):
                    j = NBLIST[j_indx]
                    if ( ELN[j] == 1 ):
                        CHARN[j] = n_groups
                        CG_SET[j] = 0
			
            if( ATYPE[i] == 'C' ):  # 
                n_groups = n_groups + 1
                CG_SET[i] = 0
                CHARN[i] = n_groups
                for j_indx in range( N_o,N_f+1):
                    j = NBLIST[j_indx]
                    if ( ELN[j] == 8 ):
                        CHARN[j] = n_groups
                        CG_SET[j] = 0

            if( ATYPE[i] == 'ON' ):  # nitroxide
                n_groups = n_groups + 1
                CG_SET[i] = 0
                CHARN[i] = n_groups
                for j_indx in range( N_o,N_f+1):
                    j = NBLIST[j_indx]
                    if ( CG_SET[j]  ):
			CHARN[j] = n_groups
			CG_SET[j] = 0

    # check unset heavy atom groups
    for i in range(NA):
        if( CG_SET[i] ):
            NNAB = calc_nnab(i,NBLIST,NBINDEX)
            ELCNT = calc_elcnt(i,ELN,NBLIST,NBINDEX)
            N_o = NBINDEX[i]
            N_f = NBINDEX[i+1] - 1
	    
	    
            if( ELN[i] == 6 or ELN[i] == 8 ):  #
		q_g = 0 
                for j_indx in range( N_o,N_f+1):
		    j = NBLIST[j_indx]
                    if ( not CG_SET[j]  ):
			q_g = CHARN[j]
			break
		if( q_g  == 0 ):
			
		    n_groups = n_groups + 1
		    CG_SET[i] = 0
		    CHARN[i] = n_groups
		else:
		    
		    CG_SET[i] = 0
		    CHARN[i] = q_g
		    
		#print " unset carbon set to ",CHARN[i]
			

            if( ATYPE[i] == 'CA' ):  # conjugate
                n_groups = n_groups + 1
                CG_SET[i] = 0
                CHARN[i] = n_groups
 
                    

    # set light atom groups 
    for i in range(NA):
        if( CG_SET[i] ):
            NNAB = calc_nnab(i,NBLIST,NBINDEX)
            ELCNT = calc_elcnt(i,ELN,NBLIST,NBINDEX)
            N_o = NBINDEX[i]
            N_f = NBINDEX[i+1] - 1
            if ( ELN[i] == 1 or  ELN[i] == 9 ):
                j_set = 0
                for j_indx in range( N_o,N_f+1):
                    j = NBLIST[j_indx]
                    if( CG_SET[j] == 0 ):
                        n_groups = CHARN[j]
                        j_set = 1
                if( j_set ):
                    CHARN[i] = n_groups
                    CG_SET[i] = 0
                else :
                    print ASYMB[i],' connected to unset atom '
                    print ASYMB[i], " - ",ASYMB[j]
                    print ATYPE[j],ASYMB[j]
                    for i in range(NA):
                        print i,ASYMB[i],ATYPE[i]
                    sys.exit('charge group error')
                    
            if ( ASYMB[i] == 'LP' ):
                j_set = 0
                for j_indx in range( N_o,N_f+1):
                    j = NBLIST[j_indx]
                    if( CG_SET[j] == 0 ):
                        n_groups = CHARN[j]
                        j_set = 1
                if( j_set ):
                    CHARN[i] = n_groups
                    CG_SET[i] = 0
                else :
                    print ' lone pair connected to unset atom '
                    print ATYPE[j],ASYMB[j]
                    sys.exit('charge group error')
                    

		
    debug = 0
    if( debug ):
        for i in range(NA):
             print i,ASYMB[i],CG_SET[i],CHARN[i],  RING_NUMB[i]
	     
	     
    # Check for unset atoms 
    for i in range(NA):
        if( CG_SET[i] ):
            NNAB = calc_nnab(i,NBLIST,NBINDEX)
            ELCNT = calc_elcnt(i,ELN,NBLIST,NBINDEX)
            N_o = NBINDEX[i]
            N_f = NBINDEX[i+1] - 1    
            print i,ASYMB[i],ATYPE[i],NNAB,ELCNT[6],ELCNT[1],' not set '
	    sys.exit('charge group error')
            #else :
            # print ASYMB[i],ATYPE[i],NNAB,ELCNT[6],ELCNT[1],GTYPE[i],' set '
   
               
    # Check for unset atoms
    if ( verbose ):
	print "    Charge groups ",max(CHARN)
	a_max = 0
	a_min = 1e16
	
	dr_x_max = -1e16
	dV_max = -1e16

	for q_g in range( 1,max(CHARN)+1 ):
	    atom_cnt = 0
	    
	    r_x_min = 1e16 
	    r_y_min = 1e16
	    r_z_min = 1e16
	    r_x_max = -1e16 
	    r_y_max = -1e16
	    r_z_max = -1e16
	    
	    for i in range(NA):
		if( CHARN[i] == q_g ):
		    atom_cnt += 1
		    
		    #print R[i]
		    
		    if( R[i][0] < r_x_min  ): r_x_min = R[i][0]
		    if( R[i][1] < r_y_min  ): r_y_min = R[i][1] 
		    if( R[i][2] < r_z_min  ): r_z_min = R[i][2] 
		    if( R[i][0] > r_x_max  ): r_x_max = R[i][0] 
		    if( R[i][1] > r_y_max  ): r_y_max = R[i][1] 
		    if( R[i][2] > r_z_max  ): r_z_max = R[i][2] 
		    
	    
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
	
	
    for i in range(NA):
        if( CG_SET[i] ):
            NNAB = calc_nnab(i,NBLIST,NBINDEX)
            ELCNT = calc_elcnt(i,ELN,NBLIST,NBINDEX)
            N_o = NBINDEX[i]
            N_f = NBINDEX[i+1] - 1    
            print i,ASYMB[i],ATYPE[i],NNAB,ELCNT[6],ELCNT[1],' not set '
	    sys.exit('charge group error')
            #else :
            # print ASYMB[i],ATYPE[i],NNAB,ELCNT[6],ELCNT[1],GTYPE[i],' set '
   
    # set others
    debug = 0
    if(debug): sys.exit('set_chargegroups')

    return CHARN

def check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX):
    import sys, numpy
    
    ATYPE_NNAB = numpy.zeros( len(ATYPE_REF ) )

    NA = len(GTYPE)
    debug = 0
    error_on_bonding = 0
    for a_ind in range(  len( ATYPE_REF )  ):
        for i in range( NA ):
            if( ATYPE_IND[i] == a_ind ):
                NNAB = calc_nnab(i,NBLIST,NBINDEX)
                if( ATYPE_NNAB[a_ind] == 0 ):
                    ATYPE_NNAB[a_ind] = NNAB
                    if( debug): print ' For atom type ',ATYPE_REF[a_ind],a_ind,' has ',NNAB,' nieghbors ',GTYPE[i]
                elif (  ATYPE_NNAB[a_ind] != NNAB ):
                    print ' Inconsistent number of nieghbors for ',GTYPE[i],NNAB
                    if( error_on_bonding): sys.exit(' error in bonding ')
                    
    if(debug): sys.exit('check_types')

    return ATYPE_NNAB

def atom_parameters(itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES):
    import sys

    global GRO_SIG,KJ_KCAL,NM_ANG

    # parameters
    ATYPE_EP = []
    ATYPE_SIG = []

    ff_type = "oplsaa"
    
    if( ff_type == "oplsaa" ): combrule = 3 
    
    debug = 0
    #
    # Find atom type parameters
    #
    for ind in range( len(ATYPE_REF) ):
        if( debug): print ATYPE_REF[ind]
        AT_i =  ATYPE_REF[ind]

        check = 1
        for ff_i in range (len(FF_ATOMTYPES)):
            FF_i = FF_ATOMTYPES[ff_i] #.split() 
            if ( AT_i == FF_i[0] ):
		if( combrule == 2 ):
                    ATYPE_SIG.append( float( FF_i[5]) * GRO_SIG )
    		if( combrule == 3 ):
		    ATYPE_SIG.append( float( FF_i[5]) * NM_ANG )
                ATYPE_EP.append( float( FF_i[6]) *  KJ_KCAL )
                check = 0
                if( debug):  print  float( FF_i[5]) ,float( FF_i[6])         
                break
        if check :
            print ' atom type ',AT_i,'  not in ff ;  index ',ind
            print ' grep "',AT_i,'" ',itp_file
	    print 'unknow ff specification '
            sys.exit(' top.atom_parameters ')

    return (ATYPE_EP, ATYPE_SIG)

def bond_parameters(itp_file,BTYPE_IND , BTYPE_REF,FF_BONDTYPES):
    import sys
    global GRO_SIG,KJ_KCAL,NM_ANG

    # parameters
    BONDTYPE_F = []
    BONDTYPE_R0 = []
    BONDTYPE_K = []

    debug = 0

    
    # Find bond type parameters
    for ind in range( len( BTYPE_REF ) ):
        AT_i =  BTYPE_REF[ind][0]
        AT_j =  BTYPE_REF[ind][1]
        check = 1
        for ff_i in range (len(FF_BONDTYPES)):
            FF_l = FF_BONDTYPES[ff_i] #.split() 
            if ( AT_i == FF_l[0] and  AT_j == FF_l[1]   ):
                check = 0
                break
            if (  AT_i == FF_l[1] and  AT_j == FF_l[0]  ):
                check = 0
                break
            
        if check :
            print ' bond ',ind,' ',AT_i,' - ',AT_j,'  not in ff'
            print ' grep "',AT_i,'" ',itp_file,' | grep "',AT_j , ' " '
	    print 'unknow ff specification '	    
            sys.exit(' top.bond_parameters ')
        else:
            func_type = int(FF_l[2])
            BONDTYPE_F.append( func_type )
            if( func_type == 1 ):
                BONDTYPE_R0.append( float(FF_l[3]) * NM_ANG )
                BONDTYPE_K.append( float(FF_l[4]) *KJ_KCAL/NM_ANG/NM_ANG/2 )
            else :
                print ' Bond type ', func_type , ' unknown '
		print 'unknow ff specification '	    
		sys.exit(' top.bond_parameters ')
            
    return ( BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K )

def angle_parameters(itp_file,ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES):
    import sys
    global GRO_SIG,KJ_KCAL,NM_ANG

    # parameters
    ANGLETYPE_F = []
    ANGLETYPE_R0 = []
    ANGLETYPE_K = []

    debug = 0

    # Find angle type parameters
    for ind in range( len( ANGTYPE_REF ) ):
        AT_i =  ANGTYPE_REF[ind][0]
        AT_j =  ANGTYPE_REF[ind][1]
        AT_k =  ANGTYPE_REF[ind][2]
        check = 1
        for ff_i in range (len(FF_ANGLETYPES)):
            FF_l = FF_ANGLETYPES[ff_i] #.split() 
            if ( AT_i == FF_l[0] and  AT_j == FF_l[1] and  AT_k == FF_l[2]   ):
                check = 0
                break
            if ( AT_k == FF_l[0] and  AT_j == FF_l[1] and  AT_i == FF_l[2]   ):
                check = 0
                break
            
            
        if check :
            print ' angle ',ind,' ',AT_i,' - ',AT_j,' - ',AT_k,'  not in ff'
            print ' grep "',AT_i,'" ',itp_file,' | grep "',AT_j ,'" | grep "',AT_k,'"'
            print ' grep "',AT_i,'" ',itp_file,' | grep "',AT_j ,'" | grep "',AT_k,'"'
	    print 'unknow ff specification '	
            sys.exit(' top.angle_parameters ')
        else:
            func_type = int(FF_l[3])
            ANGLETYPE_F.append( func_type )
            if( func_type == 1 ):
                ANGLETYPE_R0.append( float(FF_l[4]) )
                ANGLETYPE_K.append( float(FF_l[5])*KJ_KCAL/2 )
            else :
                print ' Angle type ', func_type , ' unknown '
	        print 'unknow ff specification '	
		sys.exit(' top.angle_parameters ')
            
    return ( ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K )

def dih_parameters( itp_file,norm_dihparam,DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB ):
    import sys
    global GRO_SIG,KJ_KCAL,NM_ANG

    # parameters
    DIHTYPE_F  = []
    DIHTYPE_PHASE = []
    DIHTYPE_K = []
    DIHTYPE_PN = []
    DIHTYPE_C = []

    debug = 0 
    if(debug): print "   norm_dihparam = ",norm_dihparam

    # Find dihedral type parameters
    for ind in range( len( DTYPE_REF ) ):
        AT_i =  DTYPE_REF[ind][0]
        AT_j =  DTYPE_REF[ind][1]
        AT_k =  DTYPE_REF[ind][2]
        AT_l =  DTYPE_REF[ind][3]
        check = 1
        for ff_i in range (len(FF_DIHTYPES)):
            FF_l = FF_DIHTYPES[ff_i] #.split()

            if ( FF_l[0]  == AT_i and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  FF_l[3] == AT_l ):
                check = 0
                break
            if ( FF_l[0]  == AT_l and  AT_k == FF_l[1] and  AT_j == FF_l[2] and  FF_l[3] == AT_i ):
                check = 0
                break
	    
	if( check ):
	    for ff_i in range (len(FF_DIHTYPES)):
		FF_l = FF_DIHTYPES[ff_i] #.split()
		
		if ( FF_l[0]  == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  FF_l[3] == AT_l ):
		    check = 0
		    break
		if ( FF_l[0]  == AT_i and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  FF_l[3] == 'X' ):
		    check = 0
		    break
		if ( FF_l[0]  == 'X' and  AT_j == FF_l[1] and  AT_k == FF_l[2] and  FF_l[3] == 'X' ):
		    check = 0
		    break
		
		if ( FF_l[0]  == 'X' and  AT_k == FF_l[1] and  AT_j == FF_l[2] and  FF_l[3] == AT_i ):
		    check = 0
		    break
		if ( FF_l[0]  == AT_l and  AT_k == FF_l[1] and  AT_j == FF_l[2] and  FF_l[3] == 'X' ):
		    check = 0
		    break
		if ( FF_l[0]  == 'X' and  AT_k == FF_l[1] and  AT_j == FF_l[2] and  FF_l[3] == 'X' ):
		    check = 0
		    break
	    

        if check :
            print ' dih ',ind,' ',AT_i,' - ',AT_j,' - ',AT_k,' - ',AT_l,'  not in ff',itp_file 
            # print  FF_DIHTYP[ff_i]
	    print 'unknow ff specification '
            sys.exit(' top.dih_parameters ')
        else:
            
            dihen_norm = 1.0
	    
            if( norm_dihparam ):
		print " Normalizing dihedral potential "
                # normalize by number of nieghbors
                if(debug): print " finding types for ",AT_j,AT_k
                for ind in range( len(ATYPE_REF) ):
                    if( debug): print ATYPE_REF[ind]
                    AT_ref =  ATYPE_REF[ind].strip()
                    if(debug): print " checking ",AT_ref
                    if( AT_ref == AT_j ):
                        if(debug): print "  found j",AT_j,ind
                        NNAB_i = ATYPE_NNAB[ind] - 1
                    if( AT_ref == AT_k ):
                        if(debug): print "  found  k",AT_k,ind
                        NNAB_j = ATYPE_NNAB[ind] - 1
		
                dihen_norm = float( NNAB_i + NNAB_j) #/2.0
		
		if(debug): print " dihen_norm ",dihen_norm
                
            func_type = int(FF_l[4])
            DIHTYPE_F.append( func_type )
            if( func_type == 1 or func_type == 9 ):  # Harmonic
                DIHTYPE_PHASE.append( float(FF_l[5]))        # phi ( deg)
                DIHTYPE_K.append( float(FF_l[6]) *KJ_KCAL/dihen_norm  )   # KJ/mol to kcal/mol 
                DIHTYPE_PN.append( float(FF_l[7]) )           # multiplicity 
            elif( func_type == 2 ): # improper
                DIHTYPE_PHASE.append( float(FF_l[5]) )        # (deg)
                DIHTYPE_K.append( float(FF_l[6]) *KJ_KCAL/dihen_norm  )   # KJ/mol/rad to KJ/mol/rad 
            elif( func_type == 3 ): # Ryckert-Bellemans 
                C0 = float(FF_l[5])*KJ_KCAL/dihen_norm 
                C1 = float(FF_l[6])*KJ_KCAL/dihen_norm 
                C2 = float(FF_l[7])*KJ_KCAL/dihen_norm 
                C3 =  float(FF_l[8])*KJ_KCAL/dihen_norm 
                C4 =  float(FF_l[9])*KJ_KCAL/dihen_norm 
                F1 = -1.0*( 2.0*C1 + 3.0*C3/2.0)
                F2 = -1.0*( C2 + C4)
                F3 = -0.5*C3
                F4 = -0.25*C4
                DIHTYPE_C.append( [ F1,F2,F3,F4 ]  )        # (deg)
            else :
                print ' Dihedral type ', func_type , ' unknown '
		print 'unknow ff specification '
		sys.exit(' top.dih_parameters ')
	    
        if(debug): print ' dih ',ind,' ',func_type,AT_i,' - ',AT_j,' - ',AT_k,' - ',AT_l
    
    
    if(debug):
	sys.exit(' dih_parameters ')
    
    return ( DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C )


def imp_parameters(itp_file):
    import sys
    global GRO_SIG,KJ_KCAL,NM_ANG

    IMPTYPE_F = []

    return IMPTYPE_F

def pass_i(N_i,ELN_i,ASYMB_i,R_i,ATYPE_i,GTYPE_i,CHARGES_i,CHARN_i,AMASS_i,RESID_i,RESN_i,BONDS_i,ANGLES_i,DIH_i,options):
    import elements
    import numpy as np
    
    debug = 1


    ELN_j = []
    ASYMB_j = []
    R_j = []
    ATYPE_j = []
    GTYPE_j = []
    RESID_j = []
    RESN_j = []
    CHARN_j = []
    AMASS_j = []
    CHARGES_j = []
    REF_j = []
    GHOST_j = []

    
    BONDS_j = []
    ANGLES_j = []
    DIH_j = []
    IMPS_j = []
    
    FIX_j = []

    r_cov = elements.covalent_radi()
    
    # Set reference for atom i
    REF_i = []
    for atom_i in range( len(ELN_i) ):
	REF_i.append( atom_i )
    
    j_cnt = -1
    for atom_i in N_i:
	j_cnt += 1
	REF_i[ atom_i ] = j_cnt 
	REF_j.append( atom_i )
	ELN_j.append( ELN_i[atom_i] )
	ASYMB_j.append( ASYMB_i[atom_i] )
	R_j.append( R_i[atom_i] )
	ATYPE_j.append( ATYPE_i[atom_i] )
	GTYPE_j.append( GTYPE_i[atom_i] )
	CHARN_j.append( CHARN_i[atom_i] )
	AMASS_j.append( AMASS_i[atom_i] )
	CHARGES_j.append( CHARGES_i[atom_i] )
	RESID_j.append( RESID_i[atom_i] )
	RESN_j.append( RESN_i[atom_i] )
	GHOST_j.append( 0 )
	
	FIX_j.append( 0 )
	

    if( debug):
	print "  grp atoms ",j_cnt
	
    for b_indx in range( len(BONDS_i)):
        a_i = BONDS_i[b_indx][0] 
        a_j = BONDS_i[b_indx][1]
	
	i_check = 0
	j_check = 0
	for atom_i in N_i:
	    if( a_i == atom_i ):i_check = 1
	    if( a_j == atom_i ):j_check = 1
	    
	if( j_check and i_check ):
	    BONDS_j.append(  [REF_i[a_i]  ,REF_i[a_j]  ]   )
	    if( debug):
		print "    bonds ", REF_i[a_i]  ,REF_i[a_j]  
	    
	if( options.h_term ):
	    hterm_i = -1
	    if( i_check  and j_check == 0 ):
		hterm_i = a_j
		fix_i = a_i
	    if( j_check  and i_check == 0 ):
		hterm_i = a_i
		fix_i = a_j
	    if( hterm_i >=  0 and ASYMB_i[hterm_i] != 'VS' ):
	    
		a_fix = REF_i[fix_i]
		REF_i[fix_i] 

		r_bond =   R_i[hterm_i] - R_i[fix_i] 
		bond_ij = np.linalg.norm(r_bond)
		cr_i = r_cov[ ELN_i[fix_i] ]
		cr_j = r_cov[1]
		cr_ij = cr_i + cr_j

		# Scale factor to get correct covelent bond length 
		r_scale = cr_ij/bond_ij
		r_xh = r_scale*r_bond
		R_h =  R_i[fix_i] + r_xh
	    
		j_cnt += 1
		REF_i[ hterm_i ] = j_cnt
			
		REF_j.append( hterm_i )
		ELN_j.append( 1 )
		ASYMB_j.append( "H")
		R_j.append( R_h )
		ATYPE_j.append( "HC" )
		GTYPE_j.append( "H" )
		CHARN_j.append( CHARN_i[hterm_i] )
		AMASS_j.append( 1.008 )
		CHARGES_j.append(  0.06 )
		RESID_j.append( RESID_i[hterm_i] )
		RESN_j.append( RESN_i[hterm_i] )
		GHOST_j.append( 1 )
		
		FIX_j.append( 0 )
		FIX_j[a_fix] = 1
			
			
		BONDS_j.append( [a_fix  ,j_cnt  ] )
		
		# Add to list
		N_i.append( hterm_i)

		if( debug):
		    print "     atom ",ATYPE_i[hterm_i]," not in group "
		    print "     hterm i R(i)",hterm_i,R_i[hterm_i],GTYPE_i[hterm_i]
		    print "     bonded to i R(i)",fix_i,R_i[fix_i],GTYPE_i[fix_i]
		    print "     bond ",a_fix,j_cnt
		    print "     h-term bonds ", a_fix  ,j_cnt 

    for b_indx in range( len(ANGLES_i)):
        a_k = ANGLES_i[b_indx][0] 
        a_i = ANGLES_i[b_indx][1] 
        a_j = ANGLES_i[b_indx][2]
	
	k_check = 0
	i_check = 0
	j_check = 0
	for atom_i in N_i:
	    if( a_k == atom_i ):k_check = 1
	    if( a_i == atom_i ):i_check = 1
	    if( a_j == atom_i ):j_check = 1
	    
	if( k_check and j_check and i_check ):
	    ANGLES_j.append(  [REF_i[a_k]  ,REF_i[a_i]  ,REF_i[a_j]  ]   )
	    if( debug):
		print "    angles ",REF_i[a_k]  ,REF_i[a_i]  ,REF_i[a_j]  

    for b_indx in range( len(DIH_i)):
        a_k = DIH_i[b_indx][0]
        a_i = DIH_i[b_indx][1] 
        a_j = DIH_i[b_indx][2]
        a_l = DIH_i[b_indx][3]
	
	k_check = 0
	i_check = 0
	j_check = 0
	l_check = 0
	for atom_i in N_i:
	    if( a_k == atom_i ):k_check = 1
	    if( a_i == atom_i ):i_check = 1
	    if( a_j == atom_i ):j_check = 1
	    if( a_l == atom_i ):l_check = 1
	    
	if( k_check and i_check and j_check  and l_check ):
	    DIH_j.append(  [REF_i[a_k]  ,REF_i[a_i]  ,REF_i[a_j]  ,REF_i[a_l]  ]   )
	    if( debug):
		print "    dih ",REF_i[a_k]  ,REF_i[a_i]  ,REF_i[a_j]  ,REF_i[a_l]
		print "    dih index ",b_indx,DIH_i[b_indx]
		print "      dih types i ",ATYPE_i[a_k],ATYPE_i[a_i],ATYPE_i[a_j],ATYPE_i[a_l]
		print "      dih types j ",ATYPE_j[REF_i[a_k] ] ,ATYPE_j[REF_i[a_i] ] ,ATYPE_j[REF_i[a_j] ],ATYPE_j[REF_i[a_l] ]

    return (ELN_j,ASYMB_j,R_j,ATYPE_j,GTYPE_j,CHARGES_j,RESID_j,RESN_j,CHARN_j,AMASS_j,BONDS_j,ANGLES_j,DIH_j,IMPS_j,REF_j,GHOST_j)
	    
	    
def zero_unitq(ELN,ATYPE,CHARGES,tag2,NBINDEX,NBLIST,verbose,zero_term,zero_func):
    import numpy
    import sys 
    
    #
    # Sum exsessive charges into carbon atoms 
    #
    CHARGES_unit = numpy.zeros( len(ELN)  )			
    for atom_i in range( len(ELN) ):
	# Find each terminal hydrogen
	if( tag2[atom_i] == "X"  and zero_term  ):
	    term = atom_i 
	    # Check to be sure hydrogen
	    if( ELN[atom_i] != 1 ):
		print " Non hydrogen used as terminal group "
		sys.exit(" Code unable to process multi atom (nonhyrdogen) terminal group ")
	    if( verbose ):
		print " Terminal atom found ",atom_i+1," ",ATYPE[atom_i]
	    #
	    # Loop over nieghbors to find attached atom
	    #
	    N_o = NBINDEX[atom_i]
	    N_f = NBINDEX[atom_i+1] - 1
	    for j_indx in range( N_o,N_f+1):
		atom_j = NBLIST[j_indx]
		term_con_cnt = 0 # Terminal connections count
		if( tag2[atom_j].strip() == "T" ):
		    term_con_cnt += 1
		    if( verbose ):
			print " Terminal connection found ",atom_j+1," ",ATYPE[atom_j]
		    term_con = atom_j
	    # Check to be sure multiple atoms not found
	    if( term_con_cnt < 1 ):
		print " No terminal connections found "
		sys.exit(" Error in terminal connections ")
	    if( term_con_cnt > 1 ):
		print " Multiple terminal connections found "
		sys.exit(" Error in terminal connections ")
		
	    # Sum charges into base monomer unit
	    CHARGES[term_con] = CHARGES[term_con]  + CHARGES[term] 
	    CHARGES[term]  = 0.0
	    
	# Find each functional hydrogen
	if( tag2[atom_i] == "R"  and zero_func ):
	    term = atom_i 
	    # Check to be sure hydrogen
	    if( ELN[atom_i] != 1 ):
		print " Non hydrogen used as functional group "
		sys.exit(" Code unable to process multi atom (nonhyrdogen) functional group ")
	    if( verbose ):
		print " Functional atom found ",atom_i+1," ",ATYPE[atom_i]
	    #
	    # Loop over nieghbors to find attached atom
	    #
	    N_o = NBINDEX[atom_i]
	    N_f = NBINDEX[atom_i+1] - 1
	    for j_indx in range( N_o,N_f+1):
		atom_j = NBLIST[j_indx]
		term_con_cnt = 0 # Terminal connections count
		if( tag2[atom_j].strip() == "F" ):
		    term_con_cnt += 1
		    if( verbose ):
			print " functional connection found ",atom_j+1," ",ATYPE[atom_j]
		    term_con = atom_j
	    # Check to be sure multiple atoms not found
	    if( term_con_cnt  < 1 ):
		print " No functional connections found "
		sys.exit(" Error in functional connections ")
	    # Check to be sure multiple atoms not found
	    if( term_con_cnt > 1 ):
		print " Multiple functional connections found "
		sys.exit(" Error in functional connections ")
	    # Sum charges into base monomer unit
	    CHARGES[term_con] = CHARGES[term_con]  + CHARGES[term] 
	    CHARGES[term]  = 0.0
	    
    return CHARGES


def set_cply_tags(verbose, ELN, CTYPE,UNITNUMB ,NBLIST, NBINDEX ):
    #
    #  Set tags for new cply file 
    #     Use ctype tag and bonding enviroment
    #
    cply_tag = []
    for atom_i in range( len(ELN) ):
	cply_tag.append("")
    
    if( verbose): print "    setting cply tags "
    for atom_i in range( len(ELN) ):
	i_n = atom_i + 1

	print "CHECKING CTYPE ",atom_i,i_n, CTYPE[atom_i] ,  ELN[atom_i]
	
	# Set terminal attached carbons 
	if( CTYPE[atom_i] == "T" and ELN[atom_i] == 6 ):
	    #
	    term_con_cnt = 0 # Terminal connections count
	    #
	    # Find terminal hydrogen
	    #		    
	    N_o = NBINDEX[atom_i]
	    N_f = NBINDEX[atom_i+1] - 1
	    for j_indx in range( N_o,N_f+1):
		atom_j = NBLIST[j_indx]
		j_n = atom_j + 1
		if( CTYPE[atom_j] == "X" and ELN[atom_j] == 1 ):
		    # Check to see if it has been bonded to another unit 
		    if( UNITNUMB[atom_i] ==  UNITNUMB[atom_j] ):
			term_con_cnt += 1
			cply_tag[atom_j] = "termcap_H(" + str(j_n) + ")_on_C("+str(i_n)+")"
		      
	    if(   term_con_cnt == 1 ):
		cply_tag[atom_i] = "term_C(" + str(i_n) + ")"
	    if( term_con_cnt > 1 ):
		print " Number of terminal atoms attached to atom ",i_n," greater than 1 "
		sys.exit(" Error in terminal connections ")
	    
	# Set functional attached carbons 
	if( CTYPE[atom_i] == "F" and ELN[atom_i] == 6 ):
	    #
	    term_con_cnt = 0 # Terminal connections count
	    # 
	    # Find function hydrogen
	    # 
	    N_o = NBINDEX[atom_i]
	    N_f = NBINDEX[atom_i+1] - 1
	    for j_indx in range( N_o,N_f+1):
		atom_j = NBLIST[j_indx]
		j_n = atom_j + 1
		print " ",atom_j,j_n,CTYPE[atom_j],ELN[atom_j]
		if( CTYPE[atom_j] == "R" and ELN[atom_j] == 1 ):
		    # Check to see if it has been bonded to another unit 
		    if( UNITNUMB[atom_i] ==  UNITNUMB[atom_j] ):
			term_con_cnt += 1
			cply_tag[atom_j] = "funccap_H(" + str(j_n) + ")_on_C("+str(i_n)+")"

	    if( term_con_cnt == 1 ):
		cply_tag[atom_i] = "func_C(" + str(i_n) + ")"
		
	    if( term_con_cnt > 1 ):
		print " Number of functional atoms attached to atom ",i_n," not equal to 1 "
		sys.exit(" Error in functional connections ")
		
    return cply_tag


def print_ff_files(ff_prefix,verbose,ff_software,itp_file,ff_charges,norm_dihparam,ELN,R,CHARGES,LV):
    
    import elements , lammps ,gromacs , atom_types
    
    
    calc_id = ff_prefix
    
    
    # Read in ff file
    if( verbose ):
	print "   Read in parameters from ",itp_file
    FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(itp_file)
    
    
    if( verbose ):
	print "      Setting atomic data "
    ASYMB = elements.eln_asymb(ELN)
    AMASS = elements.eln_amass(ELN)
    RESID = initialize_resid( ELN )
    
    #   Initialize  topology values
    GTYPE = initialize_gtype( ELN )
    RESID = initialize_resid( ELN )
    RESN = initialize_resn( ELN )
    CHARN = initialize_charn( ELN )
    
    if( verbose ):
	print "      Setting bonded data ",itp_file
    
    #   Build covalent nieghbor list for bonded information 
    NBLIST, NBINDEX = build_covnablist(ELN,R)
    
    NA = len(ELN)
    BONDS = nblist_bonds(NA,NBLIST, NBINDEX)
    ANGLES = nblist_angles(NA,NBLIST, NBINDEX)
    DIH = nblist_dih(NA,NBLIST, NBINDEX)
    IMPS = nblist_imp(NA,NBLIST, NBINDEX,ELN)
			
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
	
    d_mass = 0
    d_charge = 0
    RINGLIST, RINGINDEX , RING_NUMB = find_rings(ELN,NBLIST,NBINDEX)
    
    # Asign oplsaa atom types
    ATYPE, CHARGES = atom_types.oplsaa( ff_charges,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB)
    
    ATYPE , CHARGES = atom_types.biaryl_types( ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
    #Refind inter ring types
    ATYPE , CHARGES  = atom_types.interring_types(ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
    
    
    # Identify total number of atom types for lammps output 
    ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH)
    
    # Check atom types to be sure each atom of the same type has the same number of neighbors 
    ATYPE_NNAB = check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)
    
    # Modify dih potentials
    ## FF_DIHTYPES = mod_dih( ATYPE_IND 
    
    # Print new itp file with only used atom types and interactions
    if( ff_software == "gromacs"):
	new_itp = "ff-test.itp"
	#AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(options,ASYMB,ATYPE,BONDS,ANGLES,DIH,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
	AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(new_itp,norm_dihparam,ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
    
	job_gro = calc_id +  ".gro"
	gromacs.print_gro(job_gro,GTYPE,RESID,RESN,R,LV)
    
    ATYPE_EP, ATYPE_SIG = atom_parameters(itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
    BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = bond_parameters(itp_file,BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
    ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = angle_parameters(itp_file,ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
    DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = dih_parameters(itp_file, norm_dihparam, DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )
    
    IMPTYPE_F  = imp_parameters(itp_file)
    
    if( ff_software == "gromacs"):
	#
	# Set charge groups
	#
	CHARN = set_chargegroups(verbose,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LV)
	#
	# Write top file
	#
	top_file =  calc_id +  ".top"
	const = []
	angle = []
	gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
		       ,const,angle
		       ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
		       ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
       
			
	
    elif( ff_software == "lammps" ):
	
	data_file = calc_id + ".data"
	
	lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
	      BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
	      ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
	      DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
	      RESN,ATYPE_IND,CHARGES,R , ATYPE,
	      BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV)
    
    
	