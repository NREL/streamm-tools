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


def create_top(strucC,ff_charges): # Move out of class (or derived class)
    """
    Find topology information for force-field input files 
    """

    # New options that need to be passed 
    limdih =  0
    limitdih_n = 1
    verbose = True
    debug = True

    if( debug):
        print " Bonds ",len(strucC.bondC) 
        print " Angles ",len(strucC.angleC) 
        print " Dihedrals ",len(strucC.dihC) 

    cov_nblist, cov_nbindx = build_covnablist(strucC)
    # Check bonds
    if( len(strucC.bondC) <= 0 ):
        if( verbose ):
            print "    Finding bonds "
        nblist_bonds(strucC,cov_nblist, cov_nbindx)
    # Check angles
    if( len(strucC.angleC) <= 0 ):
        if( verbose ):
            print "    Finding angles "
        nblist_angles(strucC,cov_nblist, cov_nbindx)
    # Check dihedrals
    if( len(strucC.dihC) <= 0 ):
        if( verbose ):
            print "    Finding dihedrals  "
        nblist_dih(strucC,cov_nblist, cov_nbindx)
        

def nblist_bonds(strucC,cov_nblist, cov_nbindx):
    import sys
    """
    Generate bonds from neighbor list
    """
    for pid_i, ptclObj_i  in strucC.ptclC:    
        N_o = cov_nbindx[pid_i ]
        N_f = cov_nbindx[pid_i+ 1 ] - 1
        for indx_j in range( N_o,N_f+1):
            pid_j = cov_nblist[indx_j]
            if( pid_j > pid_i):
                b_i = Bond( pid_i, pid_j )            
                strucC.bondC.put(b_i)

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
                        a_i = Angles( pid_k,pid_i, pid_j )            
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
                                            self.dihC.put(d_i)
                                            dih_ij_cnt += 1
					    atom_k_i = pid_k
					    atom_l_i = pid_l
					elif( limitdih_n !=  2 ):
                                            d_i = Dihedral( pid_k, pid_i, pid_j, pid_l )            
                                            self.dihC.put(d_i)
					    dih_ij_cnt += 1
					
				else:
                                    d_i = Dihedral( pid_k, pid_i, pid_j, pid_l )            
                                    self.dihC.put(d_i)
    

def build_covnablist(strucC):
    """
    Build covalent neighbor list from elements and positions 
    """

    debug = False 
    p_time = True
    
    n_ptcl = len( strucC.ptclC )
    maxnnab = n_ptcl*12            # Assume max nieghbors as fcc 

    cov_buffer = 1.25
    
    radi_cov =  []
    cov_nblist = np.empty( maxnnab,  dtype=int )
    cov_nbindx = np.empty( maxnnab,  dtype=int )

    NNAB = 0

    if( p_time ): t_i = datetime.datetime.now()

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


def bonded_nblist(strucC):
    """
    Create neighbor list of bonded particles
    """
    sys.exit("bonded_nblist not working !!! ")
        
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
    for p_i in range( len(nblist_py)):
        # loop over  each particle p_i and get list of neighbors nlist_i
        nlist_i = nblist_py[p_i]
        NBINDEX.append( NNAB + 1 )
        # Loop over neighbor list of each particle nlist_i and get neighbor p_j
        for p_j in  nlist_i:

            if(debug): print " p_i ",p_i," p_j ",p_j

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
