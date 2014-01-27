#! /usr/bin/env python
# Subroutines for calculating structural properties

# Dr. Travis Kemper
# Gatech
# Initial Date 2012
# travis.kemper@nrel.gov

# length - angstroms
# mass   - AMU
# volume - angstroms^3


const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na


def r_frac(r_i,LV):
    # Translate real coordinate to fractional
    
    import numpy 
    # Return fractional coordinate of real position r

    f_i = numpy.array([0.0,0.0,0.0])

    #     Set lattice vectors
    AA = LV[0][0]
    AB = LV[0][1]
    AC = LV[0][2]
    BA = LV[1][0]
    BB = LV[1][1]
    BC = LV[1][2]
    CA = LV[2][0]
    CB = LV[2][1]
    CC = LV[2][2]

    X=r_i[0]
    Y=r_i[1]
    Z=r_i[2]

    f_i[0] =  -((-BC*CB*X+BB*CC*X+BC*CA*Y-BA*CC*Y-BB*CA*Z+BA*CB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)) 
    f_i[1] =  -((AC*CB*X-AB*CC*X-AC*CA*Y+AA*CC*Y+AB*CA*Z-AA*CB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))  
    f_i[2] =  -((-AC*BB*X+AB*BC*X+AC*BA*Y-AA*BC*Y-AB*BA*Z+AA*BB*Z)/(AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)) 
    
    return f_i

def f_pbc(f_i,options):
    pass

    prop_dim = 3

    f_c_pbc  = []
    
    for d in range(prop_dim):
	f_comp = f_i[d]
	f_c_pbc.append( f_comp - round(f_comp,0) )
	
    return f_c_pbc

def mag_drij(options,r_i,r_j,LV):
    import numpy 
    # Find magnitude of dr_ij
    
    f_i = r_frac(r_i,LV)
    f_j = r_frac(r_j,LV)
    f_ij = f_j - f_i

    df_pbc = f_pbc(f_ij,options)
    dr_pbc = f_real(df_pbc,LV)
    
    #print "mag_drij"
    #print dr_pbc 
    
    mag_dr =  numpy.linalg.norm(dr_pbc)
    
    return mag_dr
    
def sq_drij(options,r_i,r_j,LV):
    import numpy 
    # Find magnitude of dr_ij
    
    f_i = r_frac(r_i,LV)
    f_j = r_frac(r_j,LV)
    f_ij = f_j - f_i

    df_pbc = f_pbc(f_ij,options)
    dr_pbc = f_real(df_pbc,LV)
    
    #print "mag_drij"
    #print dr_pbc 
    
    sq_dr = numpy.dot( dr_pbc,dr_pbc)
    
    return sq_dr
    
    
def mag_drij_c(options,r_i,r_j,LV):
    import numpy 
    # Find magnitude of dr_ij
    
    
    r_ij = r_j - r_i
    
    r_x = r_ij[0] - LV[0,0] * round( r_ij[0]/  LV[0,0] )
    r_y = r_ij[1] - LV[1,1] * round( r_ij[1]/  LV[1,1] )
    r_z = r_ij[2] - LV[2,2] * round( r_ij[2]/  LV[2,2] )
    
    dr_pbc = numpy.array( [r_x,r_y,r_z] )
    
    #print "mag_drij_c"
    #print int( r_ij[0]/  LV[0,0] )
    #print int( r_ij[1]/  LV[1,1] )
    #print int( r_ij[2]/  LV[2,2] )
    
    #print dr_pbc 
    
    mag_dr =  numpy.linalg.norm(dr_pbc)
    
    return mag_dr
    
def sq_drij_c(options,r_i,r_j,LV):
    import numpy 
    # Find magnitude of dr_ij
    
    r_ij = r_j - r_i
    
    #r_x = numpy.absolute( r_j[0] - r_i[0] )
    #r_y = numpy.absolute( r_j[1] - r_i[1] )
    #r_z = numpy.absolute( r_j[2] - r_i[2] )
    
    r_x = r_ij[0] - LV[0,0] * round( r_ij[0]/  LV[0,0] )
    r_y = r_ij[1] - LV[1,1] * round( r_ij[1]/  LV[1,1] )
    r_z = r_ij[2] - LV[2,2] * round( r_ij[2]/  LV[2,2] )
    
    #r_x = r_x - LV[0,0] * round( r_x/  LV[0,0] )
    #r_y = r_y - LV[1,1] * round( r_y/  LV[1,1] )
    #r_z = r_z - LV[2,2] * round( r_z/  LV[2,2] )
    
    dr_pbc = numpy.array( [r_x,r_y,r_z] )
    
    sq_dr = numpy.dot( dr_pbc,dr_pbc)
    
    return sq_dr
    

def pbc_r_c(dr_x,dr_y,dr_z,LV):
    
    print dr_x, LV
    
    dr_x_pbc = dr_x - LV[0,0] * round( dr_x/  LV[0,0] )
    dr_y_pbc = dr_y - LV[1,1] * round( dr_y/  LV[1,1] )
    dr_z_pbc = dr_z - LV[2,2] * round( dr_z/  LV[2,2] )
    
    return ( dr_x_pbc,dr_y_pbc,dr_z_pbc )


def f_real( f_i,LV):
    # Return real position of fractional coordinates 

    r_i = []

    #     Set lattice vectors
    AA = LV[0][0]
    AB = LV[0][1]
    AC = LV[0][2]
    BA = LV[1][0]
    BB = LV[1][1]
    BC = LV[1][2]
    CA = LV[2][0]
    CB = LV[2][1]
    CC = LV[2][2]
    
    r_i.append( f_i[0]*AA + BA*f_i[1]+CA*f_i[2] )
    r_i.append( f_i[0]*AB + BB*f_i[1]+CB*f_i[2] )
    r_i.append( f_i[0]*AC + BC*f_i[1]+CC*f_i[2] )

    return r_i

def volume( LV ):
    # Calculate volume
    import numpy
    
    latvec_a = LV[0]
    latvec_b = LV[1]
    latvec_c = LV[2]
    
    br1 = numpy.cross(latvec_a,latvec_b)
    vol = numpy.dot(br1,latvec_c)

    return vol


def total_mass(  AMASS ):
    import numpy
    
    # Sum mass, charges
    total_mass = 0.0
    for atom_i in range( len(AMASS) ):
	total_mass += AMASS[atom_i]
	
    return total_mass

def vecdensity(  AMASS,LV ):
    import numpy
    
    volume_i = volume( LV )
    
    # Sum mass, charges
    total_mass = total_mass(  AMASS )

    density_i = total_mass/volume_i/const_avo*10.0
	
    return density_i



def build_nablist(ELN,R):
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
    
    # loop over atoms 
    for atom_i in range(na-1):
        NBINDEX[atom_i] = NNAB + 1
        el_i = ELN[atom_i]
        r_i = numpy.array( R[atom_i] )
	rc_i = radi_cov[el_i]*cov_buffer
        for atom_j in range( atom_i+1, na ):
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
                        print ' atom i/j ', atom_i,atom_j,el_i,el_j
                        print ' cov radi ', radi_cov[el_i] , radi_cov[el_j]
                        print '   r_ij ',r_ij
                        print '   r_cov ',r_cov
                    
    if( p_time ):
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	print "  build_nablist dt ",dt_min,dt_sec

    if(debug): sys.exit('debug')
    
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

        sys.exit('debug')


    return (NBLIST, NBINDEX)

def sys_cent_mass(ASYMB,R,options):
    import elements, numpy 
       
    ELN = elements.asymb_eln(ASYMB)
    AMASS = elements.eln_amass(ELN)
 
    # Intialize center of mass
    total_mass = 0.0
    #r_mass.append(r_zero)
    r_mass = numpy.array( [0.0,0.0,0.0] )
    
    for atom_i in range( len(ASYMB) ):
            
            a_mass = AMASS[atom_i]
            # sum center of mass
            total_mass = total_mass + a_mass

            for d in range(options.prop_dim):
                r_mass[d] += a_mass*R[atom_i][d]
                
    # Normalize 
    for d in range(options.prop_dim):
        r_mass[d] = r_mass[d]/total_mass
        
    return r_mass


def shift_r(r_cm,R_local,options):
    # shift system
    import numpy 
    
    debug = 1
    
    R_cent = numpy.zeros( [len(R_local),options.prop_dim] )
    	
    for atom_i in range(len(R_local)):
	for d in range(options.prop_dim):
	    r_o = R_local[atom_i][d]
	    r_i =  r_o + r_cm[d]
	    R_cent[atom_i][d] = r_i
	    
    return R_cent 


def r_pbc(r_i,LV,options):
    # Apply pbc's single postion vector r_i 
    
    f_i = r_frac(r_i,LV)
    f_c_pbc = f_pbc(f_i,options)
    r_i_pbc = f_real( f_c_pbc,LV)

    return r_i_pbc
    
def atom_pbc(R,LV,options):
    import numpy 
    # Apply pbc's to all atoms with center = [0,0,0]
    
    R_pbc = numpy.zeros( [len(R),options.prop_dim] )
    
    for atom_i in range( len(R) ):
	r_i = R[atom_i]
	r_i_pbc = r_pbc(r_i,LV,options)
	R_pbc[atom_i] =  r_i_pbc
	
    
    return R_pbc



def hterm_seg(ELN,ASYMB,R,ATYPE,CHARGES,NBLIST,NBINDEX,ASYMB_l,REFNUMB_l,FIX_l,calc_id,options):
    # hydrogen terminate segment from a larger system for qm calculation
    import numpy as np 
    import top, elements 
    
    debug = 0
    
    r_cov = elements.covalent_radi()
    
    R_local = R
    
    # List of terminal hydrogens
    ELN_hterm = [] 
    ASYMB_hterm = []
    R_hterm = []
    ATYPE_hterm = []
    CHARGES_hterm = []
    REFNUMB_hterm = []
    FIX_hterm = []
    GHOST_hterm = []

    # append outputed files to reference file 
    f_ref = file(options.ref_file, "a")
            
    # Loop over all atoms in segment
    for si_indx in range( len(ASYMB_l) ):
	atom_si = REFNUMB_l[si_indx]
        N_o = NBINDEX[ atom_si ]
        N_f = NBINDEX[ atom_si + 1 ] - 1
	# Loop over neighbers of atom i
	if( debug): print " Check segment ",ASYMB_l[si_indx]
	for indx_i in range(N_o,N_f+1):
            atom_j = NBLIST[indx_i]
	    # find if neighbor in segment
	    nb_nonseg = 1
	    if( debug): print "   for nieghbor ",ATYPE[atom_j]
	    for sj_indx in range( len(ASYMB_l) ):
		atom_sj = REFNUMB_l[sj_indx]
		if( atom_sj == atom_j ):
		    nb_nonseg = 0 
	    if( nb_nonseg ):
		if(debug): print " Non segment nieghbor found",ATYPE[atom_j]
		if( ELN[atom_j] == 1 ):
		    ELN_hterm.append( ELN[atom_j] )
		    ASYMB_hterm.append( ASYMB[atom_j] )
		    R_hterm.append( R_local[atom_j] )
		    ATYPE_hterm.append( ATYPE[atom_j] )
		    CHARGES_hterm.append( CHARGES[atom_j] )
		    REFNUMB_hterm.append( atom_j )
		    FIX_hterm.append( 0 )
		    GHOST_hterm.append( 1 )

		    # Fix attached atom 
		    FIX_l[si_indx] = 1 

		    if(debug): print "      adding H ",ATYPE[atom_j]
		    # Record added atom 
		    f_ref.write( " \n %s %d " % (calc_id,atom_j))

		if( ELN[atom_j] == 8 ):
		    ASYMB_hterm.append("H")
		    
		    r_bond = R_local[atom_j] - R_local[atom_si]
		    bond_ij = np.linalg.norm(r_bond)
		    cr_i = r_cov[ ELN[atom_si] ]
		    cr_j = r_cov[1]
		    cr_ij = cr_i + cr_j

		    # Scale factor to get correct covelent bond length 
		    r_scale = cr_ij/bond_ij
		    r_xh = r_scale*r_bond
		    R_h = R_local[atom_si] + r_xh

		    ELN_hterm.append( 1 )
		    R_hterm.append( R_h )
		    ATYPE_hterm.append( "HC" )
		    CHARGES_hterm.append( 0.06 )
		    REFNUMB_hterm.append( atom_j )
		    FIX_hterm.append( 0 )
		    GHOST_hterm.append( 1 )

		    # Fix attached atom 
		    FIX_l[si_indx] = 1 

		    if(debug): print "      adding H "
		    f_ref.write( " \n %s %d " % (calc_id,atom_j))

    f_ref.close()
    
    ELECTRONS_hterm = len(ASYMB_hterm )
    
		    
    return ( ASYMB_hterm , R_hterm , ATYPE_hterm , CHARGES_hterm, REFNUMB_hterm,ELECTRONS_hterm ,FIX_hterm ,GHOST_hterm ,FIX_l  )
    
    
def multiplicity(ELN,Q):
    
    e = sum( ELN )
    e_total = e + Q
    e_unpaired =  e_total % 2
    if( e_unpaired == 0 ):
        M = 1
    else:
        M = 2
        
    return M


def getDihedral(a,b,c,d):
    import numpy 
    
    debug = 0
    if(debug):
        print " a ",a 
        print " b ",b
        print " c ",c
        print " d ",d
        
    v1 = getNormedVector(a, b)
    v2 = getNormedVector(b, c)
    v3 = getNormedVector(c, d)
    v1v2 = numpy.cross(v1,v2)
    v2v3 = numpy.cross(v2,v3)
    return getAngle(v1v2,v2v3)

def getNormedVector(a,b):
    import numpy 

    delta_ba = b-a
    return (delta_ba)/numpy.linalg.norm(delta_ba)

def getAngle(a,b):
    import numpy 

    a_norm = a/numpy.linalg.norm(a)
    b_norm = b/numpy.linalg.norm(b) 
    dot_ab = numpy.dot(a_norm,b_norm)
    cos_ang = numpy.arccos(dot_ab )
    ang_deg = numpy.rad2deg( cos_ang )
    return ang_deg
