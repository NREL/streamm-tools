#! /usr/bin/env python
"""
Subroutines for calculating structural properties
 length - angstroms
 mass   - AMU
 volume - angstroms^3
"""

# Dr. Travis Kemper
# Gatech
# Initial Date 2012
# travis.kemper@nrel.gov



const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na


def r_frac(r_i,LV):
    """
    Translate real coordinate to fractional
    """
    
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
    """
    Apply periodic boundry conditions to fractional coordinates 
    """
    prop_dim = 3

    f_c_pbc  = []
    
    for d in range(prop_dim):
	f_comp = f_i[d]
	f_c_pbc.append( f_comp - round(f_comp,0) )
	
    return f_c_pbc

def mag_drij(options,r_i,r_j,LV):
    """
    Return magnitude of distance between to vectors 
    """
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
    """
    Return square of magnitude of distance between to vectors 
    """
    import numpy 
    
    f_i = r_frac(r_i,LV)
    f_j = r_frac(r_j,LV)
    f_ij = f_j - f_i

    df_pbc = f_pbc(f_ij,options)
    dr_pbc = f_real(df_pbc,LV)
    
    #print "mag_drij"
    #print dr_pbc 
    
    sq_dr = numpy.dot( dr_pbc,dr_pbc)
    
    return sq_dr
    
    
def mag_drij_c(r_i,r_j,LV):
    """
    Reutrn magnitude of distance between to vectors  using cubic periodic boundry conditions 
    """
    
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
    
def sq_drij_c(r_i,r_j,LV):
    """
    Reutrn square magnitude of distance between to vectors  using cubic periodic boundry conditions 
    """
    
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
    """
    Apply cubic periodic boundry conditions to a vector 
    """
    
    
    print dr_x, LV
    
    dr_x_pbc = dr_x - LV[0,0] * round( dr_x/  LV[0,0] )
    dr_y_pbc = dr_y - LV[1,1] * round( dr_y/  LV[1,1] )
    dr_z_pbc = dr_z - LV[2,2] * round( dr_z/  LV[2,2] )
    
    return ( dr_x_pbc,dr_y_pbc,dr_z_pbc )


def f_real( f_i,LV):
    """
    Return real position of fractional coordinates
    """

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
    """
    Calculate volume
    """
    
    import numpy
    
    latvec_a = LV[0]
    latvec_b = LV[1]
    latvec_c = LV[2]
    
    br1 = numpy.cross(latvec_a,latvec_b)
    vol = numpy.dot(br1,latvec_c)

    return vol


def total_mass(  AMASS ):
    """
    Calculate total mass 
    """
    
    import numpy
    
    # Sum mass, charges
    total_mass = 0.0
    for atom_i in range( len(AMASS) ):
	total_mass += AMASS[atom_i]
	
    return total_mass

def vecdensity(  AMASS,LV ):
    """
    Calculate density 
    """
    
    import numpy
    
    volume_i = volume( LV )
    
    # Sum mass, charges
    sys_mass = total_mass(  AMASS )

    density_i = sys_mass/volume_i/const_avo*10.0
	
    return density_i



def build_nablist(ELN,R):
    """
    Build bonded nieghbor list from covalent radi 
    """
    
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


def shift_cent_mass(AMASS,R,r_shift):
    """
    Shift center of mass to a point 
    """

    import elements, numpy
    
    
    prop_dim = len(R[0])
       
    # Intialize center of mass
    total_mass = 0.0
    #r_mass.append(r_zero)
    r_mass = numpy.array( [0.0,0.0,0.0] )
    
    for atom_i in range( len(AMASS) ):
            
            a_mass = AMASS[atom_i]
            # sum center of mass
            total_mass = total_mass + a_mass
	    
            for d in range(prop_dim):
                r_mass[d] += a_mass*R[atom_i][d]
                
    # Normalize 
    for d in range(prop_dim):
        r_mass[d] = r_mass[d]/total_mass
        
    R_shift = []
    
      
    for atom_i in range( len(AMASS) ):
        R_shift.append(  R[atom_i]  -   r_shift - r_mass  )
	
	
    return R_shift



def sys_cent_mass(ASYMB,R,options):
    """
    Calculate system center of mass 
    """

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


def list_cent_mass(atom_list,R,options):
    """
    Calculate center of mass of a list of atoms 
    """

    import numpy 

    prop_dim = len(R[0])
              
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


def cent_mass(AMASS,R):
    """
    Calculate center of mass 
    """

    import numpy 

    prop_dim = len(R[0])
         
    # Intialize center of mass
    total_mass = 0.0
    #r_mass.append(r_zero)
    r_mass = numpy.array( [0.0,0.0,0.0] )
    
    for atom_i in range( len(AMASS) ):
            
            a_mass = AMASS[atom_i]
            # sum center of mass
            total_mass = total_mass + a_mass

            for d in range(prop_dim):
                r_mass[d] += a_mass*R[atom_i][d]
                
    # Normalize 
    for d in range(prop_dim):
        r_mass[d] = r_mass[d]/total_mass
        
    return r_mass

def shift_r(r_cm,R_local,options):
    """
    Shift coordinates 
    """

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
    """
    Apply pbc's single postion vector r_i
    """
    
    f_i = r_frac(r_i,LV)
    f_c_pbc = f_pbc(f_i,options)
    r_i_pbc = f_real( f_c_pbc,LV)

    return r_i_pbc
    
def atom_pbc(R,LV,options):
    import numpy
    """
    Apply pbc's to all atoms with center = [0,0,0]
    """
    
    R_pbc = numpy.zeros( [len(R),options.prop_dim] )
    
    for atom_i in range( len(R) ):
	r_i = R[atom_i]
	r_i_pbc = r_pbc(r_i,LV,options)
	R_pbc[atom_i] =  r_i_pbc
	
    
    return R_pbc



def hterm_seg(ELN,ASYMB,R,ATYPE,CHARGES,NBLIST,NBINDEX,ASYMB_l,REFNUMB_l,FIX_l,calc_id,options):
    """
    hydrogen terminate segment from a larger system for qm calculation
    """
    
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
    """
    Calculate lowest spin multiplicity 
    """
    e = sum( ELN )
    e_total = e + Q
    e_unpaired =  e_total % 2
    if( e_unpaired == 0 ):
        M = 1
    else:
        M = 2
        
    return M


def getDihedral(a,b,c,d):
    """
    Calculate dihedral angle 
    """
    import numpy 
    
    #  k - i - j -l 
    # a = r_k
    # b = r_i
    # c = r_j
    # d = r_l 
    
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
    
    
    angle_i = getAngle(v1v2,v2v3)
    #
    # Find sign of angle 
    #
    v1v3 = numpy.cross(v1,v3)
    sign_v = numpy.dot(v2,v1v3)
    
    if( sign_v < 0.0  ):
	angle_i = -1.0*angle_i
    
    return angle_i

def getNormedVector(a,b):
    """
    Normailze vector 
    """
    import numpy 

    delta_ba = b-a
    return (delta_ba)/numpy.linalg.norm(delta_ba)

def getAngle(a,b):
    """
    Calcuate angle
      k - i - j 
     a = r_ik
     b = r_ij
     cos( \theta ) = ( a dot b ) / ( |a| |b| )
    """
    #
    import numpy 

    a_norm = a/numpy.linalg.norm(a)
    b_norm = b/numpy.linalg.norm(b) 
    dot_ab = numpy.dot(a_norm,b_norm)
    
    if( dot_ab >= 1.0 ):
       ang_deg = 0.0
    elif(  dot_ab <= -1.0 ):
       ang_deg = 180.0
    else:    
        cos_ang = numpy.arccos(dot_ab )
        ang_deg = numpy.rad2deg( cos_ang )
    
    
    return ang_deg



def d_en_2p(en_1,en_3,h):
    """
    Numerical derivative of 2 points 
    """
    # en_1  eV
    # en_2  eV 
    
    EVKC=23.0605
    
    d_en = (en_3 - en_1)/float(h)/2.0  # eV/deg
    # change to meV
    #r_delta_en = int( 1000*d_en )      # meV/deg
    r_delta_en = 1000*d_en 	    # meV/deg
    
    d_en_p = 0
    if( r_delta_en > 0 ): d_en_p = 1
    if( r_delta_en < 0 ): d_en_p =  -1
	
    return d_en_p

def d2_en_3p(en_1,en_2,en_3,h):
    """
    Numerical derivative of 3 points 
    """
    # en_1  eV
    # en_2  eV 
    # en_3  eV 
    
    d2_en = ( en_3 - 2.0*en_2  + en_1 )/ float(h)**2
    
    return d2_en

def calc_d2(  numeric_func, h ):
    import numpy, sys
    
    """
    Find max and mins based on the second derivative
    """
    
    debug = 0
    
    EVKC=23.0605
    
    success = 0 
    
    min_val =  1e16 
    qm_max = -1e16
    ang_min = "nan"
    ang_max = "nan"
#    
#    print " reading ",dih_qm
#    
#    f_qm = open(dih_qm ,"r")
#    f_qm_Lines = f_qm.readlines()
#    f_qm.close()
#    
#    # Put energies in arrays 
#    tor_en = [] #numpy.zeros()
#    tor_angle = [] #numpy.zeros()
#    
#    for f_qm_line in f_qm_Lines:
#	f_qm_col = f_qm_line.split()
#    
#	if( len(f_qm_col) >= 3 and f_qm_col[0] != "#" ):
#	    tor_angle.append( float( f_qm_col[1] ) )
#	    tor_en.append( float( f_qm_col[2] ) )
#	    
#	    success = 1
#	    
    min_val = min( numeric_func)
    
    if( debug): print "  Mim ",min_val
    
    # Find inversion points in energy
    min_indx = []
    max_indx = []
    trans_list = []
    trans_indxs = []
    k_list = []
    
    # Test first point
    calc_minmaxtrans = 0 
    calc_maxmintrans = 0
    print " length ",len(numeric_func)  
    # loop over subsequent points 
    for indx_i in range(len(numeric_func) ):
	
	
	indx_m = indx_i - 1  # _m minus 
	indx_m_m = indx_m - 1 
	indx_p = indx_i + 1  # _p plus 
	indx_p_p = indx_p + 1 
	
	# apply boundry conditions 
	if( indx_m < 0 ): indx_m  =  -1*indx_m 
	if( indx_m_m < 0 ): indx_m_m  = -1*indx_m_m 
	 	
	if(debug): print "  m i p ",indx_m_m,indx_m,indx_i,indx_p,indx_p_p

	if( indx_p > len(numeric_func)  -1  ): indx_p  = indx_p - ( indx_p -  len(numeric_func) + 1)*2  # as 0.0 == 180.0 
	if( indx_p_p > len(numeric_func)  -1  ): indx_p_p  = indx_p_p -  ( indx_p_p -  len(numeric_func) + 1)*2  #  as 0.0 == 180.0 
	
	d_en_m =d_en_2p(numeric_func[indx_m] ,numeric_func[indx_i],h) 
	d_en_m_m =d_en_2p(numeric_func[indx_m_m] ,numeric_func[indx_m],h) 
	d_en_p =d_en_2p(numeric_func[indx_i] ,numeric_func[indx_p],h)
	d_en_p_p =d_en_2p(numeric_func[indx_p] ,numeric_func[indx_p_p],h)
	
	if(debug):
	    print "  m i p ",indx_m_m,indx_m,indx_i,indx_p,indx_p_p
	    print "  m i p ",numeric_func[indx_m_m] - min_val,numeric_func[indx_m] - min_val,numeric_func[indx_i] - min_val,numeric_func[indx_p] - min_val,numeric_func[indx_p_p] - min_val
	    print "     dm dp ",d_en_m_m,d_en_m,d_en_p,d_en_p_p
	    #print "     dm dp ",d_en_m,d_en_p
	
	
	if( d_en_m_m < 0 and d_en_m < 0  and d_en_p > 0  and d_en_p_p > 0 ):
	#if(  d_en_m < 0  and d_en_p > 0  ):
	    min_indx.append( indx_i )
	    min_en_i = numeric_func[indx_i]
	    calc_minmaxtrans = 1
	    min_transindx = indx_i
	    
	    if(debug): print " min found ",numeric_func[indx_i] - min_val
	    
	    if( calc_maxmintrans ):
		trans_ev = max_en_i - min_en_i
		trans_list.append(  trans_ev )
		
		calc_maxmintrans = 0 
		trans_indxs.append( [max_transindx,indx_i] )
			
		if(debug): print "   Found  max min trans ",max_en_i - min_val," -> ",min_en_i - min_val," trans = ",trans_ev
			
	    d2_en_i = d2_en_3p(numeric_func[indx_m] ,numeric_func[indx_i],numeric_func[indx_p],h)
	    k_list.append(  d2_en_i ) 
	    if(debug): print " min found ",numeric_func[indx_i], d2_en_i
	    print " min found ",numeric_func[indx_i], d2_en_i
	    
	    success = 1
	    
	    
	if( d_en_m_m > 0 and d_en_m > 0  and d_en_p < 0  and d_en_p_p < 0 ):
	#if( d_en_m > 0  and d_en_p < 0   ):
	    max_indx.append( indx_i )
	    max_en_i = numeric_func[indx_i]
	    calc_maxmintrans = 1
	    max_transindx = indx_i

	    if(debug): print " max found ",numeric_func[indx_i] - min_val
	    
	    if( calc_minmaxtrans ):
		trans_ev = max_en_i - min_en_i
		trans_list.append(  trans_ev )
		calc_minmaxtrans = 0 
		trans_indxs.append( [indx_i,min_transindx] )
		
			
		if(debug): print "   Found min max trans ",min_en_i - min_val," -> ",max_en_i - min_val," trans = ",trans_ev
		
	    success = 1
			

	    
    for inv_indx in range( len(min_indx) ):
	indx = min_indx[inv_indx]
	if(debug): print "  Min ",inv_indx," found at ",indx," w energy ",numeric_func[indx]
        
    for inv_indx in range( len(max_indx) ):
	indx = max_indx[inv_indx]
	if(debug): print "  Max ",inv_indx," found at ",indx," w energy ",numeric_func[indx]
    
    for inv_indx in range( len(trans_list) ):
	bar = trans_list[inv_indx]
	print "  Transion ",inv_indx," =  ",bar*EVKC," kcal/mol"
        
	
    if(debug):    sys.exit(" inversion testing ")
	    
    return (success,  numeric_func,min_indx,max_indx,trans_list,trans_indxs,k_list )

