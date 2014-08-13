
def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    # Atom type options 
    parser.add_option("-i","--id_i", dest="id_i", type="string", default="", help="Atom types (FFtype,GROMACStype) of group i ")
    parser.add_option("-j","--id_j", dest="id_j", type="string", default="", help="Atom types (FFtype,GROMACStype) of group j ")


    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")

    parser.add_option("--rdf_out", dest="rdf_out", type="string", default="rdf.dat", help="Output rdf file ")
 
    parser.add_option("--r_cut", dest="r_cut", type=float, default=10.0, help=" Cut off radius in angstroms ")
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.10, help=" Bin size in angstroms")

    parser.add_option("--cubic",dest="cubic", default=False,action="store_true", help="Use cubic pbc's for speed up ")

    parser.add_option("--mol_inter",dest="mol_inter", default=False,action="store_true", help="Use only inter molecular rdf's ")
    parser.add_option("--mol_intra",dest="mol_intra", default=False,action="store_true", help="Use only intra molecular rdf's")
    
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=0, help=" Initial frame to read")
    
    parser.add_option("--frame_sufx", dest="frame_sufx", type="string", default=".gro", help=" sufix of frame file ")

    parser.add_option("--n_den", dest="n_den", type=float, default=1.0, help=" Reference number density N/A^3")
    parser.add_option("--nb_list",dest="nb_list", default=False,action="store_true", help="Use neighbor list ")

    
    (options, args) = parser.parse_args()
        
    return options, args
   
   
def build_nablist(options,ELN,R,LV):

    NNAB = -1
    na = len( ELN )
    sq_r_cut = options.r_cut**2
   
    # Calculate
    
    if( options.n_den < 0.1 ): options.n_den = 0.1
    
    den_buffer = 1.25 
    vol_cut = 4*numpy.pi/3*(r_cut**3)
    nbs_atom = vol_cut*options.n_den*den_buffer
	
    maxnnab = na*nbs_atom

    cov_buffer = 1.25 
    
    radi_cov =  []
    NBLIST = numpy.empty( maxnnab,  dtype=int )
    NBINDEX = numpy.empty( na+2,  dtype=int )

    for atom_i in range(na-1):
	
	NBINDEX[atom_i] = NNAB + 1
	
	add_i = 0 
	for id_indx in range( len(id_i)):
	    if( id_i[id_indx] == ATYPE[atom_i].strip() ): add_i = 1
	    if( id_i[id_indx] == GTYPE[atom_i].strip() ): add_i = 1
	if( add_i ):
		
	    r_i =  R[atom_i]
		    
	    for atom_j in range( atom_i+1, na ):
		add_j = 0 
		for id_indx in range( len(id_j)):
		    if( id_j[id_indx] == ATYPE[atom_j].strip()  ): add_j = 1
		    if( id_j[id_indx] == GTYPE[atom_j].strip()  ): add_j = 1
		if( add_j ):
		    if( atom_i != atom_j ):
			r_j =  R[atom_j] 
			sq_r_ij = prop.sq_drij_c(options,r_i,r_j,LV)
			if( sq_r_ij <= sq_r_cut ):
			    NNAB += 1
			    NBLIST[NNAB] =  atom_j
			    
    # Account for final atom, which has no connections 
    NBINDEX[atom_i+1] =  NNAB + 1
    NBINDEX[atom_i+2] =  NNAB + 1

    return NBLIST, NBINDEX 

def main():
    #
    # Caclulate rdf
    #
    
    import os, sys, numpy , math 
    import datetime
    import  gromacs, elements, xmol, prop, file_io, groups  #, vectors 
    
    
    
    options, args = get_options()
    
    
    p_time = 0
    
    if( options.verbose ): print "   Reading in files to establish intial conditions and atom types "
    #
    # Read in top file
    #
    if( len(options.in_top) ):
        if( options.verbose ): print "      Reading in ",options.in_top
        ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST = gromacs.read_top(options.in_top)
        ASYMB , ELN  = elements.mass_asymb(AMASS)
    #
    # Get coord
    #
    if( len(options.in_gro) ):
        if( options.verbose ): print "     Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options.in_gro)
#	
#    #
#    # Print initial properties 
#    #
#    if( options.verbose ):
#	print "  Initial volume of reference file ",vol_o," angstroms^3 "
#	print "  Cut off radius ",options.r_cut," angstroms"
#	print "  Bin size ",options.bin_size," angstroms"
#	print "  Number of bins ",n_bins
#        print "  Number density for neighbor list ",NUMB_DENSITY," atoms/angstroms^3 "
#	print "  Group i types "	    
#	for id_indx in range( len(id_list_i)):
#	    print "      ",id_list_i[id_indx]
#	print "  Group j types "	    
#	for id_indx in range( len(id_list_j)):
#	    print "      ",id_list_j[id_indx]
#	    
    #
    list_i = []
    #
    sum_i = 0
    #
    # Find atom indices  of group i and j
    #
    id_list_i = ["CHN"]
    
    print " Check types "
    for id_indx in range( len(id_list_i)):
	print id_list_i[id_indx]

    #sys.exit("cent test 0 ")
	
    for atom_i in range( len(ASYMB) ):
	
	add_i = 0 
	for id_indx in range( len(id_list_i)):
	    if( id_list_i[id_indx] == ATYPE[atom_i].strip() ): add_i = 1
	    if( id_list_i[id_indx] == GTYPE[atom_i].strip() ): add_i = 1
	    if( id_list_i[id_indx] == RESID[atom_i].strip() ): add_i = 1
	if( add_i ):
	    print  atom_i, ATYPE[atom_i], GTYPE[atom_i],  RESID[atom_i]
	    
	    list_i.append( atom_i )
	    sum_i += 1
	    
    print " atoms in list ",sum_i
    
    # Intialize center of mass
    prop_dim = 3
    total_mass = 0.0	    
    cent_mass = numpy.zeros(prop_dim)
    # delta_dr = numpy.zeros(prop_dim)
    
    for atom_i in list_i:
	a_mass = AMASS[atom_i]
	total_mass += a_mass
        for d in range(prop_dim):
            cent_mass[d] += a_mass*R[atom_i][d]
	    
            
	    
    # Normalize 
    for d in range(prop_dim):
        cent_mass[d] = cent_mass[d]/total_mass
    
    print " Intial cent of mass ",cent_mass
    
    
    dr_o = cent_mass
    
    R_cent = numpy.zeros( [len(ASYMB),prop_dim] )

    r_c = numpy.zeros(prop_dim)
    for atom_i in range( len(ASYMB) ):
        for d in range(prop_dim):
	    #R_cent.append( r_c )
	    r_i = R[atom_i][d]  - dr_o[d]
	    R_cent[atom_i][d]  = r_i - LV[d,d] * round( r_i/  LV[d,d] )
	
	if( atom_i < 4 ):
	    print atom_i, R[atom_i] ,  dr_o , r_c, R_cent[atom_i]
	
    out_xyz = "R_cent.xyz"
    file_xmol ="R_cent.xmol"
    xmol.write_xyz(ASYMB,R_cent,out_xyz)
    

    if( p_time ): t_i = datetime.datetime.now()
#    
#    #
#    # Create neighbor list
#    #
#    if( options.nb_list ):
#	print "  Creating neighbor list"
#	NBLIST, NBINDEX = build_nablist(options,ELN,R,LV)
    #
    #
    #
    #
    
    #
    for frame_i in range(options.frame_o,options.frame_f+1):
	frame_id = "n"+str(frame_i)+options.frame_sufx
	if( file_io.file_exists( frame_id ) ):
	    GTYPE, R_f, VEL, LV = gromacs.read_gro(frame_id)
	    
	    cent_mass_i = numpy.zeros(prop_dim)
	    # Find center of mass for frame i 
	    for atom_i in list_i:
		a_mass = AMASS[atom_i]
		for d in range(prop_dim):
		    cent_mass_i[d] += a_mass*R_f[atom_i][d]
		    
	    # Normalize 
	    for d in range(prop_dim):
		cent_mass_i[d]= cent_mass_i[d]/total_mass
		# Limit dr to reduce jumping
		delta_dr  =  cent_mass_i[d]-  dr_o[d]
		if( delta_dr > 1.0 ):
		    cent_mass_i[d] = dr_o[d] + 1.0 
		if( delta_dr < -1.0 ):
		    cent_mass_i[d] = dr_o[d] - 1.0 
	    
	    dr_o = numpy.zeros(prop_dim)
	    dr_o = cent_mass_i
	    
	    R_c_i = numpy.zeros( [len(ASYMB),prop_dim] )
	
	    for atom_i in range( len(ASYMB) ):
		for d in range(prop_dim):
		    #R_cent.append( r_c )
		    r_i = R_f[atom_i][d]  - cent_mass_i[d]
		    R_c_i[atom_i][d]  = r_i - LV[d,d] * round( r_i/  LV[d,d] )
		
	    xmol.print_xmol(ASYMB,R_c_i,file_xmol)
		
	    
	
	
if __name__=="__main__":
    main()
   
