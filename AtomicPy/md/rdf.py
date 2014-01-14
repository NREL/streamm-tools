
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
    
    debug = 0
    debug_2 = 0
    p_time = 1
    
    options, args = get_options()
    
    # Check options
    if( options.mol_inter and options.mol_intra ):
	print " Options --mol_inter and --mol_intra are mutually exclusive "
	sys.exit("Error is specified options ")
    
    if( options.verbose ): print "   Reading in files to establish intial conditions and atom types "
    #
    # Read in top file
    #
    if( len(options.in_top) ):
        if( options.verbose ): print "      Reading in ",options.in_top
        ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST = gromacs.read_top(options,options.in_top)
        ASYMB , ELN  = elements.mass_asymb(AMASS)
    #
    # Get coord
    #
    if( len(options.in_gro) ):
        if( options.verbose ): print "     Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options,options.in_gro)
    #
    # Calculate square of cut off
    #
    sq_r_cut = options.r_cut**2
    #
    # Calculate intial volume and total density 
    #
    vol_o = prop.volume( LV )
    NP = len( ELN )
    NUMB_DENSITY = float(NP)/vol_o
    options.n_den = NUMB_DENSITY
    #
    # Intialize lsits and counts
    #
    n_bins = int(options.r_cut/options.bin_size)
    #
    # Define rdf groups 
    #
    id_list_i  = options.id_i.split()
    id_list_j  = options.id_j.split()
    
    #
    # Print initial properties 
    #
    if( options.verbose ):
	print "  Initial volume of reference file ",vol_o," angstroms^3 "
	print "  Cut off radius ",options.r_cut," angstroms"
	print "  Bin size ",options.bin_size," angstroms"
	print "  Number of bins ",n_bins
        print "  Number density for neighbor list ",NUMB_DENSITY," atoms/angstroms^3 "
	print "  Group i types "	    
	for id_indx in range( len(id_list_i)):
	    print "      ",id_list_i[id_indx]
	print "  Group j types "	    
	for id_indx in range( len(id_list_j)):
	    print "      ",id_list_j[id_indx]
	    
    #
    list_i = []
    list_j = []
    #
    sum_i = 0
    sum_j = 0 
    #
    # Find atom indices  of group i and j
    #
    for atom_i in range( len(ASYMB) ):
	add_i = 0 
	for id_indx in range( len(id_list_i)):
	    if( id_list_i[id_indx] == ATYPE[atom_i].strip() ): add_i = 1
	    if( id_list_i[id_indx] == GTYPE[atom_i].strip() ): add_i = 1
	if( add_i ):
	    list_i.append( atom_i )
	    sum_i += 1 
	add_j = 0 
	for id_indx in range( len(id_list_j)):
	    if( id_list_j[id_indx] == ATYPE[atom_i].strip()  ): add_j = 1
	    if( id_list_j[id_indx] == GTYPE[atom_i].strip()  ): add_j = 1
	if( add_j ):
	    list_j.append( atom_i )
	    sum_j += 1
	    
	    
    if( p_time ): t_i = datetime.datetime.now()
    #
    # Create neighbor list
    #
    if( options.nb_list ):
	print "  Creating neighbor list"
	NBLIST, NBINDEX = build_nablist(options,ELN,R,LV)

    #
    # Loop over frames
    #
    rdf_cnt = numpy.zeros(n_bins+1)    
    frame_cnt = 0
    volume_i = [] #numpy.array()
    #
    for frame_i in range(options.frame_o,options.frame_f+1):
	frame_id = "n"+str(frame_i)+options.frame_sufx
	if( file_io.file_exists( frame_id ) ):
	    GTYPE, R_f, VEL, LV = gromacs.read_gro(options,frame_id)
	    volume_i.append(  prop.volume( LV ) )
	    frame_cnt += 1 
	    # xmol.print_xmol(ASYMB,R_i,file_xmol)
	    if( options.verbose ):
		print "reading ",frame_id
	    #R_frames.append( R_f )   
		
    	    #
	    # Loop over lists
	    #
	    #
	    for atom_i in list_i:
		r_i = R_f[atom_i]
		
		if( options.nb_list ):
		    list_j = []
		    N_o = NBINDEX[atom_i]
		    N_f = NBINDEX[atom_i+1] - 1
		    
		    print atom_i
		    print N_o,N_f
		    
		    for indx in range( N_o,N_f+1):
			atom_j = NBLIST[indx]
				    
			add_ij = 0
			for id_indx in range( len(id_j)):
			    if( id_j[id_indx] == ATYPE[atom_j].strip()  ): add_ij = 1
			    if( id_j[id_indx] == GTYPE[atom_j].strip()  ): add_ij = 1
			if( add_ij ):
			    list_j.append(atom_j)
		
		for atom_j in list_j:
		    add_ij = 1
		    if( atom_i >= atom_j ):
			add_ij = 0 
		    if( options.mol_inter ):
			if( MOLNUMB[atom_i] == MOLNUMB[atom_j] ): add_ij = 0 
		    if( options.mol_intra ):
			if( MOLNUMB[atom_i] != MOLNUMB[atom_j] ): add_ij = 0 
		    if( add_ij ):
			
			r_j = R_f[atom_j]
			
			if( options.cubic ):
			    sq_r_ij = prop.sq_drij_c(options,r_i,r_j,LV)
			else:
			    sq_r_ij = prop.sq_drij(options,r_i,r_j,LV)
			    
			if( sq_r_ij <= sq_r_cut ):
			    m_ij = numpy.sqrt(sq_r_ij)
			    bin_index = int( round( m_ij/options.bin_size) )
			    rdf_cnt[bin_index] += 2
			    
	
	else:
	    print " Error frame ",frame_id," does not exist "
	

    if( p_time ): 
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
	print "  rdf ti t_f ",t_i,t_f
	print "  rdf time ",dt_min,dt_sec
	
			
    if( options.verbose ):
	print "      Finding averages "
    #
    # Find averages
    #
    box_vol_ave = numpy.average( volume_i )
    vol_cut = 4.0*math.pi/3.0*options.r_cut**3
    n_shperes = float(sum_i)*float(frame_cnt)
    total_cnts = numpy.sum( rdf_cnt)
    sphere_den_j = float(total_cnts)/vol_cut/n_shperes #/2.0  # N_B A^-3
    box_den_i = float(sum_i )/float(box_vol_ave)
    box_den_j = float(sum_j )/float(box_vol_ave)
    
    if( options.verbose ):
	print "   N_i ",sum_i
	print "   N_j ",sum_j
	print "   Frames ",frame_cnt
	print "   Total counts ",total_cnts
	print "   Average box volume ",box_vol_ave
	print "   Volume of cut-off sphere ",vol_cut
	print "   Average box density i ",box_den_i," atoms/angstrom^3 "
	print "   Average box density j ",box_den_j," atoms/angstrom^3 "
	print "   Average cut-off sphere density ",sphere_den_j," atoms/angstrom^3 "
	
    # Write output 
    #
    rdf_file = options.rdf_out
    F_out = open(rdf_file,"w")
    F_out.write("# RDF frames %d %d " %  (options.frame_o,options.frame_f))
    F_out.write("\n#    Bin-size %f  " % (options.bin_size))
    F_out.write("\n#    Cut-off %f  " % (options.r_cut))
    F_out.write("\n#    Frames %d  " % (frame_cnt))
    F_out.write("\n#    Total_cnts %d  " % (total_cnts))
    F_out.write("\n#    N_i %d " % (sum_i ))
    F_out.write("\n#    N_j %d " % (sum_j ))
    F_out.write("\n#    Average Box Volume %f " % ( box_vol_ave) )
    F_out.write("\n#    Box density i %f N A^-3 " % (box_den_i ))
    F_out.write("\n#    Box density j %f N A^-3 " % (box_den_j ))
    F_out.write("\n#    Sphere volume  %f A^3 " % (vol_cut ))
    F_out.write("\n#    Average Sphere density  %f N A^3 " % (sphere_den_j ))
    F_out.write("\n#    ")
    F_out.write("\n# bin index ; r     ; count_ave/frame ; dr vol ;  dr vol(aprox) ; g_sphere ; g_boxs  ")
    #                bin_index , r_val , dr_cnt_norm      , dr_vol,  dr_vol_apx,     sphere_g, box_g
    
    for bin_index in range( 1,n_bins):
	    
	r_val = options.bin_size*float(bin_index)
	r_in = r_val - options.bin_size*0.5
	r_out = r_val + options.bin_size*0.5

	dr_cnt = float( rdf_cnt[bin_index] )
	dr_cnt_norm =     dr_cnt    /float(frame_cnt)
	dr_vol = 4.0*math.pi/3.0*( r_out**3 - r_in**3 )
	dr_vol_apx = 4.0*math.pi*(  r_val**2 )
	
	dr_rho = dr_cnt_norm/dr_vol
	sphere_g = dr_rho/sphere_den_j/float( sum_i )
	box_g = dr_rho/box_den_j/float( sum_i )
	
	F_out.write("\n  %d %f %f %f %f %f %f " % (bin_index,r_val,dr_cnt_norm,dr_vol,dr_vol_apx,sphere_g,box_g) )
	
    F_out.close()
	
	
if __name__=="__main__":
    main()
   
