
def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    # Atom type options 
    parser.add_option("-k","--id_k", dest="id_k", type="string", default="", help="Atom types (FFtype,GROMACStype) of group k ")
    parser.add_option("-i","--id_i", dest="id_i", type="string", default="", help="Atom types (FFtype,GROMACStype) of group i ")
    parser.add_option("-j","--id_j", dest="id_j", type="string", default="", help="Atom types (FFtype,GROMACStype) of group j ")
    parser.add_option("-l","--id_l", dest="id_l", type="string", default="", help="Atom types (FFtype,GROMACStype) of group l ")


    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")

    parser.add_option("--hist_out", dest="hist_out", type="string", default="dih.hist", help="Output historgam file ")
 
 
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.10, help=" Bin size in degrees ")
    
    parser.add_option("--cubic",dest="cubic", default=False,action="store_true", help="Use cubic pbc's for speed up ")

    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=0, help=" Initial frame to read")
    
    parser.add_option("--frame_sufx", dest="frame_sufx", type="string", default=".gro", help=" sufix of frame file ")

    parser.add_option("--n_den", dest="n_den", type=float, default=1.0, help=" Reference number density N/A^3")
    parser.add_option("--nb_list",dest="nb_list", default=True,action="store_true", help="Use neighbor list ")

    
    (options, args) = parser.parse_args()
        
    return options, args
   
def main():
    #
    # Caclulate angles   k-i-j-l
    #
    
    import os, sys, numpy , math 
    import datetime
    import  gromacs, elements, xmol, prop, file_io, groups, top   #, vectors 
    
    debug = 1
    p_time = 1
    
    options, args = get_options()
    
    # Check options
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
    # Define rdf groups 
    #
    id_list_k  = options.id_k.split()
    id_list_i  = options.id_i.split()
    id_list_j  = options.id_j.split()
    id_list_l  = options.id_l.split()
    # 
    # Print initial properties 
    # 
    if( options.verbose ):
	print "  Group k types "	    
	for id_indx in range( len(id_list_k)):
	    print "      ",id_list_k[id_indx]
	print "  Group i types "	    
	for id_indx in range( len(id_list_i)):
	    print "      ",id_list_i[id_indx]
	print "  Group j types "	    
	for id_indx in range( len(id_list_j)):
	    print "      ",id_list_j[id_indx]  
	print "  Group l types "	    
	for id_indx in range( len(id_list_l)):
	    print "      ",id_list_l[id_indx]

    #
    # Create neighbor list
    #
    if( options.nb_list ):
	print "  Creating neighbor list"
	NBLIST, NBINDEX = top.build_covnablist(ELN,R)	    
    #
    # Find  atom groups k-i-j
    #
    angle_list = []
    #
    sum_angles = 0
    #
    # Find atom indices  of group i and j
    #
    for atom_i in range( len(ASYMB) ):
	add_i = 0 
	for id_indx in range( len(id_list_i)):
	    if( id_list_i[id_indx] == ATYPE[atom_i].strip() ): add_i = 1
	    if( id_list_i[id_indx] == GTYPE[atom_i].strip() ): add_i = 1
	if( add_i ):
	    
	    N_o = NBINDEX[atom_i]
	    N_f = NBINDEX[atom_i+1] - 1
	    
	    for indx_j in range( N_o,N_f+1):
		atom_j = NBLIST[indx_j]
		add_j = 0 
		for id_indx in range( len(id_list_j)):
		    if( id_list_j[id_indx] == ATYPE[atom_j].strip() ): add_j = 1
		    if( id_list_j[id_indx] == GTYPE[atom_j].strip() ): add_j = 1
		if( add_j ):
		

		    Nj_o = NBINDEX[atom_j]
		    Nj_f = NBINDEX[atom_j+1] - 1
		    
					
		    for indx_k in range( N_o,N_f+1):
			atom_k = NBLIST[indx_k]
			add_k = 0 
			for id_indx in range( len(id_list_k)):
			    if( id_list_k[id_indx] == ATYPE[atom_k].strip() ): add_k = 1
			    if( id_list_k[id_indx] == GTYPE[atom_k].strip() ): add_k = 1
			if( add_k ):
			    
			    
			    for indx_l in range( Nj_o,Nj_f+1):
				atom_l = NBLIST[indx_l]
				add_l = 0 
				for id_indx in range( len(id_list_l)):
				    if( id_list_l[id_indx] == ATYPE[atom_l].strip() ): add_l = 1
				    if( id_list_l[id_indx] == GTYPE[atom_l].strip() ): add_l = 1
				if( add_l ):
					    
				    if(debug):
				        angle_list.append( [atom_k,atom_i,atom_j,atom_l] )
					print " found angle ",atom_k,atom_i,atom_j,atom_l
	
	    
    if( p_time ): t_i = datetime.datetime.now()

    #
    # Intialize lsits and counts
    #
    angle_max = 180.0
    n_bins = int(angle_max/options.bin_size)
    hist_cnt = numpy.zeros(n_bins+1)    
    #
    # Loop over frames
    #
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
	    for indx_kij in angle_list:
		
		a_k = indx_kij[0]
		a_i = indx_kij[1]
		a_j = indx_kij[2]
		
		r_k = R_f[a_k]
		r_i = R_f[a_i]
		r_j = R_f[a_j]
		
		r_ik = r_k - r_i
		r_ij = r_j - r_i
		
		# dih 
		a_l = indx_kij[3]
		r_l = R_f[a_l]
		r_il = r_l - r_j
		
		angle_i = prop.getDihedral(r_k,r_i,r_j,r_l)
		#angle_kij =  prop.getAngle(r_ik,r_ij)
		
		
		bin_index = int( round( angle_i/options.bin_size) )
		hist_cnt[bin_index] += 1
		
		if( options.verbose ):
		    #print " # ",a_k,a_i,a_j,angle_i
		    print " # ",a_k,a_i,a_j,a_l,angle_i
		
	else:
	    print " Error frame ",frame_id," does not exist "
	
    if( p_time ): 
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
	print "    ti t_f ",t_i,t_f
	print "    time ",dt_min,dt_sec
	
			
    if( options.verbose ):
	print "      Finding averages "
    #
    # Find averages
    #
    total_cnts = numpy.sum( hist_cnt)
    box_vol_ave = numpy.average( volume_i )
    
    if( options.verbose ):
	print "   Frames ",frame_cnt
	print "   Total counts ",total_cnts
	print "   Average box volume ",box_vol_ave
	
    # Write output 
    #
    hist_file = options.hist_out
    F_out = open(hist_file,"w")
    F_out.write("# Frames %d %d " %  (options.frame_o,options.frame_f))
    F_out.write("\n#    Bin-size %f  " % (options.bin_size))
    F_out.write("\n#    Frames %d  " % (frame_cnt))
    F_out.write("\n#    Total_cnts %d  " % (total_cnts))
    F_out.write("\n#    Average Box Volume %f " % ( box_vol_ave) )
    F_out.write("\n#    ")
    F_out.write("\n# bin index ; cnt    ; cnt/frame  ")
    
    for bin_index in range( 0,n_bins):
	    
	hist_val = options.bin_size*float(bin_index)
	
	val_cnt = float( hist_cnt[bin_index] )
	
	cnt_fnorm =     val_cnt    /float(frame_cnt)
	
	F_out.write("\n  %d %f %f %f  " % (bin_index,hist_val,val_cnt,cnt_fnorm) )
	
    F_out.close()
	
	
if __name__=="__main__":
    main()
   
