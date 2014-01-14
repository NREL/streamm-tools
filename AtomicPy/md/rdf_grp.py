#! /usr/bin/env python
# Convert gromacs structure file into NWCHEM files for VAB calculation

# Dr. Travis Kemper
# NREL
# Initial Date 12/18/2013
# Email travis.kemper@nrel.gov
# Version 2.00 

def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] "
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False, help="Verbose output ")

    parser.add_option("--nproc", dest="nproc", type=int, default=24, help=" Number of processors for mpi run ")

    parser.add_option("--user", dest="user", type="string",  default="tkemper", help=" User name  ")

    parser.add_option("--r_cut", dest="r_cut", type=float, default=10.0, help=" Cut off radius in angstroms ")
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.10, help=" Bin size in angstroms")

    parser.add_option("--cubic",dest="cubic", default=False,action="store_true", help="Use cubic pbc's for speed up ")


    parser.add_option("--mol_inter",dest="mol_inter", default=False,action="store_true", help="Use only inter molecular rdf's ")
    parser.add_option("--mol_intra",dest="mol_intra", default=False,action="store_true", help="Use only intra molecular rdf's")
    
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=0, help=" Initial frame to read")
    
    parser.add_option("--frame_sufx", dest="frame_sufx", type="string", default=".gro", help=" sufix of frame file ")

    parser.add_option("--rdf_out", dest="rdf_out", type="string", default="rdf.dat", help="Output rdf file ")





    parser.add_option("--cluster_host", dest="cluster_host",type="string", default="peregrine",help=" name of cluster ")

    # System properties
    parser.add_option("--prop_dim", dest="prop_dim", type="int",default="3", help="Spacial dimension of the system")

    # Input output files     
    parser.add_option("-r","--ref_file", dest="ref_file", type="string", default="sim.ref", help=" Reference file used for supesquent reading of data ")
 
    # Gromacs
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file ")
 
    # FF settings
    parser.add_option("--ff_software", dest="ff_software", type="string", default="gromacs", help=" Software package to use for Force Field calculations  ")
    parser.add_option("--atom_types", dest="atom_types", type="string", default="", help="Read atom types that will replace default elements ")
    parser.add_option("--nlist_bonds", dest="nlist_bonds", default=True, help="Build neighbor list from bonds")
       
    # Gromacs
    parser.add_option("--gromacs_sufix", dest="gromacs_sufix", type="string", default="", help=" Sufix for gromacs calculations such as _mpi ")
    parser.add_option("--gromacs_dir", dest="gromacs_dir", type="string", default="", help=" gromacs dir ")
    parser.add_option("--g_center", dest="g_center",  default=False, help=" center minimization ")
    parser.add_option("--load_gromacs", dest="load_gromacs", type="string",  default="", help=" module comand to load gromacs ")

    # Groups
    parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
    parser.add_option("--rm_VS", dest="rm_VS", default=False, help=" remove virtual sites ")
    parser.add_option("--groups_com_cut", dest="groups_com_cut", type="float", default="10.0", help="Cut off for neighbors in QM output ")
    #parser.add_option("--h_term", dest="h_term", type="float", default="10.0", help="Cut off for neighbors in QM output ")

    parser.add_option("--debug", dest="debug", default=False, help=" Debug only produce sim.ref ")
        
    (options, args) = parser.parse_args()

    return options, args

def main():
    import os, sys, numpy , math 
    import gromacs,elements , top    , groups , prop, xmol , file_io 

    options, args = get_options()
    
    if( len(options.in_gro) ):
        if( options.verbose ): print "  Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options,options.in_gro)
        
    # Read in gro file
    if( len(options.in_top) ):
        if( options.verbose ): print "  Reading in ",options.in_top
        ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST  = gromacs.read_top(options,options.in_top)
        ASYMB , ELN  = elements.mass_asymb(AMASS)
    
    # Retype special atom types and replace element
    if( len( options.atom_types)): 
        if( options.verbose ): print "  Reading in ",options.atom_types
        ASYMB , ELN  = top.special_types(ATYPE,ASYMB , ELN , options.atom_types)
        
    # Print system information
    if( options.verbose ):
	print " prop "
        #prop.print_prop( AMASS,ELN,LV,CHARGES )
        ##top.print_prop( BONDS,ANGLES,DIH )
        
    # Create neighbor list
    if( options.nlist_bonds ):
        if( options.verbose ): print "  Creating neighbor list with bonds"
        NBLIST,NBINDEX  = groups.build_nablist_bonds(ELN,BONDS)

    #Find rings
    #if( options.find_rings ):
    #   RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)

    # Find groups
    if( options.group_ptma ):
        if( options.verbose ): print "  Find groups of TEMPO "
        group_index_i,group_list_i,group_numb,group_cnt = groups.tempo(  ATYPE,ELN,NBLIST,NBINDEX, options )

    #group_cnt = len(group_index_i) - 1
    # Find center of mass of groups
    if( len(group_index_i) > 0 ):
        group_cent = groups.cent_mass( group_index_i,group_list_i, R, AMASS,options )
        
    # Create neighbor list 
    g_nb_ind,g_nb_list,g_r_list = groups.nablist(group_cnt,group_cent,LV,options)

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
    
    sum_i = group_cnt
    sum_j = group_cnt
    
    
    debug = 1
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
	    group_cent = groups.cent_mass( group_index_i,group_list_i, R_f, AMASS,options )

	    volume_i.append(  prop.volume( LV ) )
	    frame_cnt += 1 
	    # xmol.print_xmol(ASYMB,R_i,file_xmol)
	    if( options.verbose ):
		print "reading ",frame_id
	    #R_frames.append( R_f )
		    
	    for g_i in range(group_cnt ):
		
		N_o = group_index_i[g_i]
		N_f = group_index_i[g_i + 1 ] - 1
		NA_g = N_f - N_o + 1
		if( options.verbose ):
		    print "        Printing input for group_i ",g_i," with atoms ",NA_g
		    
		g_o = g_nb_ind[g_i]
		g_f = g_nb_ind[g_i+1 ] 
		
		#for a_indx in range( N_o,N_f):
		#    i = group_list_i[a_indx]
		#    #if( int(ELN[i]) > 0 ):
	
		#print " group i has nb grps ",g_o,g_f
		
		for j_indx in range(g_o,g_f):
		    g_j = g_nb_list[ j_indx]
		    mag_dr = g_r_list[ j_indx]
		    
		    Nj_o = group_index_i[g_j]
		    Nj_f = group_index_i[g_j + 1 ] - 1
		    NAj_g = Nj_f - Nj_o + 1
            	    
		    if( options.verbose ): 
			print "          group_j ",g_j,mag_dr,NAj_g
			
		    # to prevent double counting with symmetric neutral geometries
		    #if( g_j > g_i ):
		    #n_id_j = n_calcid[g_j]
		    #c_id_j = c_calcid[g_j]
	
		    calc_id_j = "g_" + str(g_j)
		
		    bin_index = int( round( mag_dr/options.bin_size) )
		    rdf_cnt[bin_index] += 2
		    
		    if( debug):
			print "            group ",g_i,g_j ,mag_dr
    
	

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
	
	dr_rho = 2.0*dr_cnt_norm/dr_vol
	sphere_g = dr_rho/sphere_den_j/float( sum_i )
	box_g = dr_rho/box_den_j/float( sum_i )
	
	F_out.write("\n  %d %f %f %f %f %f %f " % (bin_index,r_val,dr_cnt_norm,dr_vol,dr_vol_apx,sphere_g,box_g) )
	
    F_out.close()
	
	    
	    
	    
if __name__=="__main__":
    main()
