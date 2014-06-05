#! /usr/bin/env python
"""

Radious of gyration  code

 length - Angstroms
 mass   - AMU
 volume - Angstroms^3
 
 R_g = ( sum_i ( || r_i -r_mol_cms ||^2 m_i  ) / sum_i ( m_i ) )^(1/2)
 
 r_i - is the position of atom i
 r_mol_cms - the center of mass of the atoms i in the considered molecule
 
"""

# Dr. Travis Kemper
# Initial Date April 2014
# travis.kemper@nrel.gov

def get_options():
    """
    Set options
    """
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-p","--ptime", dest="ptime", default=False,action="store_true", help="Print performance information  ")
    
    # Atom type options 
    parser.add_option("-i","--id_i", dest="id_i", type="string", default="", help="Atom types (FFtype,GROMACStype) of group i ")

    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--in_itp", dest="in_itp",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")

    parser.add_option("--rdf_out", dest="rdf_out", type="string", default="rdf.dat", help="Output rdf file ")
 
    parser.add_option("--r_cut", dest="r_cut", type=float, default=10.0, help=" Cut off radius in angstroms ")
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.10, help=" Bin size in angstroms")

    parser.add_option("--cubic",dest="cubic", default=False,action="store_true", help="Use cubic pbc's for speed up ")

    parser.add_option("--mol_inter",dest="mol_inter", default=False,action="store_true", help="Use only inter molecular rdf's ")
    parser.add_option("--mol_intra",dest="mol_intra", default=False,action="store_true", help="Use only intra molecular rdf's")
    
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=0, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    
    parser.add_option("--frame_sufx", dest="frame_sufx", type="string", default=".gro", help=" sufix of frame file ")

    parser.add_option("--n_den", dest="n_den", type=float, default=1.0, help=" Reference number density N/A^3")
    parser.add_option("--nb_list",dest="nb_list", default=False,action="store_true", help="Use neighbor list ")
    parser.add_option("--add_dr",dest="add_dr",type=float,  default=0 , help="Add length to the || r_i - r_mol_cms || value ")

    
    (options, args) = parser.parse_args()
        
    return options, args
   

def main():
    """
    Calculate radial distribution based on specified atom types 
    """
    import os, sys, numpy , math 
    import datetime
    import  gromacs, elements, xmol, prop, file_io, groups  #, vectors     
    import mpiNREL


    debug = 0
    
    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()
    
    options, args = get_options()
    
    # Check options
    if( options.mol_inter and options.mol_intra and rank == 0  ):
	print " Options --mol_inter and --mol_intra are mutually exclusive "
	sys.exit("Error is specified options ")
    
    if( options.verbose and rank == 0  ):
	print "   Reading in files to establish intial conditions and atom types "
    #
    # Read in top file
    #
    if( len(options.in_top) ):
        if( options.verbose and rank == 0 ): print "      Reading in ",options.in_top
        ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST = gromacs.read_top(options,options.in_top)
        ASYMB , ELN  = elements.mass_asymb(AMASS)
    #
    # Get coord
    #
    if( len(options.in_gro) ):
        if( options.verbose and rank == 0 ): print "     Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options,options.in_gro)
    #
    # Get lammps data file 
    #
    if( len(options.in_data) ):
        if( options.verbose ): print  "     - Reading in ",options.in_data
        ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,RESN,ATYPE_IND,CHARGES,R , ATYPE, BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV = lammps.read_data(options.in_data)
        
        AMASS = []
        for atom_i in range(len(ATYPE_IND)):
            type_ind  = ATYPE_IND[atom_i]
            AMASS.append( ATYPE_MASS[type_ind])
        ASYMB ,ELN = elements.mass_asymb(AMASS)
        
        #if(  options.ff_software == "gromacs"  ):
        GTYPE = []
        RESID = []
        CHARN = []
        for i in range( len(ELN) ):
            GTYPE.append(ASYMB[i])
            RESID.append("MOL")
            CHARN.append(RESN[i])
            
    #
    # Get lammps xyz file 
    #
    if( len(options.in_lammpsxyz) ):
        if( options.verbose ): print  "     - Reading in ",options.in_lammpsxyz

        lammpsxyz_F = open(options.in_lammpsxyz , 'r' )
        lammpsxyz_lines = lammpsxyz_F.readlines()
        lammpsxyz_F.close()        
        n_frames = int( float( len(lammpsxyz_lines) )/ float( len(ASYMB) + 2) )

	line_cnt = -1
	
	for frame_i in range(n_frames):
	    
	    line_cnt += 1
	    line_cnt += 1
	    
	    if( options.verbose ):
		print " reading frame ",frame_i," starting at line ",line_cnt-1," with comment ",lammpsxyz_lines[line_cnt] 
	    
	    for atom_i in range(len(ASYMB) ):   
		line_cnt += 1
		
		if( line_cnt > len(lammpsxyz_lines)-1):
		    print " frame is missing some atoms ",atom_i," not found "
		    # sys.exit("read in e)
		    
		col =  lammpsxyz_lines[line_cnt].split()
		if( len(col) >= 4 ):
		    type_i = int(col[0]) - 1
		    r_x = float(col[1])
		    r_y = float(col[2])
		    r_z = float(col[3])
		    
	sys.exit(" this function is still being implemented, sorry yo")
	
    #
    # Define rdf groups 
    #
    id_list_i  = options.id_i.split()
    p.barrier()       
    #
    # Find atom indices  of group i 
    #
    ind_i = []
    list_i = []
    cnt_i = -1
    molindx_list = []
    #
    sum_i = 0
    #
    # Loop over all the molecule find specified types each molecule
    #   place in neighbor list 
    #
    MOL_CNT = max(MOLNUMB)
    max_length = 0.0
    
    #
    # Square of length to add to r_dmass
    #
    sq_add_dr = options.add_dr*options.add_dr
    
    debug = 0
    if(debug):
	print "  Looping over ",MOL_CNT," molecules"
    
    for mol_i in range(MOL_CNT):
        M_o = MOLPNT[mol_i]
        M_f = MOLPNT[mol_i+1] - 1
	cnt_intial = cnt_i
	ind_i.append( cnt_i + 1 )
	if(debug):
	    print " Mol ",mol_i," has  ",M_f - M_o + 1
	    
	for indx in range( M_o,M_f+1):
	    atom_i = MOLLIST[indx]
	    
	    add_i = 0 
	    for id_indx in range( len(id_list_i)):
		if( id_list_i[id_indx] == ATYPE[atom_i].strip() ): add_i = 1
		if( id_list_i[id_indx] == GTYPE[atom_i].strip() ): add_i = 1
	    if( add_i ):
		cnt_i += 1 
		list_i.append( atom_i )
		sum_i += 1
		r_i = R[atom_i]
		#
		# loop over all the other atoms in the molecule
		#   to find maximum inter-molecular atomic speration 
		#
		#
		#consider_mol = cnt_i - cnt_intial
		#if( consider_mol > 1 ):
		#    N_o = ind_i[mol_i]
		#    N_f = ind_i[cnt_i] 
		#    for indi_j in range( N_o,N_f):
		#	atom_j = list_i[indi_j]
		#	r_j = R[atom_j]
		#	sq_dr_massi = prop.sq_drij_c(r_i,r_j,LV)
		#	if( sq_dr_massi > max_length ):
		#	    print " new max length between ",atom_i,atom_j,numpy.sqrt(sq_dr_massi)
		#	
		#	
		
	consider_mol = cnt_i - cnt_intial
	if( consider_mol > 0 ):
	    
	    molindx_list.append(mol_i) # for mpi splitting
	    
	    if( debug):
		print "      with  ",consider_mol," considered atoms "
		

    #
    # Intialize lsits and counts
    #
    # !!! Hack !!! 
    max_length = 100.0 # hard set for now
    
    n_bins = int(max_length/options.bin_size)		
    #
    # If multi-core split the number of molecules onto each core
    #
    #pointIndices = range( len(molindx_list)  )
    if( debug ): print rank, size," splitOnProcs "
    # Create a list of atomic indices for each processor 
    myChunk_i  = p.splitListOnProcs(molindx_list)
    debug = 0
    if(debug):                
        print " cpu ",rank ," has atoms ",myChunk_i[0]," - ",myChunk_i[len(myChunk_i)-1],"  \n"
	for atom_i in range(len(list_i) ):
	    print myChunk_i[atom_i],list_i[atom_i]
	sys.exit(" debug myChunk ")
	
    p.barrier()
    
    #
    # Print info 
    #
    if( rank == 0  ):
        print "   - Initial properties and options "
	
    if( options.ptime ): t_i = datetime.datetime.now()
    #
    # Loop over frames
    #
    r_gy_hist = numpy.zeros(n_bins+1)    
    frame_cnt = 0
    
    rg_t_file  = "rg.dat"
    tf = open(rg_t_file,"w")
    tf.write(" # frame , rg ")
    if( options.ptime ): 
	t_i = datetime.datetime.now()    
    #
    for frame_i in range(options.frame_o,options.frame_f+1,options.frame_step):
	frame_id = "frames/n"+str(frame_i)+options.frame_sufx
	if( file_io.file_exists( frame_id ) ):
	    GTYPE, R_f, VEL, LV = gromacs.read_gro(options,frame_id)
	    frame_cnt += 1

	    sum_dr_massi_mass = 0 
	    sum_mass = 0
		    
	    if( options.verbose and rank == 0 ):
		print "    - reading ",frame_id
    	    #
	    # Loop over lists
	    #
	    #
	    for mol_i in molindx_list:
		
		
		M_o = ind_i[mol_i]
		M_f = ind_i[mol_i+1] - 1
		
		print "    mol ",mol_i,M_o,M_f
		
		
		# r_mass = prop.list_cent_mass2(list_i,mol_i,R_f,AMASS)
	    
		prop_dim = len(R_f[0])
	    
		# Intialize center of mass
		total_mass = 0.0
		#r_mass.append(r_zero)
		r_mass = numpy.array( [0.0,0.0,0.0] )
			    
		for indx in range( M_o,M_f+1):
		    atom_i = list_i[indx]
		    a_mass = AMASS[atom_i]
		    r_i = R_f[atom_i]
				
		    # sum center of mass
		    total_mass = total_mass + a_mass
	    
		    for d in range(prop_dim):
			r_mass[d] += a_mass*R_f[atom_i][d]
			    
		# Normalize 
		for d in range(prop_dim):
		    r_mass[d] = r_mass[d]/total_mass
		    
		    
		#r_mass = prop.list_cent_mass2(MOLLIST,mol_i,R_f,AMASS)
		print " center of mass of mol ",mol_i,r_mass
		
		for indx in range( M_o,M_f+1):
		    atom_i = list_i[indx]
		    a_mass = AMASS[atom_i]
		    r_i = R_f[atom_i]
		    
		    sq_dr_massi = prop.sq_drij_c(r_i,r_mass,LV)
		    # 
		    
		    if( options.add_dr > 0 ):
			dr_massi = numpy.sqrt(sq_dr_massi)
		     	dr_massi += options.add_dr
			sq_dr_massi = dr_massi*dr_massi
			
		    sum_dr_massi_mass += sq_dr_massi*a_mass
		    sum_mass += a_mass
		    
	    r_gy =  numpy.sqrt( sum_dr_massi_mass / sum_mass )
	    
	    # print frame_i,r_gy
	    tf.write("\n %f %f " % (frame_i,r_gy ) )
	    
	    bin_index = int( round( r_gy/options.bin_size) )
	    r_gy_hist[bin_index] += 1
		    
	else: 
	    if( rank == 0 ):
		print " Error frame ",frame_id," does not exist "
	    
	p.barrier() # Barrier for MPI_COMM_WORLD
    
    
    tf.close()
   
    if( options.ptime and rank == 0 ): 
	t_f = datetime.datetime.now()
	dt_sec  = t_f.second - t_i.second
	dt_min  = t_f.minute - t_i.minute
	if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
	print "  rdf ti t_f ",t_i,t_f
	print "  rdf time ",dt_min,dt_sec
	
			
    if( options.verbose and rank == 0 ):
	print "      Finding averages "
    #
    # Find averages
    #
    debug = 0

    total_cnts = sum( r_gy_hist )

    if( rank == 0 ):
	if( options.verbose ):
		
	    print "   Frames ",frame_cnt
	    print "   Total counts ",total_cnts
		
	rg_h_file  = "rg.hist"
	hf = open(rg_h_file,"w")
	hf.write("# RDF frames %d %d " %  (options.frame_o,options.frame_f))
	hf.write("\n#    Bin-size %f  " % (options.bin_size))
	hf.write("\n#    Frames %d  " % (frame_cnt))
	hf.write("\n#    Total_cnts %d  " % (total_cnts))
	hf.write("\n#    N_i %d " % (sum_i ))
	hf.write("\n#    rg , prop ")
	
	for bin_index in range( 1,n_bins):
	    r_val = options.bin_size*float(bin_index)
	    r_gy_prob = float( r_gy_hist[bin_index] )/float(total_cnts)    
	    hf.write("\n %f %f " % (r_val,r_gy_prob ) )
	    
	hf.close()
	
if __name__=="__main__":
    main()
   
