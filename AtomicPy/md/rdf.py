#! /usr/bin/env python
"""

Radial distribution  code

 length - Angstroms
 mass   - AMU
 volume - Angstroms^3
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
    parser.add_option("-o","--output_id", dest="output_id", default="rdf",type="string",help=" prefix for output files  ")
    
    # Atom type options 
    parser.add_option("-i","--id_i", dest="id_i", type="string", default="", help="Atom types (FFtype,GROMACStype) of group i ")
    parser.add_option("-j","--id_j", dest="id_j", type="string", default="", help="Atom types (FFtype,GROMACStype) of group j ")


    parser.add_option("--filter_eln_i", dest="filter_eln_i", type="string", default="", help=" filter atoms by atomic number ")
    parser.add_option("--filter_fftype_i", dest="filter_fftype_i", type="string", default="", help=" filter atoms by force field type ")
    parser.add_option("--filter_residue_i", dest="filter_residue_i", type="string", default="", help=" filter atoms by residue name ")
    parser.add_option("--filter_unit_i", dest="filter_unit_i", type="string", default="", help=" filter atoms by unit name ")
    parser.add_option("--filter_cord_i", dest="filter_cord_i", type="string", default="", help=" filter atoms by cordination ")
    parser.add_option("--filter_mol_i", dest="filter_mol_i", type="string", default="", help=" filter atoms by molecule number  ")
    parser.add_option("--filter_lmptype_i", dest="filter_lmptype_i", type="string", default="", help=" filter atoms by lammps type ")

    parser.add_option("--filter_eln_j", dest="filter_eln_j", type="string", default="", help=" filter atoms by atomic number ")
    parser.add_option("--filter_fftype_j", dest="filter_fftype_j", type="string", default="", help=" filter atoms by force field type ")
    parser.add_option("--filter_residue_j", dest="filter_residue_j", type="string", default="", help=" filter atoms by residue name ")
    parser.add_option("--filter_unit_j", dest="filter_unit_j", type="string", default="", help=" filter atoms by unit name ")
    parser.add_option("--filter_cord_j", dest="filter_cord_j", type="string", default="", help=" filter atoms by cordination ")
    parser.add_option("--filter_mol_j", dest="filter_mol_j", type="string", default="", help=" filter atoms by molecule number  ")
    parser.add_option("--filter_lmptype_j", dest="filter_lmptype_j", type="string", default="", help=" filter atoms by lammps type ")

    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--in_itp", dest="in_itp",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--read_gros", dest="read_gros",   default=False,action="store_true" ,help="read serieries of gromacs .gro files ")
    
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
			sq_r_ij = prop.sq_drij_c(r_i,r_j,LV)
			if( sq_r_ij <= sq_r_cut ):
			    NNAB += 1
			    NBLIST[NNAB] =  atom_j
			    
    # Account for final atom, which has no connections 
    NBINDEX[atom_i+1] =  NNAB + 1
    NBINDEX[atom_i+2] =  NNAB + 1

    return NBLIST, NBINDEX

def get_lmp_frame(lammpsxyz_lines,NA_i,frame_i):

    R_i = []
    for atom_i in range( NA_i ):   
        line_cnt += 1
                    
        if( line_cnt > len(lammpsxyz_lines)-1):
            print " frame is missing some atoms ",atom_i," not found "
                        
            col =  lammpsxyz_lines[line_cnt].split()
            if( len(col) >= 4 ):
                type_i = int(col[0]) - 1
                r_x = float(col[1])
                r_y = float(col[2])
                r_z = float(col[3])
                R_i.append( numpy.array( [r_x,r_y,r_z] ) )
    return R_i

def main():
    """
    Calculate radial distribution based on specified atom types 
    """
    import os, sys, numpy , math 
    import datetime
    import gromacs, elements, xmol, prop, file_io, groups, lammps  #, vectors     
    import mpiNREL


    debug = 0
    
    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()
    
    options, args = get_options()
    
    #
    # Open output files 
    #
    # if( rank == 0 ):

    log_file = options.output_id + ".log"
    log_out = open(log_file,"w") 

    dat_file = options.output_id + ".dat"
    dat_out = open(dat_file,"w") 
    dat_out.write("#   Input ")

    p.barrier()
     
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
        if( rank == 0 ):
            if( options.verbose ): print  "     - Reading in ",options.in_data
            dat_out.write("\n#   Lammps data file %s  "% (options.in_data))
            
        ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,MOLNUMB,ATYPE_IND,CHARGES,R , ATYPE, BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV = lammps.read_data(options.in_data)
        
        AMASS = []
        for atom_i in range(len(ATYPE_IND)):
            type_ind  = ATYPE_IND[atom_i]
            AMASS.append( ATYPE_MASS[type_ind])
        ASYMB ,ELN = elements.mass_asymb(AMASS)
        
        #if(  options.ff_software == "gromacs"  ):
        GTYPE = []
        RESID = []
        RESN = []
        CHARN = []
        VEL = []
        CTYPE = []
        UNITNUMB = []
        UNITTYPE = []
        RING_NUMB  = []
        for i in range( len(ELN) ):
            GTYPE.append(ASYMB[i])
            RESID.append("MOL")
            CHARN.append(MOLNUMB[i])
            RESN.append(MOLNUMB[i])
            VEL.append( numpy.array( [0.0 ,0.0 ,0.0]) )
            CTYPE.append(MOLNUMB[i])
            UNITNUMB.append(MOLNUMB[i])
            UNITTYPE.append(MOLNUMB[i])
            RING_NUMB.append(MOLNUMB[i])
            #
        N_MOL,MOLPNT,MOLLIST = groups.molecule_list(MOLNUMB)
        NBLIST,NBINDEX = groups.build_nablist_bonds(ELN,BONDS)

        p.barrier()

    #
    # Get lammps xyz file 
    #
    if( len(options.in_lammpsxyz) ):

        if( rank == 0 ):
            dat_line = "\n#   Lammps xyz file %s  "% (options.in_lammpsxyz)
            if( options.verbose ):
                print  "     - Reading in ",options.in_lammpsxyz
            dat_out.write(dat_line)

        # Some issues with this read in
        #   will read in on processor 0 and broadcast
        n_frames =  options.frame_f
        p.barrier()
        if( rank == 0 ):
            lammpsxyz_F = open(options.in_lammpsxyz , 'r' )
            lammpsxyz_lines = lammpsxyz_F.readlines()
            lammpsxyz_F.close()
            n_frames = int( float( len(lammpsxyz_lines) )/ float( len(ASYMB) + 2) )
            lammpsxyz_line_cnt = -1

            if( options.verbose ):
                print  "       with ",n_frames," frames "

    # modify based on lammpsxyz read in 
    p.barrier()
    n_frames  = p.bcast( n_frames)
    p.barrier()

    if( options.frame_f == -1 ):
        options.frame_f  = n_frames - 1
        if( options.verbose and  rank == 0 ):
            print "  modify frame_f to ",options.frame_f," rank ",rank 

    # Calculate square of cut off
    #
    sq_r_cut = options.r_cut**2
    #
    # Calculate intial volume and total density 
    #
    vol_o = prop.volume( LV )
    NP = len( ELN )
    N_MOL = max(MOLNUMB) + 1 
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
    p.barrier()   
    
    #
    # Find atom indices  of group i and j
    #
    list_i, sum_i = groups.filter_atoms(options.filter_eln_i,options.filter_fftype_i,options.filter_residue_i,options.filter_unit_i,options.filter_cord_i,options.filter_mol_i,ASYMB,ELN,ATYPE,RESID,UNITTYPE,MOLNUMB,NBLIST,NBINDEX)

    
    list_j, sum_j = groups.filter_atoms(options.filter_eln_j,options.filter_fftype_j,options.filter_residue_j,options.filter_unit_j,options.filter_cord_j,options.filter_mol_j,ASYMB,ELN,ATYPE,RESID,UNITTYPE,MOLNUMB,NBLIST,NBINDEX)

    #
    # For lammps read in 
    #

    if(  len(options.filter_lmptype_i) ):

        #
        # Find atom indices  of group i and j
        #
        list_i = []
        sum_i = 0

        for atom_i in range(len(ASYMB) ):   
            add_atom = 0
            for f_id in options.filter_lmptype.split():
                lmp_t = int( f_id )
                if( ATYPE_IND[atom_i]+1 == lmp_t ):
                    add_atom = 1
            if( add_atom ):
                list_i.append( atom_i )
                sum_i += 1
    
    if(  len(options.filter_lmptype_j) ):

        #
        # Find atom indices  of group i and j
        #
        list_j = []
        sum_j = 0

        for atom_j in range(len(ASYMB) ):   
            add_atom = 0
            for f_id in options.filter_lmptype.split():
                lmp_t = int( f_id )
                if( ATYPE_IND[atom_j]+1 == lmp_t ):
                    add_atom = 1
            if( add_atom ):
                list_j.append( atom_i )
                sum_j += 1
    
    #
    # Print info 
    #
    if( rank == 0  ):

        t_i = datetime.datetime.now()
        log_line = "\n Start time " + str(t_i)
        log_out.write(log_line)

        dat_out.write("\n#   Date "+str(t_i))
        dat_out.write("\n#   Number of processors  "+str(size))
        dat_out.write("\n#   System ")
        dat_out.write("\n#     Molecules  %d "%(N_MOL))
        dat_out.write("\n#     Particles  %d "%(NP))
        dat_out.write("\n#   Options ")
        dat_out.write("\n#     bin size  %f "%(options.bin_size))
        dat_out.write("\n#     bins  %d "%(n_bins))
        dat_out.write("\n#     box length for max  %f "%(LV[0][0]))
        dat_out.write("\n#     Initial frame  %d "%(options.frame_o))
        dat_out.write("\n#     Step frame  %d "%(options.frame_step))
        dat_out.write("\n#     Final frame  %d "%(options.frame_f))
        dat_out.write("\n#   Output ")
        dat_out.write("\n#    Frame count; Frame number ; Average length (A); Standard deviation (A), box length (A)")

        print "   - Properties and options "
	print "     Cut off radius ",options.r_cut," angstroms"
	print "     Bin size ",options.bin_size," angstroms"
	print "     Number of bins ",n_bins
        print "     Number density for neighbor list ",NUMB_DENSITY," atoms/angstroms^3 "
	print "     N_i ",sum_i
	print "     N_j ",sum_j
	print "     Group i types "    
	for id_indx in range( len(id_list_i)):
	    print "         ",id_list_i[id_indx]
	print "     Group j types "
	for id_indx in range( len(id_list_j)):
	    print "         ",id_list_j[id_indx]
	       
    #
    # If multi-core split the number of atoms in the molecule onto each core
    #
    pointIndices = range( len(list_i)  )
    if( debug ): print rank, size," splitOnProcs "
    # Create a list of atomic indices for each processor 
    myChunk_i  = p.splitListOnProcs(pointIndices)
    debug = 0
    if(debug):                
        print " cpu ",rank ," has atoms ",myChunk_i[0]," - ",myChunk_i[len(myChunk_i)-1],"  \n"
	for atom_i in range(len(list_i) ):
	    print myChunk_i[atom_i],list_i[atom_i]
	sys.exit(" debug myChunk ")
	
    p.barrier()
     
    if( options.ptime ): t_i = datetime.datetime.now()
    #
    # Create neighbor list
    #
    if( options.nb_list ):
	print "  Creating neighbor list"
	NBLIST, NBINDEX = build_nablist(options,ELN,R,LV)

    #
    # Loop over frames
    #
    rdf_cnt_i = numpy.zeros(n_bins+1)    
    frame_cnt = 0
    volume_i = [] #numpy.array()
    NA_i = len(ELN)

    if( options.ptime ): 
	t_i = datetime.datetime.now()


    #
    for frame_i in range(options.frame_o,options.frame_f+1,options.frame_step):
        frame_read = False
        if( options.read_gros ):

            frame_id = "frames/n"+str(frame_i)+options.frame_sufx
            if( file_io.file_exists( frame_id ) ):
                if( options.verbose and rank == 0 ):
                    print "    - reading ",frame_id
                GTYPE, R_f, VEL, LV = gromacs.read_gro(options,frame_id)
                frame_read = True
                
        if( len(options.in_lammpsxyz) ):
            # R_f = get_lmp_frame(lammpsxyz_lines,len(ASYMB_sys),frame_i)

            R_f = []

            # Read in R_f from xyz file
            #   This is done on processor 0 to avoid broadcasting the very large lammpsxyz_lines around
            #   it would be preferable to read this in chunk by chunk 
            if( rank == 0 ):
                lammpsxyz_line_cnt = frame_i*(NA_i + 2 ) 
                lammpsxyz_line_cnt += 1
                if( options.verbose ):
                    log_line  = "\n    - reading frame %s  starting at line %d with comment %s " % (frame_i,lammpsxyz_line_cnt-1,lammpsxyz_lines[lammpsxyz_line_cnt] )
                    print log_line
                    log_out.write(log_line)

                r_max = -1000000.0 
                r_min = 1000000.0 
                for atom_i in range( NA_i ):   
                    lammpsxyz_line_cnt += 1

                    if( lammpsxyz_line_cnt > len(lammpsxyz_lines)-1):
                        print " frame is missing some atoms ",atom_i," not found "

                    col =  lammpsxyz_lines[lammpsxyz_line_cnt].split()
                    if( len(col) >= 4 ):
                        type_i = int(col[0]) - 1
                        
                        if( type_i != ( ATYPE_IND[atom_i] ) ):
                            print " Reference file not commpatable with ",options.in_lammpsxyz
                            print " atom ",atom_i+1," in reference type ",ATYPE_IND_sys[atom_i] +1 ," and ",type_i," in ",options.in_lammpsxyz
                            sys.exit(" Lammps atom types do not agree ")


                        r_x = float(col[1])
                        r_y = float(col[2])
                        r_z = float(col[3])
                        R_f.append( numpy.array( [r_x,r_y,r_z] ) )
                        if( r_x > r_max ): r_max = r_x
                        if( r_y > r_max ): r_max = r_y
                        if( r_z > r_max ): r_max = r_z
                        if( r_x < r_min ): r_min = r_x
                        if( r_y < r_min ): r_min = r_y
                        if( r_z < r_min ): r_min = r_z
                        
                # Estimate bax size based on max/min
                l_box = r_max -r_min
                LV[0][0] = l_box
                LV[1][1] = l_box
                LV[2][2] = l_box
               
                if( options.verbose  and rank == 0  ):
                    log_line  = "\n        with box length %f Angstroms estimated from max/min"%(l_box)
                    print log_line
                    log_out.write(log_line)

            # Broadcast R_f and LV  to all processors  
            p.barrier()
            R_f = p.bcast(R_f)
            LV = p.bcast(LV)
            p.barrier()

            frame_read = True
            
        if( frame_read ):

            volume_i.append(  prop.volume( LV ) )
            frame_cnt += 1 
            # xmol.print_xmol(ASYMB,R_i,file_xmol)
            #
            # Loop over lists
            #
            #
            for l_indx_i in myChunk_i:
                atom_i = list_i[l_indx_i]
                r_i = R_f[atom_i]
                #
                # using nieghbor list has no clear pereforance advantages if update every frame
                #  and not update every frame is beyond the intelegence of this program
                #
                #if( options.nb_list ):
                #    list_j = []
                #    N_o = NBINDEX[atom_i]
                #    N_f = NBINDEX[atom_i+1] - 1
                #    
                #    print atom_i
                #    print N_o,N_f
                #    
                #    for indx in range( N_o,N_f+1):
                #	atom_j = NBLIST[indx]
                #		    
                #	add_ij = 0
                #	for id_indx in range( len(id_j)):
                #	    if( id_j[id_indx] == ATYPE[atom_j].strip()  ): add_ij = 1
                #	    if( id_j[id_indx] == GTYPE[atom_j].strip()  ): add_ij = 1
                #	if( add_ij ):
                #	    list_j.append(atom_j)
                #
                for atom_j in list_j:
                    add_ij = 1
                    if( atom_i >= atom_j ):
                        add_ij = 0
                    #
                    # Check for intra vs. inter molecular
                    #
                    if( options.mol_inter ):
                        if( MOLNUMB[atom_i] == MOLNUMB[atom_j] ): add_ij = 0 
                    if( options.mol_intra ):
                        if( MOLNUMB[atom_i] != MOLNUMB[atom_j] ): add_ij = 0 
                    if( add_ij ):

                        r_j = R_f[atom_j]

                        if( options.cubic ):
                            sq_r_ij = prop.sq_drij_c(r_i,r_j,LV)
                        else:
                            sq_r_ij = prop.sq_drij(options,r_i,r_j,LV)

                        if( sq_r_ij <= sq_r_cut ):
                            m_ij = numpy.sqrt(sq_r_ij)
                            bin_index = int( round( m_ij/options.bin_size) )
                            rdf_cnt_i[bin_index] += 2
	
	else: 
	    if( rank == 0 ):
		print " Error frame ",frame_i," does not exist "
	    
	p.barrier() # Barrier for MPI_COMM_WORLD
	
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
    rdf_cnt = numpy.zeros(n_bins+1)   
    for bin_index in range( n_bins+1):
	# Sum rdf_cnt of each bin on each processor 
	rdf_cnt[bin_index] = p.allReduceSum(rdf_cnt_i[bin_index])
	
	p.barrier() # Barrier for MPI_COMM_WORLD
	
	if( debug ):
	    print "processor ",rank," bin index ",bin_index," has ",rdf_cnt_i[bin_index]
	
	
    
    total_cnts = numpy.sum( rdf_cnt)

    box_vol_ave = numpy.average( volume_i )
    vol_cut = 4.0*math.pi/3.0*options.r_cut**3
    n_shperes = float(sum_i)*float(frame_cnt)
    sphere_den_j = float(total_cnts)/vol_cut/n_shperes #/2.0  # N_B A^-3
    box_den_i = float(sum_i )/float(box_vol_ave)
    box_den_j = float(sum_j )/float(box_vol_ave)
    
    if( rank == 0 ):
	if( options.verbose ):
		
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
   
