#! /usr/bin/env python
"""

Polymer length  code

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
    parser.add_option("-o","--output_id", dest="output_id", default="poly_length",type="string",help=" prefix for output files  ")
    
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

    
    (options, args) = parser.parse_args()
        
    return options, args
   
   
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
    import gromacs, elements, xmol, prop, file_io, groups, lammps #, vectors     
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
   
    log_file = options.output_id + ".log"
    log_out = open(log_file,"w") 
     
    dat_file = options.output_id + ".dat"
    dat_out = open(dat_file,"w") 
    dat_out.write("#   Input ")
     
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
        dat_out.write("\n#   Gromacs top file %s  "% (options.in_top))
    #
    # Get coord
    #
    if( len(options.in_gro) ):
        if( options.verbose and rank == 0 ): print "     Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options,options.in_gro)
        dat_out.write("\n#   Gromacs gro file %s  "% (options.in_gro))
    #
    # Get lammps data file 
    #
    if( len(options.in_data) ):
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
        for i in range( len(ELN) ):
            GTYPE.append(ASYMB[i])
            RESID.append("MOL")
            CHARN.append(MOLNUMB[i])
            RESN.append(MOLNUMB[i])
            #
        N_MOL,MOLPNT,MOLLIST = groups.molecule_list(MOLNUMB)


    #
    # Get lammps xyz file 
    #
    if( len(options.in_lammpsxyz) ):

        if( options.verbose ):
            print  "     - Reading in ",options.in_lammpsxyz
        dat_out.write("\n#   Lammps xyz file %s  "% (options.in_lammpsxyz))

        lammpsxyz_F = open(options.in_lammpsxyz , 'r' )
        lammpsxyz_lines = lammpsxyz_F.readlines()
        lammpsxyz_F.close()
        n_frames = int( float( len(lammpsxyz_lines) )/ float( len(ASYMB) + 2) )
        lammpsxyz_line_cnt = -1

        if( options.verbose ):
            print  "       with ",n_frames," frames "

        
    #
    # Intialize lsits and counts
    #
    n_bins = int(LV[0][0]/options.bin_size)
    polyl_bin = numpy.zeros(n_bins+1)    
    N_MOL = max(MOLNUMB) + 1 
    NA_i = len(ELN)
    p.barrier()

    #
    # If multi-core split the molecules onto each core
    #
    pointIndices = range( N_MOL )
    if( debug ): print rank, size," splitOnProcs "
    # Create a list of atomic indices for each processor 
    myChunk_i  = p.splitListOnProcs(pointIndices)

    #
    # Write input information 
    #
    t_i = datetime.datetime.now()
    log_line = " Start time " + str(t_i)
    log_out.write(log_line)

    dat_out.write("\n#   Date "+str(t_i))
    dat_out.write("\n#   System ")
    dat_out.write("\n#     Molecules  %d "%(N_MOL))
    dat_out.write("\n#     Atoms  %d "%(NA_i))
    dat_out.write("\n#   Options ")
    dat_out.write("\n#     bin size  %f "%(options.bin_size))
    dat_out.write("\n#     bins  %d "%(n_bins))
    dat_out.write("\n#     box length for max  %f "%(LV[0][0]))
    dat_out.write("\n#     Initial frame  %d "%(options.frame_o))
    dat_out.write("\n#     Step frame  %d "%(options.frame_o))
    dat_out.write("\n#     Final frame  %d "%(options.frame_o))
    dat_out.write("\n#   Output ")
    dat_out.write("\n#    Frame count; Frame number ; Average length (A); Standard deviation (A)")
    
    #
    # Loop over frames 
    #
    frame_cnt = 0
    for frame_i in range(options.frame_o,options.frame_f+1,options.frame_step):

        if( len(options.in_lammpsxyz) ):
            # R_f = get_lmp_frame(lammpsxyz_lines,len(ASYMB_sys),frame_i)

            R_f = []
            lammpsxyz_line_cnt += 2
            
            if( options.verbose ):
                log_line  = "\n    - reading frame %s  starting at line %d with comment %s " % (frame_i,lammpsxyz_line_cnt-1,lammpsxyz_lines[lammpsxyz_line_cnt] )
                print log_line
                log_out.write(log_line)
                               
            for atom_i in range( NA_i ):   
                lammpsxyz_line_cnt += 1

                if( lammpsxyz_line_cnt > len(lammpsxyz_lines)-1):
                    print " frame is missing some atoms ",atom_i," not found "

                col =  lammpsxyz_lines[lammpsxyz_line_cnt].split()
                if( len(col) >= 4 ):
                    type_i = int(col[0]) - 1
                    r_x = float(col[1])
                    r_y = float(col[2])
                    r_z = float(col[3])
                    R_f.append( numpy.array( [r_x,r_y,r_z] ) )
                    
        frame_cnt += 1 
        #
        # Loop over all the molecule
        #
        debug = 0
        if(debug):
            print "  Looping over ",N_MOL," molecules"
            print "   length of MOLPNT ",len(MOLPNT)

        poly_len_i = numpy.zeros((N_MOL))
        poly_std_i = numpy.zeros((N_MOL))
        
        for mol_i in myChunk_i:
            M_o = MOLPNT[mol_i]
            M_f = MOLPNT[mol_i+1] - 1
            n_atoms_mol = M_f - M_o + 1
            if(debug):
                print " Mol ",mol_i," has  ",n_atoms_mol," atoms "

            #
            # Find the maximum separation between atoms of a give molecule 
            #
            r_sq_moli =  numpy.zeros((n_atoms_mol)) 
            a_cnt = -1
            for indx_i in range( M_o,M_f+1):
                atom_i = MOLLIST[indx_i]
                if( MOLNUMB[atom_i] == mol_i + 1  ):
                    log_line =  "\n \n Error atom  is in  not listed \n \n" % (atom_i,MOLNUMB[atom_i] ,mol_i )
                    print log_line
                    sys.exit(" Error in molecule list")
                    
                a_cnt += 1 
                r_i = R_f[atom_i]

                for indx_j in range( M_o,M_f+1):
                    atom_j = MOLLIST[indx_j]
                    if( MOLNUMB[atom_j] == mol_i + 1  ):
                        log_line =  "\n \n Error atom  is in  not listed \n \n" % (atom_j,MOLNUMB[atom_j] ,mol_i )
                        print log_line
                        sys.exit(" Error in molecule list")
                    r_j = R_f[atom_j]

                    # compute the square of the speration 
                    sq_r_ij = prop.sq_drij_c(r_i,r_j,LV)
                    # store in numpy array 
                    r_sq_moli[a_cnt] = sq_r_ij
                    

            # Find maximum speration and use that as the molecular length 
            r_sq_max = numpy.amax(r_sq_moli)
            r_max = numpy.sqrt(r_sq_max)
            # store in numpy array 
            poly_len_i[mol_i] = r_max
            # Compute the standard deviation of the intra-molecular separations
            #   Not need currently 
            #   poly_std_i[mol_i] =   numpy.sqrt( numpy.std(r_sq_moli) )

            if( options.verbose ):
                log_line="\n        - Molecule %s has a length of %f Angstroms " %( mol_i,r_max)
                print log_line
                log_out.write(log_line)
                
	p.barrier() # Barrier for MPI_COMM_WORLD
        # Gather list form each processor 
        poly_len = p.gatherList(poly_len_i)
        # Calculate average and std
        polyl_ave = numpy.mean(poly_len)
        polyl_std = numpy.std(poly_len)
        
        if( options.verbose ):
            log_line="\n      - Average length %f with standard diveation %f Angstroms "%(polyl_ave,polyl_std)
            print log_line
            log_out.write(log_line)
        dat_line="\n %d %d %f %f "%(frame_cnt,frame_i,polyl_ave,polyl_std)
        dat_out.write(dat_line)


    t_f = datetime.datetime.now()
    dt_sec  = t_f.second - t_i.second
    dt_min  = t_f.minute - t_i.minute
    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
    log_line="\n  Finished time  " + str(t_f)
    log_out.write(log_line)
    log_line="\n  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "
    log_out.write(log_line)

    
    dat_out.close()
    log_out.close()
    
if __name__=="__main__":
    main()
   
