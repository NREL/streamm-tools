#! /usr/bin/env python
"""
Create a supercell of randomly placed and rotated molecules
 length - Angstroms
 mass   - AMU
 volume - Angstroms^3
"""

# Dr. Travis Kemper
# Initial Date April 2014
# travis.kemper@nrel.gov

const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

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
    parser.add_option("-j","--in_json", dest="in_json", type="string", default="", help="Input json file")

    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--mult_i", dest="mult_i", type=int, default=1, help=" Number of molecules  ")

    parser.add_option("--sol_top", dest="sol_top", type="string", default="", help="Input gromacs topology file for solvent (.top) ")
    parser.add_option("--sol_gro", dest="sol_gro", type="string", default="", help="Input gromacs structure file for solvent (.gro) ")
    parser.add_option("--mult_s", dest="mult_s",  type=int, default=1, help=" Number of solvent molecules  ")
    parser.add_option("--lat_s", dest="lat_s",  default=False,action="store_true",  help=" Add solvent onto lattice to prevent the intra solvent overlap calculation  ")
    parser.add_option("--buf_s", dest="buf_s",  type=float, default=10.0, help=" intra solvent buffer " )

    parser.add_option("--box_l", dest="box_l", type=float, default=200.0, help=" Box size")
    
    parser.add_option("--in_itp", dest="in_itp",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")

    parser.add_option("--atomic_cut", dest="atomic_cut", type=float, default=2.5, help="Minimum distance between atoms of molecules ")

    parser.add_option("--den_target", dest="den_target", type=float, default=0.01, help="Target density g/cm^3 ")
    parser.add_option("--atoms_target", dest="atoms_target", type=int, default=100000, help="Target number of atoms ")
    parser.add_option("--max_mol_place", dest="max_mol_place", type=float, default=50, help="Maximum attempts to place a molecule  ")
    parser.add_option("--max_sys", dest="max_sys", type=float, default=3, help="Maximum system recreations at a certain lattice constant ")
    parser.add_option("--lc_expand", dest="lc_expand", type=float, default=0.100, help="Fraction of the box size to increase system size after max_sys is excieded ")

    parser.add_option("--out_gro", dest="out_gro", type="string", default="mol_system.gro", help="gromacs output file ")

    # Force field generation options     
    parser.add_option("--ff_software", dest="ff_software",type="string",default="lammps",help=" what software to use for the ff calculations   ")
    parser.add_option("--norm_dihparam", dest="norm_dihparam",default=0, help="Normalize dihedral potential terms if single dihedral is specified in itp file  ")
    
    (options, args) = parser.parse_args()
        
    return options, args
   
def ran_rot_shift(R_moli_c,lv_i ,ang_acc):
    import numpy , math , random
    
    prop_dim =3
    
    R_rs = []  # List of numpy arrays of rotated and shifted coordinates 
    #
    # Rotate molecule randomly 
    #
    rot_angle_i = float(random.randrange(0,ang_acc))*numpy.pi/float(ang_acc) # Random angle 1
    rot_angle_j = float(random.randrange(0,ang_acc))*numpy.pi/float(ang_acc) # Random angle 2
    
    cy = math.cos(rot_angle_i)
    sy = math.sin(rot_angle_i)
    cz = math.cos(rot_angle_j)
    sz = math.sin(rot_angle_j)
    
    for atom_i in range( len(R_moli_c) ):
        xd = R_moli_c[atom_i][0]
        yd = R_moli_c[atom_i][1]
        zd = R_moli_c[atom_i][2]
        
        
        r_x =  cy*cz*xd - sz*cy*yd + sy*zd 
        r_y =  sz*xd    + cz*yd            
        r_z = -sy*cz*xd + sy*sz*yd + cy*zd
        
        r_i =  numpy.array( [r_x,r_y,r_z] )
        R_rs.append(  r_i )
    
    
    #
    # Shift molecule to random point 
    #
    mol_origin = numpy.zeros(prop_dim)
    for d in range(prop_dim):
        mol_origin[d] = random.randrange(0,int(lv_i))

    for atom_i in range( len(R_moli_c) ):
        R_rs[atom_i] =  R_rs[atom_i] +  mol_origin 
        
    return R_rs

def main():
    """
    Read in structure information and determine number of molecules and volume to achieve desired density with specified
    number of atoms

    Returns: None
    """
    
    import os, sys, numpy , math , random
    import datetime
    import time    
    import gromacs, elements, xmol, prop, file_io, groups,lammps , top #, vectors
    import mpiNREL

    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    debug = 1
    
    #
    # Load information onto all processors 
    #
    options, args = get_options()
    prop_dim = 3
    ang_acc = 1000  # number of digets in random angle 
    
    #
    # Set debug options 
    #
    debug = 0       # Print debug statements 

 
    #
    # Read in json file
    #
    if( len(options.in_json) ):
        if( options.verbose ):
            print  "     - Reading in ",options.in_json
	json_data,json_success = jsonapy.read_jsondata(options.in_json)
	if(  json_success ):
	    #
            mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
            #
            # Need meta data to proceed 
            #      		    
            if( metadata_found ):
                ELN_i,ASYMB_i,CTYPE_i,CHARGES_i,UNITNUMB_i,UNITTYPE_i,R_i,VEL_i,ATYPE_i,AMASS_i,MOLNUMB_i,RING_NUMB_i,RESID_i,RESN_i,CHARN_i,json_atomicdata  = jsonapy.read_atomic(json_data)
                #
                # 
                #
                GTYPE_i = []
                for i in range( len(ELN_i) ):
                    GTYPE_i.append(ASYMB_i[i])
                
    #
    # Read in top file
    #
    if( len(options.in_top) ):
        if( options.verbose ): print  "     - Reading in ",options.in_top
        ATYPE_i,RESN_i,RESID_i,GTYPE_i,CHARN_i,CHARGES_i,AMASS_i,BONDS_i,ANGLES_i,DIH_i,MOLNUMB_i,MOLPNT_i,MOLLIST_i = gromacs.read_top(options,options.in_top)
        ASYMB_i,ELN_i  = elements.mass_asymb(AMASS_i)
    #
    # Get gro file 
    #
    if( len(options.in_gro) ):
        if( options.verbose ): print  "     - Reading in ",options.in_gro
        GTYPE_i,R_i,VEL_i,LV_i = gromacs.read_gro(options,options.in_gro)

    if( debug):
        print " Molecule "
        print "    ATYPE_i ",len(ATYPE_i)
        print "    ELN_i ",len(ELN_i)
        print "    R_i ",len(R_i)
         
    #
    # Read in top file
    #
    if( len(options.sol_top) ):
        if( options.verbose ): print  "     - Reading in ",options.sol_top
        ATYPE_s,RESN_s,RESID_s,GTYPE_s,CHARN_s,CHARGES_s,AMASS_s,BONDS_s,ANGLES_s,DIH_s,MOLNUMB_s,MOLPNT_s,MOLLIST_s = gromacs.read_top(options,options.sol_top)
        ASYMB_s,ELN_s  = elements.mass_asymb(AMASS_s)
    #
    # Get gro file 
    #
    if( len(options.sol_gro) ):
        if( options.verbose ): print  "     - Reading in ",options.sol_gro
        GTYPE_s,R_s,VEL_s,LV_s = gromacs.read_gro(options,options.sol_gro)            

    if( debug):
        print " solvent  "
        print "    ATYPE_s ",len(ATYPE_s)
        print "    ELN_s ",len(ELN_s)
        print "    R_s ",len(R_s)
         
    #
    # Read in ff file
    #
    if( len(options.in_itp) ):
        if( options.verbose ): print  "     - Reading in ",options.in_itp    
        FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(options.in_itp)
        if(  options.ff_software == "lammps"  ):
            # Identify total number of atom types for lammps output 
            ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND_i , BTYPE_REF, ANGTYPE_IND_i , ANGTYPE_REF, DTYPE_IND_i , DTYPE_REF = lammps.lmp_types(ELN_i,ATYPE_i,AMASS_i,BONDS_i,ANGLES_i,DIH_i)

            #   Build covalent nieghbor list for bonded information 
            NBLIST, NBINDEX = top.build_covnablist(ELN_i,R_i)
            
            # Check atom types to be sure each atom of the same type has the same number of neighbors 
            ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE_i,NBLIST,NBINDEX)
        
            ATYPE_EP, ATYPE_SIG = top.atom_parameters(options.in_itp,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
            BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(options.in_itp,BTYPE_IND_i , BTYPE_REF,FF_BONDTYPES)
            ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(options.in_itp,ANGTYPE_IND_i , ANGTYPE_REF,FF_ANGLETYPES)
            DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(options.in_itp, options.norm_dihparam, DTYPE_IND_i , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )
            
            IMPTYPE_F  = top.imp_parameters(options.in_itp)
    elif(len(options.in_data) == 0 ):
        if(  options.ff_software == "lammps"  ):
            print " An itp file specified with the --in_itp option is needed to create a lammps input file "
            sys.exit(" Read in error")
    #
    # Get lammps data file 
    #
    if( len(options.in_data) ):
        if( options.verbose ): print  "     - Reading in ",options.in_data
        ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,DIH_i,DTYPE_IND_i,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,RESN_i,ATYPE_IND_i,CHARGES_i,R_i , ATYPE_i, BONDS_i ,BTYPE_IND_i, ANGLES_i ,ANGTYPE_IND_i, LV_i = lammps.read_data(options.in_data)
        
        AMASS_i = []
        for atom_i in range(len(ATYPE_IND_i)):
            type_ind  = ATYPE_IND_i[atom_i]
            AMASS_i.append( ATYPE_MASS[type_ind])
        ASYMB_i ,ELN_i = elements.mass_asymb(AMASS_i)
        
        #if(  options.ff_software == "gromacs"  ):
        GTYPE_i = []
        RESID_i = []
        CHARN_i = []
        for i in range( len(ELN_i) ):
            GTYPE_i.append(ASYMB_i[i])
            RESID_i.append("MOL")
            CHARN_i.append(RESN_i[i])
            
        #CHARN_i = top.set_chargegroups(options,verbose,CG_SET,CHARN,ATYPE_i,ASYMB_i,ELN,R,NBLIST,NBINDEX, RING_NUMB,LAT_CONST)


    
    #
    # Shift molecule to have center of mass at origin 
    # 
    r_shift = numpy.array( [0.0,0.0,0.0] )
    R_moli_c = prop.shift_cent_mass(AMASS_i,R_i,r_shift)
    
    if( debug ):
        print  " molecule coordinates ", len(ELN_i) 
        for atom_i in range( len(ELN_i) ):
            print atom_i,R_i[atom_i], " -> ",R_moli_c[atom_i] 
    #
    # Calculate need molecules to achieve specified density and total number atoms 
    #
            
    target_density_amuang = options.den_target*const_avo/10.0
    mol_mass_amu = prop.total_mass( AMASS_i )
    mass_amu = mol_mass_amu*options.mult_i
    volume_target_ang = mass_amu/target_density_amuang    
    len_target_ang = volume_target_ang**(1.0/3.0)

            
    #
    # Set lattice vectors for pbc 
    #
    LV = numpy.zeros([3,3])
    
    LV[0,0] = options.box_l #len_target_ang
    LV[1,1] = options.box_l #len_target_ang
    LV[2,2] = options.box_l #len_target_ang

    #
    #  Find size of box for lattice method of solvent addtion 
    #
    if(options.lat_s ):

        #
        # Shift molecule to have center of mass at origin 
        # 
        r_shift = numpy.array( [0.0,0.0,0.0] )
        R_solvmol_c = prop.shift_cent_mass(AMASS_s,R_s,r_shift)
            
        # Find length of solvent molecule
        sq_mol_l = -1000000.0 
        for atom_i in range( len(ELN_s) ):
            r_i = R_solvmol_c[atom_i]
            
            for atom_j in range( len(ELN_s) ):
                r_j = R_solvmol_c[atom_j]
                sq_r_ij = prop.sq_drij_c(r_i,r_j,LV_s)
                if( sq_r_ij > sq_mol_l):
                    sq_mol_l = sq_r_ij
                    max_dr_i = atom_i
                    max_dr_j = atom_j
                    max_r_i = r_i 
                    max_r_j = r_j

        mol_l = numpy.sqrt(sq_mol_l)
        # Add a seperation to each molecule
        mol_l_lj = mol_l + options.buf_s   
        sq_mol_l_lj = mol_l_lj*mol_l_lj
        # Number molecules per side of box cube
        #  Leave extra space for non solvent molecules 
        mol_box_side = math.ceil( options.mult_s**(1.0/3.0) ) 
        sq_mol_box_side = mol_box_side*mol_box_side
        solv_box_l = mol_box_side*mol_l_lj #+ mol_l_lj

        log_line = " Solvent molecules  %d "%(options.mult_s)
        print log_line
	
        log_line = " Solvent molecules per cubic box side  %d "%(mol_box_side)
        print log_line
	
        log_line = " Reseting box size %f to %f for solvent lattice method "%(options.box_l,solv_box_l)
        print log_line
	

        LV[0,0] = solv_box_l
        LV[1,1] = solv_box_l
        LV[2,2] = solv_box_l
	
	#sys.exit(" box size debug  ")
	
    #
    # Initialize system
    #
    ASYMB_sys = []
    ELN_sys  = []
    ATYPE_sys = []
    RESN_sys = []
    RESID_sys = []
    GTYPE_sys = []
    CHARN_sys = []
    CHARGES_sys = []
    AMASS_sys = []
    GTYPE_sys = []
    R_sys = []
    VEL_sys = []
    BONDS_sys = []
    ANGLES_sys = []
    DIH_sys = []
    ATYPE_IND_sys = []
    BTYPE_IND_sys = []
    ANGTYPE_IND_sys = []
    DTYPE_IND_sys = []
        
    #
    # Set cut-off squared 
    #

    cut_ij_sq = options.atomic_cut* options.atomic_cut

    
    if( rank == 0  ):
        print "   - Tragets "
        print "     Input molecule has ",len(ELN_i)," atoms and mass of ",mol_mass_amu," AMU "
        print "     System target is ",options.atoms_target," atoms "
        print "     Input molecule will be multiplied ",options.mult_i," times "
        print "     The target density of ",options.den_target," g/cnm^3 ",target_density_amuang," AMU Angstroms^-3"
        print "     Cubic unit cell of ",len_target_ang," Angstroms will be used "
        # print "     Giving a molecular volume of ",mol_vol," mol/Angstrom^3 and a molecular unit length of ",mol_unit_l," Angstorms "
        print "   - Options "
        print "     ff_software ",options.ff_software
        print "     in_top ",options.in_top
        print "     in_gro ",options.in_gro
        print "     atomic_cut ",options.atomic_cut
        print "     den_target ",options.den_target
        print "     atoms_target ",options.atoms_target
        print "     max_mol_place ",options.max_mol_place
        print "     max_sys ",options.max_sys
        print "     lc_expand ",options.lc_expand
        print "     out_gro ",options.out_gro
        print "     ptime ",options.ptime
    #
    # If multi-core split the number of atoms in the molecule onto each core
    #
    p.barrier()
    #
    # Place the atomic indices into list 
    # 
    pointIndices = range( len(ELN_i)  )
    if( debug ): print rank, size," splitOnProcs "
    # Create a list of atomic indices for each processor 
    myChunk  = p.splitListOnProcs(pointIndices)
    if(debug):                
        print " cpu ",rank ," has atoms ",myChunk[0]," - ",myChunk[len(myChunk)-1],"  \n"
                    
    p.barrier()

    sys_mol_n = 0    # Number of molecules add to the system
    sys_attempts = 0  # Number of times the system has been reset
    #
    # Start adding molecules to the system
    #
    add_mol = 1
    while ( add_mol and options.mult_i > 0 ):
        #
        # Initialize 
        #
        add_mol = 1
        overlap_sum = 1
        moladd_atempts = 0
        
        # Record intial time for pereformance testing 
        if( options.ptime ):t_i = datetime.datetime.now()
        
        while ( overlap_sum ):
            moladd_atempts += 1
            
            # Declare on all processors
            R_shift_o = []
                
            if ( rank == 0 ):
                #
                # Get coordinates of randomly rotated and shifted molecule
                #   on processor 0
                #
                R_shift_o = ran_rot_shift(R_moli_c,len_target_ang,ang_acc )
            
            #
            # Broadcast molecular position from processor 0 to all other processors 
            #
            R_shift = p.bcast(R_shift_o)
            #
            # Check molecules do not overlap
            #
            overlap = 0
            if( len(ELN_sys) > 0 ):
                for atom_i in myChunk:
                    # Loop over all the atoms of the molecule being placed
                    r_i = R_shift[atom_i]
                    for sys_atom in range(len(ELN_sys)):
                        # Loop over all the atoms of the system 
                        r_j = R_sys[sys_atom]
                        r_ij_sq = prop.sq_drij_c(r_i,r_j,LV)
                        if( r_ij_sq < cut_ij_sq ):
                            overlap = 1
            
                p.barrier() # Barrier for MPI_COMM_WORLD
                
            #
            # Reduce sum the overlap variable from all the processors
            #   if it is zero everywhere there was no overlap detected 
            #
            # overlap_sum = mpi.all_reduce(world,overlap, lambda x,y: x + y)
            overlap_sum = p.allReduceSum(overlap)
            
            if( overlap_sum ==  0 ):
                # If no overlap detected add molecule to the system 
                sys_mol_n += 1
                
                for atom_i in range( len(ELN_i) ):
                    ASYMB_sys.append( ASYMB_i[atom_i])
                    ELN_sys .append( ELN_i[atom_i])
                    ATYPE_sys.append( ATYPE_i[atom_i])
                    RESN_sys.append( RESN_i[atom_i])
                    RESID_sys.append( RESID_i[atom_i])
                    GTYPE_sys.append( GTYPE_i[atom_i])
                    CHARN_sys.append( CHARN_i[atom_i])
                    CHARGES_sys.append( CHARGES_i[atom_i])
                    AMASS_sys.append( AMASS_i[atom_i])
                    R_sys.append( R_shift[atom_i])
                    if(  options.ff_software == "lammps"  ):
                        ATYPE_IND_sys.append( ATYPE_IND_i[atom_i])
                    
                
                if( options.verbose ):
                    if( rank == 0  ):
                        print "      -  Molecule ",sys_mol_n," has been added to the system after ",moladd_atempts," placment attempts "

                        if( options.ptime ):
                            t_f = datetime.datetime.now()
	                    dt_sec  = t_f.second - t_i.second
	                    dt_min  = t_f.minute - t_i.minute
	                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                            print "        - with placement time ",dt_min," min ",dt_sec," seconds "
                            

            if( moladd_atempts >= options.max_mol_place ):
                # If attempts to place molecule into the system exceed max set by max_mol_place
                #   reset system and star over 
                if( options.verbose ):
                    if( rank == 0  ):
                        print "        -  Attempts to add molecule ",sys_mol_n," has exceeded max attempts ",options.max_mol_place," system will be reset for the ",sys_attempts," time "
                    
                sys_mol_n = 0                
                sys_attempts += 1 
                ASYMB_sys = []
                ATYPE_IND_sys = []
                ELN_sys  = []
                ATYPE_sys = []
                RESN_sys = []
                RESID_sys = []
                GTYPE_sys = []
                CHARN_sys = []
                CHARGES_sys = []
                AMASS_sys = []
                R_sys = []
                VEL_sys = []

                if( sys_attempts >= options.max_sys  ):
                    # If the system has been reset over max_sys times expand the box size by lc_expand
                    da = LV[0,0]*options.lc_expand
                    db = LV[1,1]*options.lc_expand
                    dc = LV[2,2]*options.lc_expand
                    LV[1,1] = LV[1,1] + da
                    LV[1,1] = LV[1,1] + db
                    LV[2,2] = LV[2,2] + dc
                    if( options.verbose ):

                                            
                        if( rank == 0  ):
                            print '          - Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                            print '          - Lattice vectors will be expanded by (option lc_expand)'
                            print '              -  length added ',da,db,dc
                            
            p.barrier() # Barrier for MPI_COMM_WORLD
                    
                                    
        if( sys_mol_n ==  options.mult_i ):
            # If all the molecule have been added exit while loop and print system 
            add_mol = 0
            p.barrier() # Barrier for MPI_COMM_WORLD
            if( options.verbose and rank == 0  ):
                print " All molecules have been added "

    p.barrier()

    # Save initial configuration for restart 

    ASYMB_p = []
    ATYPE_IND_p = []
    ELN_p  = []
    ATYPE_p = []
    RESN_p = []
    RESID_p = []
    GTYPE_p = []
    CHARN_p = []
    CHARGES_p = []
    AMASS_p = []
    R_p = []
    VEL_p = []

    for atom_sys in range( len(ELN_sys) ):
        ASYMB_p.append( ASYMB_sys[atom_sys])
        ELN_p.append( ELN_sys[atom_sys])
        ATYPE_p.append( ATYPE_sys[atom_sys])
        RESN_p.append( RESN_sys[atom_sys])
        RESID_p.append( RESID_sys[atom_sys])
        GTYPE_p.append( GTYPE_sys[atom_sys])
        CHARN_p.append( CHARN_sys[atom_sys])
        CHARGES_p.append( CHARGES_sys[atom_sys])
        AMASS_p.append( AMASS_sys[atom_sys])
        R_p.append( R_sys[atom_sys])
        if(  options.ff_software == "lammps"  ):
            ATYPE_IND_p.append( ATYPE_IND_sys[atom_sys])
                    

    #
    # Shift solvent molecule to have center of mass at origin 
    # 
    r_shift = numpy.array( [0.0,0.0,0.0] )
    R_moli_c = prop.shift_cent_mass(AMASS_s,R_s,r_shift)
    
    if( debug ):
        print  " molecule coordinates ", len(ELN_i) 
        for atom_i in range( len(ELN_i) ):
            print atom_i,R_i[atom_i], " -> ",R_moli_c[atom_i] 
    #
    # Calculate need molecules to achieve specified density and total number atoms 
    #
            
    #
    # Place the atomic indices into list 
    # 
    pointIndices = range( len(ELN_s)  )
    if( debug ): print rank, size," splitOnProcs "
    # Create a list of atomic indices for each processor 
    myChunk  = p.splitListOnProcs(pointIndices)
    if(debug):                
        print " cpu ",rank ," has atoms ",myChunk[0]," - ",myChunk[len(myChunk)-1],"  \n"
                    
    p.barrier()

    sys_mol_n = 0    # Number of molecules added to the system
    sys_attempts = 0  # Number of times the system has been reset
    
    #
    # Start adding solvent to the system
    #
    add_mol = 1
    randomsolv = 0
    while ( add_mol and randomsolv and options.mult_s > 0 ):
        #
        # Initialize 
        #
        add_mol = 1
        overlap_sum = 1
        moladd_atempts = 0
        
        # Record intial time for pereformance testing 
        if( options.ptime ):t_i = datetime.datetime.now()
        
        while ( overlap_sum ):
            moladd_atempts += 1
            
            # Declare on all processors
            R_shift_o = []
                
            if ( rank == 0 ):
                #
                # Get coordinates of randomly rotated and shifted molecule
                #   on processor 0
                #
                R_shift_o = ran_rot_shift(R_moli_c,len_target_ang,ang_acc )
            
            #
            # Broadcast molecular position from processor 0 to all other processors 
            #
            R_shift = p.bcast(R_shift_o)
            #
            # Check molecules do not overlap
            #
            overlap = 0
            if( len(ELN_sys) > 0 ):
                for atom_i in myChunk:
                    # Loop over all the atoms of the molecule being placed
                    r_i = R_shift[atom_i]
                    for sys_atom in range(len(ELN_sys)):
                        # Loop over all the atoms of the system 
                        r_j = R_sys[sys_atom]
                        r_ij_sq = prop.sq_drij_c(r_i,r_j,LV)
                        if( r_ij_sq < cut_ij_sq ):
                            overlap = 1
            
                p.barrier() # Barrier for MPI_COMM_WORLD
                
            #
            # Reduce sum the overlap variable from all the processors
            #   if it is zero everywhere there was no overlap detected 
            #
            # overlap_sum = mpi.all_reduce(world,overlap, lambda x,y: x + y)
            overlap_sum = p.allReduceSum(overlap)
            
            if( overlap_sum ==  0 ):
                # If no overlap detected add molecule to the system 
                sys_mol_n += 1
                for atom_i in range( len(ELN_s) ):
                    ASYMB_sys.append( ASYMB_s[atom_i])
                    ELN_sys .append( ELN_s[atom_i])
                    ATYPE_sys.append( ATYPE_s[atom_i])
                    RESN_sys.append( RESN_s[atom_i])
                    RESID_sys.append( RESID_s[atom_i])
                    GTYPE_sys.append( GTYPE_s[atom_i])
                    CHARN_sys.append( CHARN_s[atom_i])
                    CHARGES_sys.append( CHARGES_s[atom_i])
                    AMASS_sys.append( AMASS_s[atom_i])
                    R_sys.append( R_shift[atom_i])
                    if(  options.ff_software == "lammps"  ):
                        ATYPE_IND_sys.append( ATYPE_IND_s[atom_i])
                    
                
                if( options.verbose ):
                    if( rank == 0  ):
                        print "      -  Molecule ",sys_mol_n," has been added to the system after ",moladd_atempts," placment attempts "

                        if( options.ptime ):
                            t_f = datetime.datetime.now()
	                    dt_sec  = t_f.second - t_i.second
	                    dt_min  = t_f.minute - t_i.minute
	                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
	                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                            print "        - with placement time ",dt_min," min ",dt_sec," seconds "
                            

            if( moladd_atempts >= options.max_mol_place ):
                # If attempts to place molecule into the system exceed max set by max_mol_place
                #   reset system and star over 
                if(  rank == 0  ):
                    if(  options.verbose ):
                        print "        -  Attempts to add molecule ",sys_mol_n," has exceeded max attempts ",options.max_mol_place," system will be reset for the ",sys_attempts," time "
                    
                sys_mol_n = 0                
                sys_attempts += 1
                
                ASYMB_sys = []
                ATYPE_IND_sys = []
                ELN_sys  = []
                ATYPE_sys = []
                RESN_sys = []
                RESID_sys = []
                GTYPE_sys = []
                CHARN_sys = []
                CHARGES_sys = []
                AMASS_sys = []
                R_sys = []
                VEL_sys = []

                # Add saved configuration back to system 
                for atom_p in range( len(ELN_p) ):
                    ASYMB_sys.append( ASYMB_p[atom_p])
                    ELN_sys.append( ELN_p[atom_p])
                    ATYPE_sys.append( ATYPE_p[atom_p])
                    RESN_sys.append( RESN_p[atom_p])
                    RESID_sys.append( RESID_p[atom_p])
                    GTYPE_sys.append( GTYPE_p[atom_p])
                    CHARN_sys.append( CHARN_p[atom_p])
                    CHARGES_sys.append( CHARGES_p[atom_p])
                    AMASS_sys.append( AMASS_p[atom_p])
                    R_sys.append( R_p[atom_p])
                    if(  options.ff_software == "lammps"  ):
                        ATYPE_IND_sys.append( ATYPE_IND_p[atom_p])
                              

                if( sys_attempts >= options.max_sys  ):
                    # If the system has been reset over max_sys times expand the box size by lc_expand
                    da = LV[0,0]*options.lc_expand
                    db = LV[1,1]*options.lc_expand
                    dc = LV[2,2]*options.lc_expand
                    LV[1,1] = LV[1,1] + da
                    LV[1,1] = LV[1,1] + db
                    LV[2,2] = LV[2,2] + dc
                    if( options.verbose ):

                                            
                        if( rank == 0  ):
                            print '          - Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                            print '          - Lattice vectors will be expanded by (option lc_expand)'
                            print '              -  length added ',da,db,dc
                            
            p.barrier() # Barrier for MPI_COMM_WORLD
                    
                                    
        if( sys_mol_n ==  (options.mult_s - 1) ):
            # If all the molecule have been added exit while loop and print system 
            add_mol = 0
            p.barrier() # Barrier for MPI_COMM_WORLD
            if( options.verbose and rank == 0  ):
                print " All molecules have been added "

    #
    # Start adding solvent to the system
    #
    add_mol = 1
    while ( add_mol and options.lat_s ):
        #
        # Initialize 
        #
        add_mol = 1
        overlap_sum = 1
        moladd_atempts = 0

        
        if( options.verbose):
            log_line = "  Solvent molecle is %f Angstroms long between atoms %d %d  "%(mol_l,max_dr_i,max_dr_j)
            print log_line
            log_line = "    molecules %d per side of a cube box to give %d total molecules in the volume  "%(mol_box_side,options.mult_s)
            print log_line
            log_line = "    Supercell length  %f "%(solv_box_l)
            print log_line

        # Record intial time for pereformance testing 
        if( options.ptime ):t_i = datetime.datetime.now()

        # loop over lattice positions
        #  Add 
        sol_cnt = 0 
        for x_pos in range(0,int(mol_box_side)):
	    if( sol_cnt ==  (options.mult_s -1)  ): break
            for y_pos in range(0,int(mol_box_side)):
		if( sol_cnt ==  (options.mult_s -1)  ): break
                for z_pos in range(0,int(mol_box_side)):
                    lat_pos = numpy.array( [float(x_pos)*mol_l_lj,float(y_pos)*mol_l_lj,float(z_pos)*mol_l_lj] )
                    # Make sure there is no overlap with the added molecules
                    overlap = 0
                    for atom_p in range( len(ELN_p) ):
                        r_p = R_p[atom_p]
                        r_ij_sq = prop.sq_drij_c(r_p,lat_pos,LV)

                        if( r_ij_sq < sq_mol_l_lj ):
                            overlap = 1
                            
                    p.barrier() # Barrier for MPI_COMM_WORLD
		    # Shift molecule to lattice point 
		    R_shift = prop.shift_r(R_moli_c,lat_pos)
		    
                    #
                    # Reduce sum the overlap variable from all the processors
                    #   if it is zero everywhere there was no overlap detected 
                    #
                    # overlap_sum = mpi.all_reduce(world,overlap, lambda x,y: x + y)
                    overlap_sum = p.allReduceSum(overlap)

                    if( overlap_sum ==  0 ):
			# Shift molecule to lattice point 
                        # If no overlap detected add molecule to the system 
                        sol_cnt += 1
                        for atom_i in range( len(ELN_s) ):
                            ASYMB_sys.append( ASYMB_s[atom_i])
                            ELN_sys .append( ELN_s[atom_i])
                            ATYPE_sys.append( ATYPE_s[atom_i])
                            RESN_sys.append( RESN_s[atom_i])
                            RESID_sys.append( RESID_s[atom_i])
                            GTYPE_sys.append( GTYPE_s[atom_i])
                            CHARN_sys.append( CHARN_s[atom_i])
                            CHARGES_sys.append( CHARGES_s[atom_i])
                            AMASS_sys.append( AMASS_s[atom_i])
                            R_sys.append( R_shift[atom_i])
                            if(  options.ff_software == "lammps"  ):
                                ATYPE_IND_sys.append( ATYPE_IND_s[atom_i])

                        if( options.verbose ):

                            if( rank == 0  ):
                                print "      -  Solvent molecule ",sol_cnt," has been added to the system after ",moladd_atempts," placment attempts "

                                if( options.ptime ):
                                    t_f = datetime.datetime.now()
                                    dt_sec  = t_f.second - t_i.second
                                    dt_min  = t_f.minute - t_i.minute
                                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                    print "        - with placement time ",dt_min," min ",dt_sec," seconds "
                    if( sol_cnt == (options.mult_s -1) ): break
                        

        if( sol_cnt >  (options.mult_s -1)  ):
            # Failed to fit all the needed solvent molecules in the box
            # Rest system and increase buffer
            
            if( rank == 0  ):
                if( options.verbose ):
                    log_line = "        -  %d solvent molecules add to box with length %f "%(sol_cnt,mol_box_side)
                    print log_line
                    
            sys_mol_n = 0                
            sys_attempts += 1

            ASYMB_sys = []
            ATYPE_IND_sys = []
            ELN_sys  = []
            ATYPE_sys = []
            RESN_sys = []
            RESID_sys = []
            GTYPE_sys = []
            CHARN_sys = []
            CHARGES_sys = []
            AMASS_sys = []
            R_sys = []
            VEL_sys = []

            # Add saved configuration back to system 
            for atom_p in range( len(ELN_p) ):
                ASYMB_sys.append( ASYMB_p[atom_p])
                ELN_sys.append( ELN_p[atom_p])
                ATYPE_sys.append( ATYPE_p[atom_p])
                RESN_sys.append( RESN_p[atom_p])
                RESID_sys.append( RESID_p[atom_p])
                GTYPE_sys.append( GTYPE_p[atom_p])
                CHARN_sys.append( CHARN_p[atom_p])
                CHARGES_sys.append( CHARGES_p[atom_p])
                AMASS_sys.append( AMASS_p[atom_p])
                R_sys.append( R_p[atom_p])
                if(  options.ff_software == "lammps"  ):
                    ATYPE_IND_sys.append( ATYPE_IND_p[atom_p])

            # expand the buffer distance between solvent molecules 
            options.buf_s = options.buf_s + options.buf_s*options.lc_expand
            
            p.barrier() # Barrier for MPI_COMM_WORLD
                    
                                    
        else:
            # If all the molecule have been added exit while loop and print system 
            add_mol = 0
            p.barrier() # Barrier for MPI_COMM_WORLD
            if( options.verbose and rank == 0  ):
                print " All molecules have been added "
                
    if( rank == 0 ):
        # Write xyz file 
        out_xyz = "mol_system.xyz"
        xmol.write_xyz(ASYMB_sys,R_sys,out_xyz)
        
        if( options.ff_software == "gromacs" ):
            # If gromacs write gro file
            #   assume the top file can be modified by the user 
            gromacs.print_gro(options.out_gro,GTYPE_sys,RESID_sys,RESN_sys,R_sys,LV)

        elif(  options.ff_software == "lammps"  ):
            # Find topology for entire system
            #   Repeat molecular topology for all molecules added to the system
	    mol_mult = options.mult_i + options.mult_s
            BONDS_sys,BTYPE_IND_sys = top.replicate_bonds(mol_mult,ELN_i,BONDS_i,BTYPE_IND_i) 
            ANGLES_sys,ANGTYPE_IND_sys = top.replicate_angles(mol_mult,ELN_i,ANGLES_i,ANGTYPE_IND_i) 
            DIH_sys,DTYPE_IND_sys = top.replicate_dih(mol_mult,ELN_i,DIH_i,DTYPE_IND_i) 

            debug = 0
            if( debug ):
                print " BONDS ", len(BONDS_sys)
                print "     ",BONDS_sys[0]
                for bond_indx in range(len(BONDS_sys)):
                    print BONDS_sys[bond_indx][0]+1, BONDS_sys[bond_indx][1]+1
                    
                print " ANGLES_sys ", len(ANGLES_sys)
                print "     ",ANGLES_sys[0]
                print " DIH_sys ", len(DIH_sys)
                print "     ",DIH_sys[0]
                
                sys.exit("debug 6 ")

            data_file = "mol_system.data" 
            lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
              BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
              ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
              DIH_sys,DTYPE_IND_sys,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
              RESN_sys,ATYPE_IND_sys,CHARGES_sys,R_sys , ATYPE_sys,
              BONDS_sys ,BTYPE_IND_sys, ANGLES_sys ,ANGTYPE_IND_sys, LV)

        else:
            print " Unknown ff software ",ff_software.options
        

    
if __name__=="__main__":
    main()
   
