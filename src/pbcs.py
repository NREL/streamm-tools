#! /usr/bin/env python
"""
Vector opperations for particles in a box 
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/07/2014
# travis.kemper@nrel.gov

import numpy as np

def delta_r_c(r_i,r_j,latticevec):

    # Find magnitude of dr_ij
    r_ij  = r_j - r_i
    
    r_x = r_ij[0] - latticevec[0][0] * round( r_ij[0]/  latticevec[0][0] )
    r_y = r_ij[1] - latticevec[1][1] * round( r_ij[1]/  latticevec[1][1] )
    r_z = r_ij[2] - latticevec[2][2] * round( r_ij[2]/  latticevec[2][2] )
    
    dr_pbc = np.array( [r_x,r_y,r_z] )

    return dr_pbc

def sq_drij_c(r_i,r_j,latticevec):
    """
    Reutrn square magnitude of distance between to vectors  using cubic periodic boundry conditions 
    """

    dr_pbc = delta_r_c(r_i,r_j,latticevec)
    
    sq_dr = np.dot( dr_pbc,dr_pbc)
    
    return sq_dr
    
    

def norm_r_ij(r_i,r_j,latticevec):
    """
    Normailze difference between two vectors 
    """

    debug = False
    
    delta_ij = delta_r_c(r_i,r_j,latticevec)

    if( debug):
        print "delta_ij ",delta_ij
        print " mag ",np.linalg.norm(delta_ij)
    
    return (delta_ij)/np.linalg.norm(delta_ij)

def getAngle(r_i,r_j):
    """
    Calcuate angle
      k - i - j 
      r_i = r_ik
      r_j = r_ij
      cos( \theta ) = ( a dot b ) / ( |a| |b| )
    """
    #

    r_i_norm = r_i/np.linalg.norm(r_i,)
    r_j_norm = r_j/np.linalg.norm(r_j) 
    dot_ij = np.dot(r_i_norm,r_j_norm)
    
    if( dot_ij >= 1.0 ):
       ang_deg = 0.0
    elif(  dot_ij <= -1.0 ):
       ang_deg = 180.0
    else:    
        cos_ang = np.arccos(dot_ij )
        ang_deg = np.rad2deg( cos_ang )
    
    
    return ang_deg


def replicate(p,options,oligo_array,sol_array): 

    """
    Replicate structures

    Arguments
        p (object) mpirNREL 

        options
            verbose
            ptime
            output_id
            dir_id
            json
            top
            gro
            calc_overlap 

            sol_json
            sol_top
            sol_gro

            sol_buf
            atomic_cut
            den_target
            atoms_target
            max_mol_place
            max_sys
            lc_expand
            perc_sol
            itp_file
            norm_dihparam

        oligo_array  (StructureContainer) of oligomer  molecules
        sol_array    (StructureContainer) of solvent molecules  
        
    Returns: None

    """
    import mpiNREL
    import file_io, units 
    from structureContainer import StructureContainer
    import sys , datetime, random, math 

    debug = False 

    # If options set for fixed then set random seed (for regression testing)
    if options.fixed_rnd_seed:
        random.seed(0)

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    # Open log file 
    if( rank == 0  ):
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w")

        log_lines = " Running on %d processors  \n"%(size)
        print log_lines
        log_out.write(log_lines)

        t_i = datetime.datetime.now()

    # Set location of origin
    org = np.array( [0.0,0.0,0.0] )

    oligo_cnt     = 0      # Number of oligomers to replicate 
    oligo_nprt    = 0      # Total number of particles in all the oligomers 
    oligo_mass    = 0.0    # Total mass of particles in all the oligomers
    oligo_maxlength = 0.0 
    for oligo_i in oligo_array:
        oligo_cnt += 1
        nprt = len( oligo_i.ptclC )

        print "nprt",nprt

        oligo_nprt += nprt 
        tot_mass = oligo_i.getTotMass()
        oligo_mass += tot_mass
        nbonds = len( oligo_i.bondC )
        # Shift center of mass to origin
        oligo_i.shift_center_mass(org)

        # Record oligomer molecule length
        oligo_length = oligo_i.getlength()

        # Record solvent molecule length for grid spacing 
        oligo_length = oligo_i.getlength()
        if( oligo_length > oligo_maxlength):
            oligo_maxlength = oligo_length

        if( options.verbose and rank == 0  ):
            log_lines = ""
            log_lines += "  Oligomers %d \n"%oligo_cnt 
            log_lines +=  "    Particles  %d \n"%nprt
            log_lines +=  "    Total mass   %f \n"%tot_mass
            log_lines +=  "    Length   %f \n"%oligo_length
            log_lines +=  "    Bonds  %d \n"%nbonds
            print log_lines
            log_out.write(log_lines)


    sol_cnt     = 0   # Number of solvent to replicate 
    sol_nprt    = 0   # Total number of solvent in all the solture 
    sol_mass    = 0.0    # Total mass of solvent in all the solture
    sol_maxlength = -100000.0 
    #

    for sol_i in sol_array:
        sol_cnt += 1
        nprt = len( sol_i.ptclC )
        sol_nprt += nprt 
        tot_mass = sol_i.getTotMass()
        sol_mass += tot_mass
        nbonds = len( sol_i.bondC )

        # Shift center of mass to origin
        sol_i.shift_center_mass(org)

        # Record solvent molecule length for grid spacing 
        sol_length = sol_i.getlength()
        if( sol_length > sol_maxlength):
            sol_maxlength = sol_length

        if( options.verbose and rank == 0  ):
            log_lines = ""
            log_lines += "  Solvents %d \n"%sol_cnt 
            log_lines +=  "    Particles  %d \n"%nprt
            log_lines +=  "    Total mass   %f \n"%tot_mass
            log_lines +=  "    Length   %f \n"%sol_length
            log_lines +=  "    Bonds  %d \n"%nbonds
            print log_lines
            log_out.write(log_lines)
    #
    # Calculate the number of oligomers and  solvent molecules
    #
    if( options.perc_sol > 0.0 ):
        # Variables
        #   atoms_target - target number of atoms # options.atoms_target
        #   perc_sol - perecent solvent by mass
        #   frac_sol - fraction solvent by mass
        #   sol_nprt - number atoms in the list of solvent molecules 
        #   oligo_nprt - number atoms in the list of oligomer molecules 
        #   n_sol_l - number of solvent list replications
        #   n_olgio_l - number of oligomer list replications
        #   sol_mass - mass of all the solvents in the solvent list
        #   oligo_mass - mass of all the oligomers in the oligomer list
        # Equations
        #   Equ 1 : perc_sol = n_sol_l*sol_mass/( n_sol_l*sol_mass + n_olgio_l*oligo_mass )
        #   Equ 2 : atoms_target =  n_sol_l*sol_nprt + n_olgio_l*oligo_nprt
        # Solutions
        frac_sol = options.perc_sol/100.0
        n_sol_l= int((-1.0*frac_sol*oligo_mass*float(options.atoms_target))/(frac_sol*sol_mass*float(oligo_nprt)  - frac_sol*oligo_mass*float(sol_nprt) - sol_mass*float(oligo_nprt) ))
        n_olgio_l = int( (float(options.atoms_target) - float(sol_nprt)*float(n_sol_l))/float(oligo_nprt))
    else:
        n_sol_l = 0
        n_olgio_l = int(float(options.atoms_target)/float(oligo_nprt))
        solv_box_l = 0.0


        print " atoms_target oligo_bnprt ",float(options.atoms_target),float(oligo_nprt)

    #
    # Calculate the box size for a target density 
    #
    target_density_amuang = units.convert_gcm3_AMUA3( options.den_target) # densit in AMU/Angstrom^3
    
    #
    if( n_olgio_l + n_sol_l <= 0 ):
        print "  n_olgio_l ",n_olgio_l
        print "  n_sol_l ",n_sol_l
        print " Specified number of target atoms %d produced no replications "%(options.atoms_target)
        sys.exit("error")

    total_n = n_olgio_l*oligo_nprt + n_sol_l*sol_nprt
    total_mass = oligo_mass*n_olgio_l + n_sol_l*sol_mass
    volume_target_ang = total_mass/target_density_amuang
    len_target_ang = volume_target_ang**(1.0/3.0)

    if( options.perc_sol > 0.0 ):
        # Check to be sure the grid is large enough
        sol_box_side = int(math.ceil(n_sol_l**(1.0/3.0) ) )         # Number of solvents per box side 
        sol_length = sol_maxlength + options.sol_buf           # length of solvent 
        sol_length_sq = sol_length*sol_length                  # length of solvent squared for overlap calculation 
        vol_olgio = float(n_olgio_l)*(oligo_maxlength**3.0)    # Volume occupied by oligomer 
        len_sol_box =  float(sol_box_side)*sol_length              # Length of pure solvent box
        vol_sol = len_sol_box**3.0                                 # Volume occupied by solvent 
        tot_mol_vol = vol_olgio + vol_sol                      # total volume for oligomers and solvents
        len_tot_box = tot_mol_vol**(1.0/3.0)

        if( len_target_ang < len_tot_box ):
            if( options.verbose and rank == 0 ):
                log_lines = "  Expanding box length from %f to %f based on the length of the solvent  "%(len_target_ang,len_tot_box)
                print log_lines
                log_out.write(log_lines)

            len_target_ang = len_tot_box
            volume_target_ang = len_tot_box**3.0

        sol_grid_len = len_target_ang/float(sol_box_side)

        print "len_target_ang ",len_target_ang 
        print "sol_grid_len ",sol_grid_len 

    else:
        len_target_ang = volume_target_ang**(1.0/3.0)

    cut_ij_sq = options.atomic_cut* options.atomic_cut
    # Calculate actual final structure properties
    vol_f = len_target_ang**3.0
    den_AMU_f = total_mass/vol_f
    den_f = units.convert_AMUA3_gcm3(den_AMU_f)
    perc_sol_f = (n_sol_l*sol_mass)/(n_sol_l*sol_mass + n_olgio_l*oligo_mass)

    # Recalculate solvent molecules along box length 

    print "target_density_amuang",target_density_amuang, options.den_target
    print "total_n",total_n
    print "total_mass",total_mass
    print "len_target_ang",len_target_ang
    print n_sol_l,sol_mass,n_olgio_l,oligo_mass

    # oligomer molecule container
    oligomer_rep = StructureContainer() 
    # Set lattice vector to new box size
    latvec_list = [np.array([len_target_ang,0.0,0.0]),np.array( [0.0,len_target_ang,0.0]),np.array( [0.0,0.0,len_target_ang]) ]

    if(debug):
        print "latvec_list" , latvec_list
    oligomer_rep.setLatVec(latvec_list)

    if(debug):
        print " oligomer_rep.getLatVec() ",oligomer_rep.getLatVec()


    strucC = StructureContainer()

    # Solvent molecule container
    sol_rep = StructureContainer()  
    sol_rep.setLatVec(latvec_list)
    strucC.setLatVec(latvec_list)

    if(debug):
        print " s1 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s1 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s1 strucC. .getLatVec() ",strucC.getLatVec()


    # Print script information and settings 
    if( rank == 0  ):

        log_lines =  "  Replication settings \n"
        log_lines += "   - Tragets \n"
        log_lines += "       Total atoms %d \n"%(options.atoms_target)
        log_lines += "       Density %f g/cm^3 %f AMU/Angstrom^3 \n"%(options.den_target,target_density_amuang)
        log_lines += "       Solvent mass percentage %f \n"%(options.perc_sol)
        log_lines += "   - Oligomers \n"
        log_lines += "       Total atoms in set %d \n"%(oligo_nprt)
        log_lines += "       Set of structures will replicated %d times \n"%(n_olgio_l)
        if( options.perc_sol > 0.0 ):
            log_lines += "   - Solvents \n"
            log_lines += "       Total atoms in set %d \n"%(sol_nprt)
            log_lines += "       With a length of  %f Angstrom \n"%(sol_length)
            log_lines += "       Set of structures will replicated %d times \n"%(n_sol_l)
            log_lines += "       With an initial grid spacing of %f  \n"%(sol_grid_len)

        log_lines += "   - Final porperties  \n"
        log_lines += "       Total atoms %d \n"%(total_n)
        log_lines += "       Volume %f Angstrom^3 \n"%(volume_target_ang)
        log_lines += "       Density %f g/cm^3 %f AMU/Angstrom^3 \n"%(den_f,den_AMU_f)
        log_lines += "       Solvent mass percentage %f \n"%(perc_sol_f)
        log_lines += "       Cubic cell with length of %f Angstroms \n"%(len_target_ang)
        log_lines += "   - Placement \n"
        log_lines += "       Maximum structure placements %d \n"%(options.max_mol_place)
        log_lines += "       Maximum number of restarts before box is expanded %d \n"%(options.max_sys)
        log_lines += "       Percent of box size to add during expantion %8.2f \n"%(100*options.lc_expand)
        log_lines += "   - id's \n"
        log_lines += "       Directory %s \n"%(options.dir_id)	
        log_lines += "       Output id %s \n"%(options.output_id)
        print log_lines
        log_out.write(log_lines)

    p.barrier()


    #sys.exit(" debug 2 ")

    # Record initial time
    if( rank == 0  ): 
        t_i = datetime.datetime.now()

    sys_oligo_n = 0     # Number of oligomers add to the system
    sys_attempts = 0    # Number of times the system has been reset
    struc_add_cnt = 0   # Total number of structures added to the final structure

    max_oligo_residue_number_list = [] # Max residue number for relabeling
    for struc_i in oligo_array:
        max_i = 0 
        for pid, ptclObj in struc_i.ptclC :
            if( ptclObj.tagsDict["residue"] > max_i ):
                max_i  = ptclObj.tagsDict["residue"]
        max_oligo_residue_number_list.append(max_i)
        
    #
    # Start adding molecules to the system
    #
    add_oligo = True

    print " Adding %d  oligomers  "%n_olgio_l
    while ( add_oligo ):
        #
        # Initialize 
        #
        add_oligo = True
        overlap_found = True
        strucadd_atempts = 0


        # Record intial time for pereformance testing 
        if( rank == 0  ): 
            tadd_i = datetime.datetime.now()
        for oligo_l in range( n_olgio_l ):
            # loop over the number of times each oligomer in the oligomer list needs to be replicated
            o_cnt = -1
            for struc_i in oligo_array:
                o_cnt += 1 
                max_oligo_residue_number = max_oligo_residue_number_list[o_cnt]
                if(debug):
                    print " s12 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
                    print " s12 sol_rep.getLatVec() ",sol_rep.getLatVec()
                    print " s12 strucC. .getLatVec() ",strucC.getLatVec()


                #
                # Place the atomic indices into list 
                # 
                pointIndices = range( len(struc_i.ptclC)  )
                if( debug ):
                    print rank, size," splitOnProcs "
                # Create a list of atomic indices for each processor 
                myChunk  = p.splitListOnProcs(pointIndices)
                p_debug = False 
                if(p_debug):                
                    print " cpu ",rank ," has atoms ",myChunk[0]," - ",myChunk[len(myChunk)-1],"  \n"

                    sys.exit("P debug 1 ")

                # For each structure add to 
                while ( overlap_found ):
                    strucadd_atempts += 1

                    n_dim = 3
                    rot_angle_i_o = 0.0 
                    rot_angle_j_o = 0.0 
                    r_random_o  = np.zeros(n_dim)

                    if ( rank == 0 ):
                        #
                        #  Get random rotation angles from single processor  
                        #
                        ang_acc = 1000  # number of digets in random angle 
                        rot_angle_i_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                        rot_angle_j_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                        #
                        #  Get random translation from single processor 
                        #
                        r_random_o = np.zeros(n_dim)
                        for x_indx in range( n_dim ):
                            r_random_o[x_indx] = random.randrange(0,ang_acc*int(oligomer_rep.latvec[x_indx][x_indx]) )/float(ang_acc)

                        debug = 0
                        if( debug ):
                            print " ran ",x_indx,(oligomer_rep.latvec[x_indx][x_indx])
                            print rot_angle_i_o,rot_angle_j_o,r_random_o
                            #sys.exit(" Random # 's test 1")

                    p.barrier() # Barrier for MPI_COMM_WORLD
                    #
                    # Broadcast random rotation angles and translations to all processors 
                    #
                    rot_angle_i = p.bcast(rot_angle_i_o)
                    rot_angle_j = p.bcast(rot_angle_j_o)
                    r_random = p.bcast(r_random_o)
                    p.barrier() # Barrier for MPI_COMM_WORLD
                    #
                    # Get coordinates of randomly rotated and shifted 
                    #
                    struc_i.shift_center_mass(org)
                    struc_i.rotate(rot_angle_i,rot_angle_j)
                    struc_i.vec_shift(r_random)

                    overlap = 0
                    if( len(oligomer_rep.ptclC) > 0 ):
                        #
                        # If there are particles in the system check atoms do not overlap
                        #   
                        if( options.calc_overlap ):

                            for p_i, ptclObj_i in struc_i.ptclC(myChunk):
                                r_i = np.array( ptclObj_i.position )
                                for p_sys, ptclObj_sys in oligomer_rep.ptclC :
                                    r_sys = np.array( ptclObj_sys.position )
                                    r_ij_sq = sq_drij_c(r_i,r_sys,oligomer_rep.getLatVec() )
                                    if( r_ij_sq < cut_ij_sq ):
                                        overlap = 1

                    p.barrier() # Barrier for MPI_COMM_WORLD
                    #
                    # Reduce sum the overlap variable from all the processors
                    #   if it is zero everywhere there was no overlap detected 
                    #
                    overlap_sum = p.allReduceSum(overlap)
                    p.barrier() # Barrier for MPI_COMM_WORLD

                    if(debug):
                        print " s13 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
                        print " s13 sol_rep.getLatVec() ",sol_rep.getLatVec()
                        print " s13 strucC. .getLatVec() ",strucC.getLatVec()


                    if( overlap_sum ==  0 ):
                        # If no overlap detected add molecule to the system 
                        sys_oligo_n += 1
                        struc_add_cnt += 1
                        # Rest molecule numbers
                        for pid, ptclObj in struc_i.ptclC :
                            ptclObj.tagsDict["chain"] = struc_add_cnt
                            if( sys_oligo_n > 1 ):
                                res_numb_i = max_oligo_residue_number + ptclObj.tagsDict["residue"]
                                ptclObj.tagsDict["residue"] = res_numb_i

                            

                        # add oligomer structure to system structure
                        struc_i.setLatVec(oligomer_rep.getLatVec())
                        oligomer_rep += struc_i

                        if( options.verbose ):
                            if( rank == 0  ):
                                print "      -  Molecule ",sys_oligo_n," has been added to the system after ",strucadd_atempts," placment attempts "
                                print "         system has %d atoms and %d bonds "%(len(oligomer_rep.ptclC),len(oligomer_rep.bondC))
                                #print " Printing  oligomer_rep bonds "
                                #oligomer_rep.printbondlengths()

                                if(debug):
                                    print " s11 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
                                    print " s11 sol_rep.getLatVec() ",sol_rep.getLatVec()
                                    print " s11 strucC. .getLatVec() ",strucC.getLatVec()


                                if( options.ptime ):
                                    t_f = datetime.datetime.now()
                                    dt_sec  = t_f.second - t_i.second
                                    dt_min  = t_f.minute - t_i.minute
                                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                    print "        - with placement time ",dt_min," min ",dt_sec," seconds "

                        overlap_found = False
                    else:
                        overlap_found = True

                    if( strucadd_atempts >= options.max_mol_place ):
                        # If attempts to place molecule into the system exceed max set by max_mol_place

                        #   reset system and star over 
                        if(  rank == 0  ):
                            if(  options.verbose ):


                                print " add0 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()


                                print "        -  Attempts to add molecule ",sys_oligo_n+1," has exceeded max attempts ",options.max_mol_place," system will be reset for the ",sys_attempts," time "

                        sys_oligo_n = 0
                        struc_add_cnt = 0 
                        strucadd_atempts = 0
                        sys_attempts += 1

                        # Save lattice vectors as to no loose any expansions 
                        latvec_i = oligomer_rep.getLatVec()

                        print " saving s1 latvec_i ",latvec_i

                        # Delete system 
                        del oligomer_rep
                        oligomer_rep = StructureContainer()  # Output replicated structure
                        # Set lattice vectors 
                        oligomer_rep.setLatVec(latvec_i) 


                    if( sys_attempts >= options.max_sys  ):


                        print " exp0 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()


                        # If the system has been reset over max_sys times expand the box size by lc_expand
                        oligomer_rep.expandLatVec(options.lc_expand)

                        # Save lattice vectors as to no loose any expansions 
                        latvec_i = oligomer_rep.getLatVec()


                        print " exp1 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()

                        # Delete system 
                        del oligomer_rep
                        oligomer_rep = StructureContainer()  # Output replicated structure
                        # Set lattice vectors 
                        oligomer_rep.setLatVec(latvec_i) 

                        sys_attempts = 0

                        if( options.verbose ):                
                            if( rank == 0  ):
                                print '          - Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                                print '          - Lattice vectors will be expanded by (option lc_expand)',options.lc_expand
                                print '             v_1 ',latvec_i[0]
                                print '             v_2 ',latvec_i[1]
                                print '             v_3 ',latvec_i[2]

            p.barrier() # Barrier for MPI_COMM_WORLD


        if( sys_oligo_n ==  n_olgio_l  ):
            # If all the molecule have been added exit while loop and print system 
            add_oligo = False
            latvec_oligo = oligomer_rep.getLatVec()
            p.barrier() # Barrier for MPI_COMM_WORLD
            if( options.verbose and rank == 0  ):
                print " All oligomers  have been added "

    # Add replicated oligmers to final structure
    # strucC = StructureContainer()

    debgu_n = False
    if( debgu_n ):

        print "strucC pre add "
        print strucC.ptclC
        print strucC.bondC

        print "oligomer_rep"
        print oligomer_rep.ptclC
        print oligomer_rep.bondC


    if(debug):

        print " s21 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s21 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s21 strucC. .getLatVec() ",strucC.getLatVec()


    strucC =  oligomer_rep
    strucC.setLatVec(latvec_oligo)
    # Rest solvent lattice vectors to match oligo 
    sol_rep.setLatVec(latvec_oligo)

    if( debgu_n ):
        print "strucC"
        print strucC.ptclC
        print strucC.bondC

        sys.exit("debug ")

    if(debug):

        print " s2 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s2 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s2 strucC. .getLatVec() ",strucC.getLatVec()


    if( options.perc_sol > 0.0 ):

        sys_sol_n = 0    # Number of solvent add to the system
        sys_attempts = 0  # Number of times the system has been reset
        #
        # Start adding molecules to the system
        #
        print " Adding %d  solvent  "%n_sol_l
        add_sol = True
        while ( add_sol ):
            #
            # Initialize 
            #
            strucadd_atempts = 0

            # Record intial time for pereformance testing 
            if( rank == 0  ): 
                tadd_i = datetime.datetime.now()

            bad_grid_point = False 
            for sol_l in range( n_sol_l ):

                #
                # 
                #
                print " starting %d of %d "%(sol_l,n_sol_l)

                # loop over the number of times each solvent in the solvent list needs to be replicated
                if( bad_grid_point ): break
                if( not add_sol): break 
                for struc_i in sol_array:
                    overlap_found = True
                    if( bad_grid_point ): break
                    if( not add_sol): break 
                    while ( overlap_found ):
                        strucadd_atempts += 1
                        for x_indx in range(sol_box_side):
                            if( not add_sol): break 
                            for y_indx in range(sol_box_side):
                                if( not add_sol): break 
                                for z_indx in range(sol_box_side):

                                    if( sys_sol_n == n_sol_l*sol_cnt ):
                                        add_sol = False
                                        break


                                    #l_x =  float(lat_indx[0])*sol_length
                                    #l_y =  float(lat_indx[1])*sol_length
                                    #l_z =  float(lat_indx[2])*sol_length

                                    l_x =  float(x_indx)*sol_grid_len
                                    l_y =  float(y_indx)*sol_grid_len
                                    l_z =  float(z_indx)*sol_grid_len

                                    print " Checking overlap for solvent %d at lattice point %f %f %f "%(sys_sol_n,l_x,l_y,l_z)

                                    if( l_x > oligomer_rep.latvec[0][0] or l_y > oligomer_rep.latvec[1][1] or l_z > oligomer_rep.latvec[2][2] ):
                                        print " Lattic point beyond box %f %f %f "%(oligomer_rep.latvec[0][0],oligomer_rep.latvec[1][1], oligomer_rep.latvec[2][2])
                                        bad_grid_point = True
                                        break 

                                    lat_pos = np.array( [l_x,l_y,l_z] )

                                    # Make sure there is no overlap with the added molecules
                                    overlap = 0
                                    for p_i, ptclObj_i in oligomer_rep.ptclC :
                                        r_i = np.array( ptclObj_i.position )
                                        r_ij_sq = sq_drij_c(r_i,lat_pos,oligomer_rep.getLatVec() )
                                        if( r_ij_sq < sol_length_sq ):
                                            overlap = 1

                                    p.barrier() # Barrier for MPI_COMM_WORLD
                                    #
                                    # Reduce sum the overlap variable from all the processors
                                    #   if it is zero everywhere there was no overlap detected 
                                    #
                                    overlap_sum = p.allReduceSum(overlap)

                                    if( overlap_sum ==  0 ):
                                        # If no overlap detected add molecule to the system 
                                        sys_sol_n += 1
                                        overlap_found = False 
                                        # Shift molecule to lattice point
                                        struc_i.shift_center_mass(org)
                                        struc_i.vec_shift(lat_pos)

                                        struc_i.setLatVec(sol_rep.getLatVec())
                                        sol_rep += struc_i

                                        if( options.verbose ):
                                            if( rank == 0  ):
                                                print "      -  Molecule %d  has been added to the system at lattice point %f %f %f  "%(sys_sol_n,l_x,l_y,l_z)

                                                if( options.ptime ):
                                                    t_f = datetime.datetime.now()
                                                    dt_sec  = t_f.second - t_i.second
                                                    dt_min  = t_f.minute - t_i.minute
                                                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                                    print "        - with placement time ",dt_min," min ",dt_sec," seconds "
                                    else:
                                        overlap_found = True

                                        print "  lattice point was found to overlap "



            if( sys_sol_n == n_sol_l*sol_cnt ):
                add_sol = False
                p.barrier() # Barrier for MPI_COMM_WORLD
                if( options.verbose and rank == 0  ):
                    print " All solvents  have been added "

            else:
                #
                # If attempts to place solvent molecule into the system failed
                #
                sys_sol_n = 0                
                sys_attempts += 1
                #
                # Save lattice vectors as to not loose any expansions
                #
                latvec_i = sol_rep.getLatVec()
                #
                # Delete system
                #
                del sol_rep
                sol_rep = StructureContainer()  # Output replicated structure
                #
                # Set lattice vectors
                #
                sol_rep.setLatVec(latvec_i)
                #
                #
                if( (sol_grid_len*0.90) > sol_length):
                    # If grid spacing is larger than the solvent length shrink grid
                    sol_grid_len = sol_grid_len*0.90
                else:
                    # Otherwise increase volume and set grid spacing to solvent length
                    sol_grid_len = sol_length
                    #
                    # Expand the box size by a single solvent length 
                    #
                    #sol_rep.expandLatVec(options.lc_expand)
                    sol_rep.latvec[0] = sol_rep.latvec[0] + sol_length
                    sol_rep.latvec[1] = sol_rep.latvec[1] + sol_length
                    sol_rep.latvec[2] = sol_rep.latvec[2] + sol_length
                    #
                    # Increase the number of solvent molecules along the box by 1
                    #
                    sol_box_side += 1


                if( options.verbose ):                
                    if( rank == 0  ):
                        print '          - Lattice vectors will be expanded by solvent length %f ',sol_length

            p.barrier() # Barrier for MPI_COMM_WORLD

        # Add replicated solvents to final structure
        strucC += sol_rep
        strucC.setLatVec(sol_rep.latvec)


    if(debug):
        print " s3 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s3 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s3 strucC. .getLatVec() ",strucC.getLatVec()


    strucC.compressPtclIDs()

    print "         f_rep has %d atoms and %d bonds "%(len(strucC.ptclC),len(strucC.bondC))
    print "             lat vec 1 ",sol_rep.latvec[0]
    print "             lat vec 2 ",sol_rep.latvec[1]
    print "             lat vec 3 ",sol_rep.latvec[2]
    #strucC.printbondlengths()


    if( rank == 0 ):
        t_f = datetime.datetime.now()
        dt_sec  = t_f.second - t_i.second
        dt_min  = t_f.minute - t_i.minute
        if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
        log_line="\n  Finished time  " + str(t_f)
        log_out.write(log_line)
        log_line="\n  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "
        log_out.write(log_line)

        log_out.close()


    debgu_n = False 
    if( debgu_n ):
        print "strucC"
        print strucC.ptclC
        print strucC.bondC

        sys.exit("debug 2")

    return strucC
