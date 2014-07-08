#! /usr/bin/env python
"""
Run Molecular Dynamics simulation at low temperature NVT to get a minimized structure
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/02/2014
# travis.kemper@nrel.gov

from structureContainer import StructureContainer
import pbcs
import mpiNREL

import numpy as np 
import sys , datetime, random, math 

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
    parser.add_option("-o","--output_id", dest="output_id", default="replicate",type="string",help=" prefix for output files  ")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on ")
    parser.add_option("--sol_json", dest="sol_json", default="",type="string",help=" json files of solvents to add to act on ")

    # Solv
    parser.add_option("--sol_buf", dest="sol_buf",  type=float, default=3.0, help=" intra solvent buffer " )
    # Cut-off
    #   should be replaced by Van der Waals radi  
    parser.add_option("--atomic_cut", dest="atomic_cut", type=float, default=2.5, help="Minimum distance between atoms of molecules ")

    # Replication 
    parser.add_option("--den_target", dest="den_target", type=float, default=0.01, help="Target density g/cm^3 ")
    parser.add_option("--atoms_target", dest="atoms_target", type=int, default=100000, help="Target number of atoms ")
    parser.add_option("--max_mol_place", dest="max_mol_place", type=float, default=50, help="Maximum attempts to place a molecule  ")
    parser.add_option("--max_sys", dest="max_sys", type=float, default=3, help="Maximum system recreations at a certain lattice constant ")
    parser.add_option("--lc_expand", dest="lc_expand", type=float, default=0.100, help="Fraction of the box size to increase system size after max_sys is excieded ")
    parser.add_option("--perc_sol", dest="perc_sol", type=float, default=0.0, help="Percent solvent by mass ")
    
    (options, args) = parser.parse_args()
        
    return options, args

def struc_array_json(json_list):
    """
    Loop over list of json files and place each structure container into a list

    Arguments
      json_list (list) list of json files
      
    Returns
      struc_array (list) of structure objects
      
    """

    # Read in structure containers from json files 
    struc_array = []
    if( len(json_list)):
        # Loop over json files 
        json_files = json_list.split(',')
        if( len(json_files) > 0 ):
            # Read index files from args
            for json_file in json_files:
                struc_i = StructureContainer()
                struc_i.getsys_json(json_file)
                struc_array.append(struc_i)
                
    return struc_array


def main():
    """
    Replicate structures  

    Returns: None
    """

    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    options, args = get_options()

    # Open log file 
    if( rank == 0  ):
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w")
        
    # Read in oligomers  from json files 
    oligo_array = struc_array_json(options.json)

    # Read in solvents from json files 
    sol_array = struc_array_json(options.sol_json)
    
    # Set location of origin
    origin = np.array( [0.0,0.0,0.0] )
    
    oligo_cnt   = 0   # Number of oligomers to replicate 
    oligo_nprt    = 0   # Total number of particles in all the oligomers 
    oligo_mass  = 0.0    # Total mass of particles in all the oligomers 
    for oligo_i in oligo_array:
        oligo_cnt += 1
        nprt = len( oligo_i.ptclC )
        oligo_nprt += nprt 
        tot_mass = oligo_i.getTotMass()
        oligo_mass += tot_mass
        nbonds = len( oligo_i.bondC )
        # Shift center of mass to origin
        oligo_i.shift_center_mass(origin)
        
        # Record oligomer molecule length
        oligo_length = oligo_i.getlength()
        
        if( options.verbose ):
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
        sol_i.shift_center_mass(origin)

        # Record solvent molecule length for grid spacing 
        sol_length = sol_i.getlength()
        if( sol_length > sol_maxlength):
            sol_maxlength = sol_length

        if( options.verbose ):
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
        # Check to be sure the grid is large enough
        sol_length = sol_maxlength + options.sol_buf
        sol_box_side = math.ceil( n_sol_l**(1.0/3.0) )
        solv_box_l = sol_box_side*sol_length #+ mol_l_lj
        sol_length_sq = sol_length*sol_length
    else:
        n_sol_l = 0
        n_olgio_l = int(float(options.atoms_target)/float(oligo_nprt))
        solv_box_l = 0.0 

    #
    # Calculate the box size for a target density 
    #
    target_density_amuang = options.den_target*const_avo/10.0 # densit in AMU/Angstrom^3
    total_n = n_olgio_l*oligo_nprt + n_sol_l*sol_nprt
    total_mass = oligo_mass*n_olgio_l + n_sol_l*sol_mass
    volume_target_ang = total_mass/target_density_amuang    
    len_target_ang = volume_target_ang**(1.0/3.0)
    
    if( len_target_ang < solv_box_l ):
        # Extend box length if too small to fit 
        len_target_ang = solv_box_l
        
    cut_ij_sq = options.atomic_cut* options.atomic_cut
    # Calculate actual final structure properties
    vol_f = len_target_ang**3.0
    den_AMU_f = total_mass/vol_f
    den_f = den_AMU_f/const_avo*10.0
    perc_sol_f = (n_sol_l*sol_mass)/(n_sol_l*sol_mass + n_olgio_l*oligo_mass)

    print n_sol_l,sol_mass,n_olgio_l,oligo_mass
    
    # oligomer molecule container
    oligomer_rep = StructureContainer() 
    # Set lattice vector to new box size
    latvec_list = [np.array([len_target_ang,0.0,0.0]),np.array( [0.0,len_target_ang,0.0]),np.array( [0.0,0.0,len_target_ang]) ]    
    oligomer_rep.setLatVec(latvec_list)

    # Solvent molecule container
    sol_rep = StructureContainer()  
    sol_rep.setLatVec(latvec_list)

    
    
    
    # Print script information and settings 
    if( rank == 0  ):
        
        log_lines =  "  Replication settings \n"
        log_lines += "   - Tragets \n"
        log_lines += "       Total atoms %d \n"%(options.atoms_target)
        log_lines += "       Density %f g/cm^3 %f AMU/Angstrom^3 \n"%(options.den_target,target_density_amuang)
        log_lines += "       Solvent mass percentage %f \n"%(options.perc_sol)
        log_lines += "   - Oligomers \n"
        log_lines += "       Total atoms in set %d \n"%(sol_nprt)
        log_lines += "       Set of structures will replicated %d times \n"%(n_olgio_l)
        log_lines += "   - Solvents \n"
        log_lines += "       Total atoms in set %d \n"%(oligo_nprt)
        log_lines += "       Set of structures will replicated %d times \n"%(n_sol_l)
        log_lines += "   - Final porperties  \n"
        log_lines += "       Total atoms %d \n"%(total_n)
        log_lines += "       Density %f g/cm^3 %f AMU/Angstrom^3 \n"%(den_f,den_AMU_f)
        log_lines += "       Solvent mass percentage %f \n"%(perc_sol_f)
        log_lines += "       Cubic cell with length of %f Angstroms \n"%(len_target_ang)
        log_lines += "   - Placement \n"
        log_lines += "       Maximum structure placements %d \n"%(options.max_mol_place)
        log_lines += "       Maximum number of restarts before box is expanded %d \n"%(options.max_sys)
        log_lines += "       Percent of box size to add during expantion %8.2f \n"%(100*options.lc_expand)
        print log_lines
        log_out.write(log_lines)

    p.barrier()

    # Record initial time
    if( rank == 0  ): 
	t_i = datetime.datetime.now()

    sys_oligo_n = 0    # Number of oligomers add to the system
    sys_attempts = 0  # Number of times the system has been reset
    #
    # Start adding molecules to the system
    #
    add_oligo = True
    while ( add_oligo ):
        #
        # Initialize 
        #
        add_oligo = True
        overlap_found = True
        strucadd_atempts = 0

        print " Adding structure "
        
        # Record intial time for pereformance testing 
        if( rank == 0  ): 
            tadd_i = datetime.datetime.now()
        for oligo_l in range( n_olgio_l ):
            # loop over the number of times each oligomer in the oligomer list needs to be replicated
            for struc_i in oligo_array:
                # For each structure add to 
                while ( overlap_found ):
                    strucadd_atempts += 1

                    if ( rank == 0 ):
                        #
                        #  Get random rotation angles from single processor  
                        #
                        ang_acc = 1000  # number of digets in random angle 
                        rot_angle_i = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                        rot_angle_j = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                        #
                        #  Get random translation from single processor 
                        #
                        n_dim = 3
                        r_random = np.zeros(n_dim)
                        for x_indx in range( n_dim ):
                            r_random[x_indx] = random.randrange(0,ang_acc*int(oligomer_rep.latticevec[x_indx][x_indx]) )/float(ang_acc)

                        debug = 0
                        if( debug ):
                            print " ran ",x_indx,(oligomer_rep.latticevec[x_indx][x_indx])
                            print rot_angle_i,rot_angle_j,r_random
                            sys.exit(" Random # 's test 1")

                    #
                    # Broadcast random rotation angles and translations to all processors 
                    #
                    rot_angle_i = p.bcast(rot_angle_i)
                    rot_angle_j = p.bcast(rot_angle_j)
                    r_random = p.bcast(r_random)
                    #
                    # Get coordinates of randomly rotated and shifted 
                    #
                    struc_i.rotate(rot_angle_i,rot_angle_j)
                    struc_i.vec_shift(r_random)

                    overlap = 0
                    if( len(oligomer_rep.ptclC) > 0 ):
                        #
                        # If there are particles in the system check atoms do not overlap
                        #
                        for p_i, ptclObj_i in struc_i.ptclC :
                            r_i = np.array( ptclObj_i.position )
                            for p_sys, ptclObj_sys in oligomer_rep.ptclC :
                                r_sys = np.array( ptclObj_sys.position )
                                r_ij_sq = pbcs.sq_drij_c(r_i,r_sys,oligomer_rep.getLatVec() )
                                if( r_ij_sq < cut_ij_sq ):
                                    overlap = 1

                    p.barrier() # Barrier for MPI_COMM_WORLD
                    #
                    # Reduce sum the overlap variable from all the processors
                    #   if it is zero everywhere there was no overlap detected 
                    #
                    overlap_sum = p.allReduceSum(overlap)

                    if( overlap_sum ==  0 ):
                        # If no overlap detected add molecule to the system 
                        sys_oligo_n += 1
                        # add oligomer structure to system structure 
                        oligomer_rep += struc_i

                        if( options.verbose ):
                            if( rank == 0  ):
                                print "      -  Molecule ",sys_mol_n," has been added to the system after ",strucadd_atempts," placment attempts "

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
                                print "        -  Attempts to add molecule ",sys_oligo_n," has exceeded max attempts ",options.max_mol_place," system will be reset for the ",sys_attempts," time "

                        sys_oligo_n = 0                
                        sys_attempts += 1

                        # Save lattice vectors as to no loose any expansions 
                        latvec_i = oligomer_rep.getLatVec()
                        # Delete system 
                        del oligomer_rep
                        oligomer_rep = StructureContainer()  # Output replicated structure
                        # Set lattice vectors 
                        oligomer_rep.setLatVec(latvec_i) 


                    if( sys_attempts >= options.max_sys  ):
                        # If the system has been reset over max_sys times expand the box size by lc_expand
                        oligomer_rep.expandLatVec(options.lc_expand)
                    
                        if( options.verbose ):                
                            if( rank == 0  ):
                                print '          - Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                                print '          - Lattice vectors will be expanded by (option lc_expand)',options.lc_expand
                            
            p.barrier() # Barrier for MPI_COMM_WORLD
                    
                                    
        if( sys_oligo_n ==  (n_olgio_l - 1) ):
            # If all the molecule have been added exit while loop and print system 
            add_oligo = False
            p.barrier() # Barrier for MPI_COMM_WORLD
            if( options.verbose and rank == 0  ):
                print " All oligomers  have been added "


    sys_sol_n = 0    # Number of solvent add to the system
    sys_attempts = 0  # Number of times the system has been reset
    #
    # Start adding molecules to the system
    #
    add_sol = True
    while ( add_sol ):
        #
        # Initialize 
        #
        add_sol = True
        strucadd_atempts = 0

        print " Adding solvent "
        
        # Record intial time for pereformance testing 
        if( rank == 0  ): 
            tadd_i = datetime.datetime.now()

        # Initialize lattice indexs
        lat_indx = np.zeros(3)
        dim_indx = -1
        bad_grid_point = False 
        for sol_l in range( n_olgio_l ):
            # loop over the number of times each solvent in the solvent list needs to be replicated
            if( bad_grid_point ): break 
            for struc_i in sol_array:
                overlap_found = True
                if( bad_grid_point ): break 
                while ( overlap_found ):

                    dim_indx += 1
                    if( dim_indx > 2 ): dim_indx = -1
                    lat_indx[dim_indx]  += 1

                    l_x = float(lat_indx[0])*sol_length
                    l_y = float(lat_indx[1])*sol_length
                    l_z = float(lat_indx[2])*sol_length

                    if( l_x > oligomer_rep.lattvec[0][0] or l_y > oligomer_rep.lattvec[1][1] or l_y > oligomer_rep.lattvec[1][1] ):
                        bad_grid_point = True
                        break 
                    
                    lat_pos = numpy.array( [l_x,l_y,l_z] )

                    # Make sure there is no overlap with the added molecules
                    overlap = 0
                    for p_i, ptclObj_i in oligomer_rep.ptclC :
                        r_i = np.array( ptclObj_i.position )
                        for p_j, ptclObj_j in struc_i.ptclC :
                            r_j = np.array( ptclObj_j.position )
                            r_ij_sq = pbcs.sq_drij_c(r_i,r_j,oligomer_rep.getLatVec() )
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
                        struc_i.vec_shift(lat_pos)

                        sol_rep += struc_i
                        
                        if( options.verbose ):
                            if( rank == 0  ):
                                print "      -  Molecule ",sys_sol_n," has been added to the system after "

                                if( options.ptime ):
                                    t_f = datetime.datetime.now()
                                    dt_sec  = t_f.second - t_i.second
                                    dt_min  = t_f.minute - t_i.minute
                                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                    print "        - with placement time ",dt_min," min ",dt_sec," seconds "
                    else:
                        overlap_found = True
                        
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
            # Save lattice vectors as to no loose any expansions
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
            # Expand the box size by lc_expand
            #
            oligomer_rep.expandLatVec(options.lc_expand)
                    
            if( options.verbose ):                
                if( rank == 0  ):
                    print '          - Lattice vectors will be expanded by (option lc_expand)',options.lc_expand
                    
            p.barrier() # Barrier for MPI_COMM_WORLD
            
    if( rank == 0 ):
        log_out.close()
        
if __name__=="__main__":
    main()
   

