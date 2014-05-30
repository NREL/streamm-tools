#! /usr/bin/env python
"""
Create a supercell of randomly placed and rotated molecules
 length - angstroms
 mass   - AMU
 volume - angstroms^3
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
    
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--itp", dest="itp_file",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")

    parser.add_option("--atomic_cut", dest="atomic_cut", type=float, default=2.5, help="Minimum distance between atoms of molecules ")

    parser.add_option("--den_target", dest="den_target", type=float, default=0.10, help="Target density g/cm^3 ")
    parser.add_option("--atoms_target", dest="atoms_target", type=float, default=100000.0, help="Target number of atoms ")
    parser.add_option("--max_mol_place", dest="max_mol_place", type=float, default=1000, help="Maximum attempts to place a molecule  ")
    parser.add_option("--max_sys", dest="max_sys", type=float, default=10, help="Maximum system recreations at a certain lattice constant ")
    parser.add_option("--lc_expand", dest="lc_expand", type=float, default=2.5, help="Distance (angstroms) to increase system size after max_sys is excieded ")

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

#
# Partition list on 'parts' numbers with last chunk
# is the remainder of elements
#  * data  --> eg [1.2, 3.3........38.2]
#  * parts --> 2
#  * chunks = [ [1.2, 3.3], [.....] ]
#
def partitionList(data, partLength):
    import boost.mpi as mpi
    import math

    # Make 'parts' number of chunks
    chunks=[data[x:x+partLength] for x in xrange(0, len(data), partLength)]
    return chunks

#
# Given a list and a number of processors (given by mpi.size)
# split up list equally (as possible).
# Driver for partitionList(....)
# NOTE: no communication is performed... all data known globally
# 
def splitOnProcs(data):
    import boost.mpi as mpi
    import sys, math


    parts     = mpi.size                       # Number of procesors
    numData   = float(len(data))               # Global size of data list
    dataPProc = int(math.ceil(numData/parts))  # Est. num of data on each proc
    plist     = partitionList(data, dataPProc) # Split data into chunks
    mpi.world.barrier()

    # Error check or return results
    if len(plist) != parts:
        if mpi.rank == 0:
            print " " 
            print "Partitioning failed"
            print " ... check data length and number of processors requested"
            print " len(plist) = ", len(plist)
            print "      parts = ", parts
        sys.exit(0)
    else:
        return plist[mpi.rank]
   
def replicate_dih(N_repeats,ELN_i,DIH_i,DTYPE_IND_i):
    """
    Replicate dihedrals for given topology  to add mocules or groups to a system
    """
    
    debug = 0
    
    DIH_sys = []
    DTYPE_IND_sys = []
    
    NA_i = len( ELN_i)
        
    for unit_n in range( N_repeats ):
        
        # Repeat bonds
        for bond_indx in range(len(DIH_i)):
            k_o = DIH_i[bond_indx][0]
            i_o = DIH_i[bond_indx][1]
            j_o = DIH_i[bond_indx][2]
            l_o = DIH_i[bond_indx][3]
            k_add = k_o + ( unit_n  )*NA_i
            i_add = i_o + ( unit_n  )*NA_i
            j_add = j_o + ( unit_n  )*NA_i
            l_add = l_o + ( unit_n  )*NA_i
            DIH_sys.append( [k_add,i_add ,j_add,l_add] )
            DTYPE_IND_sys.append(DTYPE_IND_i[bond_indx] )
            if( debug ):
                print i_o,j_o, " -> ",i_add ,j_add
            
    return DIH_sys,DTYPE_IND_sys

def replicate_angles(N_repeats,ELN_i,ANGLES_i,ANGTYPE_IND_i):
    """
    Replicate angles for given topology  to add mocules or groups to a system
    """
    
    debug = 0
    
    ANGLES_sys = []
    ANGTYPE_IND_sys = []
    NA_i = len( ELN_i)
        
    for unit_n in range( N_repeats ):
        
        # Repeat bonds
        for bond_indx in range(len(ANGLES_i)):
            k_o = ANGLES_i[bond_indx][0]
            i_o = ANGLES_i[bond_indx][1]
            j_o = ANGLES_i[bond_indx][2]
            k_add = k_o + ( unit_n  )*NA_i
            i_add = i_o + ( unit_n  )*NA_i
            j_add = j_o + ( unit_n  )*NA_i
            ANGLES_sys.append( [k_add,i_add ,j_add] )
            ANGTYPE_IND_sys.append(ANGTYPE_IND_i[bond_indx] )
            
            if( debug ):
                print i_o,j_o, " -> ",i_add ,j_add
            
    return ANGLES_sys,ANGTYPE_IND_sys


def replicate_bonds(N_repeats,ELN_i,BONDS_i,BTYPE_IND_i):
    """
    Replicate bonds for given topology  to add mocules or groups to a system
    """
    
    debug = 0
    
    BONDS_sys = []
    BTYPE_IND_sys = []
    NA_i = len( ELN_i)
        
    for unit_n in range( N_repeats ):
        
        # Repeat bonds
        for bond_indx in range(len(BONDS_i)):
            i_o = BONDS_i[bond_indx][0]
            j_o = BONDS_i[bond_indx][1]
            i_add = i_o + ( unit_n  )*NA_i
            j_add = j_o + ( unit_n  )*NA_i
            BONDS_sys.append( [i_add ,j_add] )
            BTYPE_IND_sys.append(BTYPE_IND_i[bond_indx] )
            if( debug ):
                print i_o,j_o, " -> ",i_add ,j_add
            
    return BONDS_sys,BTYPE_IND_sys


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


    use_mpi = 1
    if( use_mpi ):
            
        import boost.mpi as mpi
        
        #
        # MPI startup
        #
        world = mpi.world
        mpi.world.barrier() # Barrier for MPI_COMM_WORLD
        
    #
    # Load information onto all processors 
    #
    options, args = get_options()
    prop_dim = 3
    ang_acc = 1000  # number of digets in random angle 
    
    
    debug = 0
    p_time = 1

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
    #
    # Read in parameter file 
    #


    #
    # Read in ff file
    #
    if( len(options.itp_file) ):
        if( options.verbose ): print  "     - Reading in ",options.itp_file    
        FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(options.itp_file)
        if(  options.ff_software == "lammps"  ):
            # Identify total number of atom types for lammps output 
            ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND_i , BTYPE_REF, ANGTYPE_IND_i , ANGTYPE_REF, DTYPE_IND_i , DTYPE_REF = lammps.lmp_types(ELN_i,ATYPE_i,AMASS_i,BONDS_i,ANGLES_i,DIH_i)

            #   Build covalent nieghbor list for bonded information 
            NBLIST, NBINDEX = top.build_covnablist(ELN_i,R_i)
            
            # Check atom types to be sure each atom of the same type has the same number of neighbors 
            ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE_i,NBLIST,NBINDEX)
        
            ATYPE_EP, ATYPE_SIG = top.atom_parameters(options.itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
            BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(options.itp_file,BTYPE_IND_i , BTYPE_REF,FF_BONDTYPES)
            ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(options.itp_file,ANGTYPE_IND_i , ANGTYPE_REF,FF_ANGLETYPES)
            DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(options.itp_file, options.norm_dihparam, DTYPE_IND_i , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )
            
            IMPTYPE_F  = top.imp_parameters(options.itp_file)
    elif(len(options.in_data) == 0 ):
        if(  options.ff_software == "lammps"  ):
            print " An itp file specified with the --itp_file option is needed to create a lammps input file "
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
    # Test that geometry was read in
    #
    try:
        R_i
    except NameError:
        sys.exit("Geometry read in error ")
    
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
    # Calculate need molecules to achieve specified density
    #
    mol_mult = int( options.atoms_target/float(len(ELN_i)) )
    target_density_amuang = options.den_target*const_avo/10.0
    mol_mass_amu = prop.total_mass( AMASS_i )
    mass_amu = mol_mass_amu*mol_mult
    volume_target_ang = mass_amu/target_density_amuang    
    len_target_ang = volume_target_ang**(1.0/3.0)
    
    """
    Set lattice vectors for pbc 
    """
    LV = numpy.zeros([3,3])
    
    LV[0,0] = len_target_ang
    LV[1,1] = len_target_ang
    LV[2,2] = len_target_ang

    #
    # Find target volume for single molecule
    #
    mol_vol = float(mol_mult)/volume_target_ang
    mol_unit_l =  mol_vol**(-1.0/3.0)
    
    
    
    #
    # Initialize system
    #
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
    GTYPE_sys = []
    R_sys = []
    VEL_sys = []
    BONDS_sys = []
    BTYPE_IND_sys = []
    ANGLES_sys = []
    ANGTYPE_IND_sys = []
    DIH_sys = []
    DTYPE_IND_sys = []
        
    sys_mol_n = 0    # Number of molecules add to the system
    sys_attempts = 0  # Number of times the system has been reset
    
    #
    # Set cut-off squared 
    #

    cut_ij_sq = options.atomic_cut* options.atomic_cut

    
    if( use_mpi ):
        if( mpi.rank == 0 and options.verbose ):
                
            print "   - Tragets "
            print "     Input molecule has ",len(ELN_i)," atoms and mass of ",mol_mass_amu," AMU "
            print "     To achieve ",options.atoms_target," atom systems it will be multiplied ",mol_mult
            print "     giving ",mol_vol," mol/Angstrom^3 and a molecular unit length of ",mol_unit_l," Angstorms "
            print "     For the target density of ",options.den_target," g/cnm^3 ",target_density_amuang," AMU Angstrom^-3"
            print "       a target cubic unit cell of ",len_target_ang," will be needed "
            print "   - Options "
            print "     in_top",options.in_top
            print "     in_gro",options.in_gro
            print "     atomic_cut",options.atomic_cut
            print "     den_target",options.den_target
            print "     atoms_target",options.atoms_target
            print "     max_mol_place",options.max_mol_place
            print "     max_sys",options.max_sys
            print "     lc_expand",options.lc_expand
            print "     out_gro",options.out_gro

    elif( options.verbose ):
        print "   - Tragets "
        print "     Input molecule has ",len(ELN_i)," atoms and mass of ",mol_mass_amu," AMU "
        print "     To achieve ",options.atoms_target," atom systems it will be multiplied ",mol_mult
        print "     giving ",mol_vol," mol/Angstrom^3 and a molecular unit length of ",mol_unit_l," Angstorms "
        print "     For the target density of ",options.den_target," g/cnm^3 ",target_density_amuang," AMU Angstrom^-3"
        print "       a target cubic unit cell of ",len_target_ang," will be needed "
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
        
    
    if( use_mpi ):
        

        mpi.world.barrier() # Barrier for MPI_COMM_WORLD
            
        pointIndices = range( len(ELN_i)  )
        if( debug ): print mpi.rank, mpi.size," splitOnProcs "
        myChunk  = splitOnProcs(pointIndices)
            
        if(debug):                
            # mpi.world.recv(R_shift)
            print " cpu ",mpi.rank ," has atoms ",myChunk[0]," - ",myChunk[len(myChunk)-1],"  \n"
                    
                    
                            
        mpi.world.barrier() # Barrier for MPI_COMM_WORLD

    add_mol = 1
        
    while ( add_mol ):
    
        add_mol = 1
        overlap_sum = 1
        moladd_atempts = 0
        

        while ( overlap_sum ):
            moladd_atempts += 1
            
            if( use_mpi ):
                # Declare on all processors
                R_shift_o = []
                
                if ( mpi.rank == 0 ):
                    #
                    # Get coordinates of randomly rotated and shifted molecule
                    #
                    R_shift_o = ran_rot_shift(R_moli_c,len_target_ang,ang_acc )
                    
            else:
                #
                # Get coordinates of randomly rotated and shifted molecule
                #
                R_shift = ran_rot_shift(R_moli_c,len_target_ang,ang_acc )
                

            if( use_mpi ):
    
                if( debug ):
                    print mpi.rank, mpi.size,"  initialization finished "
                    if(  mpi.rank == 0 ):
                        print " R_shift first coordinate ",R_shift_o[0]," on ",mpi.rank
                    print
                #
                # Broadcast molecular position from processor 0 to all other processors 
                #
                R_shift = mpi.broadcast(world,R_shift_o,0)
                
            
            #
            # Check molecules do not overlap
            #
            overlap = 0
            if( len(ELN_sys) > 0 ):
                if( use_mpi ):
                
                    
                    for atom_i in myChunk:
                        r_i = R_shift[atom_i]
                        for sys_atom in range(len(ELN_sys)):
                            r_j = R_sys[sys_atom]
                            r_ij_sq = prop.sq_drij_c(r_i,r_j,LV)
                            if( r_ij_sq < cut_ij_sq ):
                                overlap = 1
                
                    mpi.world.barrier() # Barrier for MPI_COMM_WORLD
                    

                else:
                        
                    for atom_i in range(len(ELN_i)):
                        r_i = R_shift[atom_i]
                        for sys_atom in range(len(ELN_sys)):
                            r_j = R_sys[sys_atom]
                            r_ij_sq = prop.sq_drij_c(r_i,r_j,LV)
                            if( r_ij_sq < cut_ij_sq ):
                                overlap = 1
                            
            
            if( use_mpi ):    
                overlap_sum = mpi.all_reduce(world,overlap, lambda x,y: x + y)
            else:
                overlap_sum = overlap
                

            if( overlap_sum ==  0 ):
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
                    ATYPE_IND_sys.append( ATYPE_IND_i[atom_i])
                    
                
                if( options.verbose ):
                    if( use_mpi ):
                        if( mpi.rank == 0 ):
                            print "      -  Molecule ",sys_mol_n," has been added to the system after ",moladd_atempts," placment attempts "
                    else:
                        print "      -  Molecule ",sys_mol_n," has been added to the system after ",moladd_atempts," placment attempts "

            if( moladd_atempts >= options.max_mol_place ):
                
                if( options.verbose ):
                    if( use_mpi ):
                        if( mpi.rank == 0 ):
                            print "        -  Attempts to add molecule ",sys_mol_n," has exceeded max attempts ",options.max_mol_place," system will be reset for the ",sys_attempts," time "
                    else:
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
                    
                    if( options.verbose ):
                            
                        if( use_mpi ):
                            if( mpi.rank == 0 ):
                                        
                                print '          - Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                                print '          - Lattice vectors will be expanded by (option lc_expand)',options.lc_expand
                                
                        else:
                            
                            print '          - Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                            print '          - Lattice vectors will be expanded by (option lc_expand)',options.lc_expand
                            
                    mol_system[0,0] = LV[0,0] + options.lc_expand
                    LV[1,1] = LV[1,1] + options.lc_expand
                    LV[2,2] = LV[2,2] + options.lc_expand
                
            if( use_mpi): mpi.world.barrier() # Barrier for MPI_COMM_WORLD
                    
                                    
        if( sys_mol_n ==  mol_mult ):
            add_mol = 0
            if( use_mpi ):
                mpi.world.barrier() # Barrier for MPI_COMM_WORLD
                if( options.verbose and mpi.rank == 0 ):
                   print " All molecules have been added "
            else:
                if( options.verbose ):
                    print " All molecules have been added "
                
                
    if( use_mpi ):
            
        if( mpi.rank == 0 ):
            out_xyz = "mol_system.xyz"
            xmol.write_xyz(ASYMB_sys,R_sys,out_xyz)
            
            if( options.ff_software == "gromacs" ):
                gromacs.print_gro(options.out_gro,GTYPE_sys,RESID_sys,RESN_sys,R_sys,LV)

            elif(  options.ff_software == "lammps"  ):
                # Find topology for entire system
                #   Repeat molecular topology for all molecules added to the system
        
                BONDS_sys,BTYPE_IND_sys = replicate_bonds(mol_mult,ELN_i,BONDS_i,BTYPE_IND_i) 
                ANGLES_sys,ANGTYPE_IND_sys = replicate_angles(mol_mult,ELN_i,ANGLES_i,ANGTYPE_IND_i) 
                DIH_sys,DTYPE_IND_sys = replicate_dih(mol_mult,ELN_i,DIH_i,DTYPE_IND_i) 

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
            

    else:
        
        out_xyz = "mol_system.xyz"
        xmol.write_xyz(ASYMB_sys,R_sys,out_xyz)
        gromacs.print_gro(options.out_gro,GTYPE_sys,RESID_sys,RESN_sys,R_sys,LV)
            
    
if __name__=="__main__":
    main()
   
