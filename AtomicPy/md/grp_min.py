#! /usr/bin/env python
# Minimize indvidual groups 

# Dr. Travis Kemper
# NREL
# Initial Date 10/28/2013
# travis.kemper@nrel.gov





def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog <file to convert> [options] "
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-r","--restart", dest="restart", default= False ,action="store_true" , help="Check for finished calculations and use previous results")

    parser.add_option("--host", dest="host",type="string", default="",help=" name of cluster ")
    parser.add_option("--pbcs",  dest="pbcs", default= True , help=" Use periodic boundry conditions ")
    parser.add_option("--npros", dest="npros", type="int",default="1", help="Number of processors to run calculations")


    # Group 
    parser.add_option("--h_term", dest="h_term", default=True, help=" Hydrogen terminate dangeling bonds of group")

    parser.add_option("--software", dest="software", type="string", default="gromacs", help=" Software package to use for Mininimization calculations  ")
    
    # FF settings 
    parser.add_option("--atom_types", dest="atom_types", type="string", default="", help="Read atom types that will replace default elements ")    
    parser.add_option("--nlist_bonds", dest="nlist_bonds", default=True, help="Build neighbor list from bonds")
    parser.add_option("--find_rings", dest="find_rings", default=False, help="Find rings ")
    #parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
    parser.add_option("--ff_file",dest="ff_file",  type="string", default="ff-new.itp", help=" itp file to read in for ff parameters ")
    parser.add_option("--zero_q", dest="zero_q", default=False, help=" Zero atomic charges ")

    # Gromacs
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file ")
    parser.add_option("--out_gro", dest="out_gro", type="string", default="out.gro", help="Output gromacs structure file ")

    parser.add_option("--gromacs_sufix", dest="gromacs_sufix", type="string", default="", help=" Sufix for gromacs calculations such as _mpi ")
    parser.add_option("--gromacs_dir", dest="gromacs_dir", type="string", default="", help=" gromacs dir ")
    parser.add_option("--g_center", dest="g_center",  default=False, help=" center minimization ")
    parser.add_option("--load_gromacs", dest="load_gromacs", type="string",  default="", help=" module comand to load gromacs ")

    #Lammps
    parser.add_option("--lammp_dir", dest="lammp_dir", type="string",  default="", help=" Lammps executable directory ")
    
    # Groups
    parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
    parser.add_option("--opt_groups", dest="opt_groups", default=True, help=" Optimize groups individually ")
    parser.add_option("--rm_VS", dest="rm_VS", default=False, help=" remove virtual sites ")
    

    # QM Settings
    parser.add_option("--qm_method", dest="qm_method", type="string",default="HF", help="Method of QM calculation ")
    parser.add_option("--qm_basis", dest="qm_basis", type="string",default="6-31G", help="Basis set of QM calculation ")
    parser.add_option("--qm_kywd", dest="qm_kywd", type="string",default="", help="Key words for QM calculation ")
    parser.add_option("--qm_mult", dest="qm_mult", type="int",default="0", help=" Shift in default spin multiplicity ( singlet,doublet) QM calculation, allows for triplets ")
    parser.add_option("--qm_npros", dest="qm_npros", type="int",default="1", help="Number of processors to run secondary calculations")
    parser.add_option("--qm_exclude", dest="qm_exclude", type="string",default="1", help="Exclude ")
    parser.add_option("--qm_charge", type="int",action="append",default="0", help="Input gaussain log file ")
    parser.add_option("--qm_triplet", dest="qm_triplet", default=False,  help="Output triplets ")
    
    (options, args) = parser.parse_args()

    return options, args

def opt_groups( group_index_i,group_list_i, ELN,ASYMB,R,ATYPE,GTYPE,CHARGES,CHARN,AMASS,RESID,RESN,BONDS,ANGLES,DIH ,LV,NBLIST,NBINDEX,options):
    import os, sys
    import  gaussian, top, gromacs ,lammps , file_io, prop, atom_types, gaussian
    import numpy as np
    import datetime

    # Set debug 
    debug = 0
    p_time = 0
    
    # Read in parameter file
    
    FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(options.ff_file)

    # Intialize optimize as positions as the input positions 
    R_opt = R
    
    R_i = []
    for atom_i in range( len(ELN) ):
	R_i.append( R[atom_i] )
    
 
    # Loop over all groups 
    for g_indx in range( len(group_index_i) - 1):
	
        # Get group index range
        N_o = group_index_i[g_indx]
        N_f = group_index_i[g_indx + 1 ] - 1
        NA_g = N_f - N_o + 1
        if( options.verbose ): print "        Printing gaussian input for group ",g_indx," with atoms ",NA_g
	
        F_id = "g_" + str(g_indx) 

        # Record id
        group_id = F_id
	
	
        # Place group into local array
	N_i = []
        ELECTRONS_i = 0 
        for a_indx in range( N_o,N_f+1):
            i = group_list_i[a_indx]
	    if( options.rm_VS  ):
		if( int(ELN[i]) > 0 ):
		    ELECTRONS_i += ELN[i]
		    N_i.append( i )
		    
		    # print i,ELN[i]
	    else:
		ELECTRONS_i += ELN[i]
		N_i.append( i )
        
	# Extract group topology information
	ELN_j,ASYMB_j,R_j,ATYPE_j,GTYPE_j,CHARGES_j,RESID_j,RESN_j,CHARN_j,AMASS_j,BONDS_j,ANGLES_j,DIH_j,IMPS_j,REF_j,GHOST_j = top.pass_i(N_i,ELN,ASYMB,R,ATYPE,GTYPE,CHARGES,CHARN,AMASS,RESID,RESN,BONDS,ANGLES,DIH,options)
	
	if( options.zero_q ):    
	    for atom_i in range( len(CHARGES_j)):
		CHARGES_j[atom_i] = 0.0 
	    

	#print_oniom(  id_name, ASYMB,R,ATYPE,CHARGES,ELECTRONS_i,Q,M,FIX,ONMTAG,options)
	if( options.software == "gromacs"):

	    # Create input files
	    g_id = "min_"+str(g_indx)
	    g_top = g_id + '.top'
	    g_gro = g_id + '.gro'
	    g_mdp = g_id + '.mdp'
	    g_log = g_id + '.log'
	    rm_bck = "rm \"#\"*"
	    s_suf = options.gromacs_sufix + '_d'

	    if( file_io.file_exists(g_top) ): os.remove(g_top)
	    if( file_io.file_exists(g_gro) ): os.remove(g_gro)
	    if( file_io.file_exists(g_mdp) ): os.remove(g_mdp)
	    if( file_io.file_exists(g_log) ): os.remove(g_log)
	    os.system(rm_bck)

	    gromacs.print_gro(g_gro,GTYPE_j,RESID_j,RESN_j,R_j,LV)

	    # Identify total number of atom types for lammps output 
	    ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ATYPE_j,AMASS_j,BONDS_j,ANGLES_j,DIH_j)

	    BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
	    ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
	    DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,options)	    
	    IMPTYPE_F  = top.imp_parameters()

	    const = []
	    angle = []

	    gromacs.print_top( g_top,ASYMB_j , ELN_j,ATYPE_j, GTYPE_j, CHARN_j , CHARGES_j, AMASS_j,RESN_j, RESID_j ,BONDS_j , ANGLES_j , DIH_j , IMPS_j
		    ,const,angle
		    ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
		    ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)

	    gromacs.print_min_nopbcs(g_mdp)

	    #sys.exit(" debug gro ")

	    # Run minimization 
	    ff_energy_f = gromacs.run_gromacs(g_gro,g_top,g_mdp,s_suf,options )

	    # Get optimize positions 
	    R_k = gromacs.get_coord(g_id,options)
	
	elif( options.software == "gaussian" ):
	    
	    # Create input files
	    g_id = "min_"+str(g_indx)
	    g_com = g_id+"/"+g_id + '.com'
	    g_chk = g_id+"/"+g_id + '.chk'
	    g_fchk = g_id+"/"+g_id + '.fchk'
	    g_log = g_id+"/"+g_id + '.log'
	    
	    run_calc = gaussian.check_fchk( g_fchk )
	    
	    # Run calculation if results are there or not 
	    if( not options.restart ): run_calc = 1
	    
	    if( run_calc  ):
		if( options.verbose ): print " running ",g_id,options.restart
    
		if( file_io.file_exists(g_com) ): os.remove(g_com)
		if( file_io.file_exists(g_chk) ): os.remove(g_chk)
		if( file_io.file_exists(g_fchk) ): os.remove(g_fchk)
		if( file_io.file_exists(g_log) ): os.remove(g_log)
		
	    
		ELECTRONS_j = 0
		for atom_j in range(len(ELN_j)):
		    ELECTRONS_j += ELN_j[atom_j] 
		
		gaussian.print_com( g_id, ASYMB_j,R_j,ATYPE_j,CHARGES_j,ELECTRONS_j,options.qm_method,options.qm_basis,options.qm_kywd,options.qm_charge,options.qm_mult)
		gaussian.run(options,g_id)
		
		#fchk_file = g_id +"/"+g_fchk
		
	    NA, ELN_k, R_k, TOTAL_ENERGY_k , Q_ESP = gaussian.parse_fchk( g_fchk )

			
	#
	#elif( options.ff_software == "lammps" ):
	#    
	#    rest_file = 'rest.in'
	#    data_file = "out.data"
	#    lammps.print_rest(rest_file,data_file,DIH_CONST_ANGLE,DIH_ATOMS[cent_indx],options)
	#    
	#    # Print lammps structure file 
	#    data_file = "out.data"
	#    lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
	#       BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
	#       ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
	#       DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
	#       RESN,ATYPE_IND,CHARGES,R , ATYPE,
	#       BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LAT_CONST)
	#    
	# Update positions of non ghost atoms 
	for atom_j in range(len(ELN_j)):
	    if( GHOST_j[atom_j] == 0 ):
		atom_i = REF_j[atom_j]
		R[atom_i] = R_k[atom_j]
			


    
	    
    # update VS positions
    if( options.rm_VS ):
	for atom_i in range( len(ELN) ):
	    if( int(ELN[atom_i]) <= 0 ):
		N_o = NBINDEX[atom_i ]
		N_f = NBINDEX[atom_i+ 1 ] - 1
		for indx_j in range( N_o,N_f+1):
		    atom_j = NBLIST[indx_j]
		    #print "  j - ",ATYPE[atom_j] ,R_i[atom_j] ,R[atom_j]
		    dr =   R[atom_j] - R_i[atom_j]
			
		    #print "    i- " ,atom_i,ATYPE[atom_i],R[atom_i] ,dr
		    R[atom_i] = R[atom_i] + dr
		    #print "      r_i new ",R[atom_i]
			
    return R
    
def main():
    import os, sys
    import gromacs, elements, top, groups

    options, args = get_options()

    # Print options
    #simpy.print_options(options)
    
    if( options.verbose ):
        print "  software ",options.software 
        print "  pbcs ",options.pbcs
        print "  "

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
    
    if( options.opt_groups ):
        if( options.verbose ): print "  Print groups into structure files"
        R_opt = opt_groups(group_index_i,group_list_i, ELN,ASYMB,R,ATYPE,GTYPE,CHARGES,CHARN,AMASS,RESID,RESN,BONDS,ANGLES,DIH ,LV,NBLIST,NBINDEX,options)

    # Update virtual site postions if 
    if( options.group_ptma and options.software == "gaussian" ):
        if( options.verbose ): print "  Update virtual site postions "
        group_index_i,group_list_i,group_numb,group_cnt = groups.tempo(  ATYPE,ELN,NBLIST,NBINDEX, options )
            
    ##os.chdir(work_dir)
    #options.out_gro = "min_groups.gro"
    if( options.out_gro ):
        if( options.verbose ): print "  Writing gro file ",options.out_gro
	gromacs.print_gro(options.out_gro,GTYPE,RESID,RESN,R_opt,LV)
        
        #pbc_whole = "echo -e \" 0 \n \" | trjconv -f " + options.out_gro -s ${RNID}.tpr -trans 1 1 1 -o  ${ID}-MIN.gro -pbc whole
        

if __name__=="__main__":
    main()