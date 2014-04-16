#! /usr/bin/env python
# Make input files for torsional potential energy surface 

# Dr. Travis Kemper
# NREL
# 12/09/2013
# travis.kemper@nrel.gov

def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    usage = usage + "Input files \n"
    usage = usage + "  specify the destination name of an option followed by the value"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")
    
    
    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")

    # How to run the needed calculations 
    parser.add_option("--submit", dest="submit",action="store_true", default=False,help=" submit calculations to the queue ")
    parser.add_option("--localrun", dest="localrun",action="store_true", default=False,help=" Run calculations locally")
    parser.add_option("--submit_command", dest="submit_command",type="string", default="qsub",help=" command used to submit script to the queue ")

    parser.add_option("--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")

    # Torsion
    parser.add_option("--cent_min", dest="cent_min", type="int", default="0",help=" Initial torsional angle ")
    parser.add_option("--cent_max", dest="cent_max", type="int", default="180",help=" Final torsional angle ")
    parser.add_option("--cent_step", dest="cent_step", type="int", default="5",help=" Step size torsional angle ")
    

    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    
    # Force field generation options     
    parser.set_defaults(ff_software="gromacs")
    parser.add_option("--ff_software", dest="ff_software",type="string",help=" what software to use for the ff calculations   ")

    parser.add_option("--itp", dest="itp_file",  type="string", default="oplsaa_biaryl.itp",help="gromacs force field parameter file")
    parser.set_defaults(set_pmma=False)

    parser.set_defaults(ff_charges=False)
    parser.add_option("--ff_charges", dest="ff_charges",action="store_true",help=" Use ff charges ")
    
    parser.set_defaults(zero_charges=False)
    parser.add_option("--zero_charges", dest="zero_charges",action="store_true",help=" Zero charges ")
    
    parser.add_option("--norm_dihparam", dest="norm_dihparam", help="Normalize dihedral potential terms if single dihedral is specified in itp file  ", default=0)

    # Gromacs related options 
    parser.add_option("--gromacs_dir", dest="gromacs_dir",type="string",default="",help=" Directory of gromacs run files   ")
    parser.add_option("--gromacs_sufix", dest="gromacs_sufix",type="string",default="",help=" sufix for gromacs such as _d or _mpi  ")
    
    # Lammps options 
    parser.set_defaults(lammp_dir="$HOME/Software/lammps/src/")
    parser.add_option("--lammp_dir", dest="lammp_dir",type="string",help=" Directory of lammps run files  ")

    parser.add_option("--pmem",dest="pmem", type="int", default="1700",help=" Memory per processor ")
    parser.add_option("--npros", dest="npros", type="int", default="4",help=" Number of processors ")
    parser.add_option("--nnodes", dest="nnodes", type="int", default="1",help=" Number of nodes ")

    parser.add_option("--high_basis", dest="high_basis", type="string", default="cc-pVTZ",help=" Basis set for hihgh level energy calculations ")
    parser.add_option("--dih_temp", dest="dih_temp", type="string", default="mp2_dih.com.template",help=" Template for Links of dihedral calculation ")

    parser.add_option("--limdih", dest="limdih", default=False,action="store_true", help="Limit the number of dihedrals per ij pair ")
    parser.add_option("--limitdih_n", dest="limitdih_n", type="int", default="1",help=" Number of dihedrals per ij pair with limdih ")

    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")

    (options, args) = parser.parse_args()
    
    return options, args

def submit_job( struct_dir, pbs_id ,options):
    import sys, os

    print " submitting job " , pbs_id
    #submit = "sbatch " + pbs_id
    submit = options.submit_command +" "+ pbs_id

    #print " sumitting ",submit
    #sys.exit(' submit_job ')
    
    os.system(submit)
 
    return 
    #sys.exit('test job')

def build_nablist(ELN,BONDS):
    import sys,numpy

    debug = 0
    NNAB = 0

    maxnnab = len(BONDS)*2 + 1
    NBLIST = numpy.empty( maxnnab,  dtype=int )
    NBINDEX = numpy.empty( maxnnab,  dtype=int )

    if(debug ):
	for b in range( len(BONDS) ) :
	    bnd_i = BONDS[b][0]
	    bnd_j = BONDS[b][1]
	    print b, bnd_i,bnd_j
	    
    for i in range(len(ELN)):
	NBINDEX[i] = NNAB + 1
	for b in range( len(BONDS) ) :
	    bnd_i = BONDS[b][0]
	    bnd_j = BONDS[b][1]
	    if ( i == bnd_i ):
		NNAB = NNAB + 1
		NBLIST[NNAB] =  bnd_j
		if(debug): print " adding bnd_i",i,b,bnd_i,bnd_j,len(NBLIST),NBINDEX[i]
	
	    if ( i == bnd_j ):
		NNAB = NNAB + 1
		NBLIST[NNAB] =  bnd_i
		if(debug): print " adding bnd_j",i,b,bnd_i,bnd_j,len(NBLIST),NBINDEX[i]
		    
    #if(debug): sys.exit('debug')
    
    # Account for final atom position
    NBINDEX[i+1] =  NNAB + 1

    debug = 0
    if ( debug ):
       print ' total nbs ',NNAB,' total bonds ',len(BONDS)
       for i in range(len(ELN)):
	    N_o = NBINDEX[ i  ]
	    N_f = NBINDEX[ i + 1 ] - 1
	    NNAB = N_f - N_o + 1
	    print ' atom ', i,' ',ELN[i],NNAB,N_o
    
	    #
	    # Find number of elements
	    #
	    #ELCNT = numpy.zeros(120)
	    for indx in range( N_o,N_f+1):
		j = NBLIST[indx]
		print ELN[j],j
		#    el_j = ELN[j]
		#    ELCNT[j] = ELCNT[j] + 1

       sys.exit('mk_ff.build_nablist')


    return (NBLIST,NBINDEX)


def read_dihlist(dlist_name):
    
    
    DIH_ID = []
    DIH_VAL = []
    DIH_ATOMS = []
    DIH_TAG = []
    
    # get lines 
    f = open(dlist_name,'r')
    Lines = f.readlines()
    f.close()
    
    
    for line in Lines: 
        col = line.split()
        if( col[0] == "dih" and len(col) >= 11 ):
            dih_indx = col[1]
            DIH_TAG.append( col[2] )
            DIH_ID.append( col[3] )
            DIH_VAL.append(  col[4] )
            a_l = int( col[5] ) - 1
            a_i = int( col[6] ) - 1
            a_j = int( col[7] ) - 1
            a_k = int( col[8] ) - 1
            DIH_ATOMS.append( [ a_l,a_i,a_j,a_k ] )

    return (  DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS )
    
def check_dihedrals( dih_a, dih_b ):
    
    ij_same = 0
    kl_same = 0
    dih_same = 0 
    if( dih_a[1] == dih_b[1]   and  dih_a[2] == dih_b[2] ):
	ij_same = 1 
    if( dih_a[1] == dih_b[2]   and  dih_a[2] == dih_b[1] ):
	ij_same = 1 
    if( ij_same ):
	    
	if( dih_a[0] == dih_b[0]   and  dih_a[3] == dih_b[3] ):
	    kl_same = 1 
	if( dih_a[0] == dih_b[3]   and  dih_a[3] == dih_b[0] ):
	    kl_same = 1
	if( kl_same ):
	    dih_same = 1
	    
    return dih_same

def main():
    import sys, os 
    import string  , numpy , json 
    from string import replace
    import jsonapy
    import file_io, gaussian, elements, gromacs, lammps , top, atom_types, xmol 
    
    options, args = get_options()
    
    #
    if( options.host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	
	if( options.ff_software == "gromacs"):
	    md_min_template = "peregrine.md_min.template"
	    f = open(md_min_template,"r")
	    md_min_templ = f.read()
	    f.close()
	elif( options.ff_software == "lammps" ):	
	    lammps_min_template = "peregrine.lmp_min.template"
	    f = open(lammps_min_template,'r')
	    md_min_templ = f.read()
	    f.close()
	
        f = open("peregrine.pbs.template",'r') # redmesa.slurm.template	
	pbs_templ = f.read()
	f.close()

    elif( options.host == "macbook" ):
	load_gaussian = ""
        user = 'tkemper'

	
    # sufix for force field calcs to run multiple ff types or excultions, such as q(i) = 0.00 
    ff_type_id = "_fit"
    
    # Read in ff file
    FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(options.itp_file)
    
    LV = numpy.zeros( (3,3) )
    
    LV[0][0] = 50.0
    LV[1][1] = 50.0
    LV[2][2] = 50.0

    # Store working dir  
    work_dir = os.getcwd()
    
    json_files = options.json.split(',')
    print json_files
    if( len(json_files) > 0 ):
	# Read index files from args
	for json_file in json_files:
	    
	    # Verbose output
	    if( options.verbose ):
		print "The molecules specified in json file ",options.json," will be read in "
    
	    json_data,json_success = jsonapy.read_jsondata(json_file)
	    if(  json_success ):
		
		mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
		dih_id_list ,cent_min_list ,cent_max_list ,cent_step_list,a_k_list, a_i_list, a_j_list, a_l_list,qmtor_found = jsonapy.read_qm_tor(json_data)
		
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found and qmtor_found  ):
		    json_atomicdata = 0
		    fchk_atomicdata = 0
		    
		    if( options.verbose ):
			print " Meta data found will use specified method and basis unless others are specified in the options "
		    #
		    # Construct file names 
		    #
		    #short_name = "acc%d_%s_n%d" % (accuracy, tag, number )
		    job_name = "acc%d_%s_n%d" % (accuracy, tag, n_units )
		    struct_dir = "%s/%s/" % (mol_dir, tag )
		    
		    calc_id = "%s/%s%s" % (struct_dir, job_name , "-ZMATOPT" )
		    
		    
				    
		    #append_qm_tor_json
		    ff_tor_data = {}
		    json_data['metadata']["ff_tor_data"] = ff_tor_data
		    #
		    #
		    ff_tor_data["cent_id"] = []
		    ff_tor_data["a_k"] =  []
		    ff_tor_data["a_i"] =  []
		    ff_tor_data["a_j"] =  []
		    ff_tor_data["a_l"] =  []
		    ff_tor_data["cent_min"] =  []
		    ff_tor_data["cent_max"] =  []
		    ff_tor_data["cent_step"] = []
		    ff_tor_data["ff_type"] = []
		    #
		    # 
		    #
		    zmat_fchk = "%s/%s%s" % ( calc_id ,job_name,"-ZMATOPT.fchk" )
		    print " Checking for complete zmatrix optimiztion ",zmat_fchk
                    zmat_finished = file_io.file_exists( zmat_fchk )
		    
		    
		    
		    if( zmat_finished  ):
			
			    
			# Read in from zmatrix optimization 
			NA, ELN, R, TOTAL_ENERGY, CHARGES   = gaussian.parse_fchk( zmat_fchk )
			ASYMB = elements.eln_asymb(ELN)
						
			logz_name = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.log" )
			comz_name = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.com" )
			
			# Get bonds from log file
			if( options.verbose ):
			    print "  Getting bonds from optimization log file ",logz_name
			BONDS, ANGLES, DIH  = gaussian.read_optlog(logz_name)
			NBLIST,NBINDEX = build_nablist(ELN,BONDS)
	
			# pars zmatirix file
			f = open(comz_name,'r')
			zmatrix_lines = f.readlines()
			f.close()
			
			# need to change to reading out of connections 
			#DIH_ID,DIH_ATOMS = get_dih_id( RING_CONNECT,RING_NUMB, zmatrix_lines )
			# Find rings		    
			RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
			
			RING_CONNECT  = top.find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB) 
    
			#DIH_ID,DIH_ATOMS,zmatrix = gaussian.get_dih_id( RING_CONNECT,RING_NUMB, zmatrix_lines )
			
			zmatrix = gaussian.com_zmatrix(comz_name)    
			DIH_ID, DIH_VAL, DIH_ATOMS = gaussian.get_dih_id( zmatrix)
			
			
			os.chdir(struct_dir)
    
			dlist_name = job_name + "_dih.list"
    
			if( options.verbose ):
			    print "  Getting dihedral information from ",dlist_name
			    
			dlist_exists = file_io.file_exists( dlist_name )
			if ( dlist_exists ):                    
			    if( options.verbose ):
				print "       Writing input files for all loop dihedrals "
				print "         Nodes ",options.nnodes
				print "         Processors  ",options.npros
			    DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS = read_dihlist(dlist_name)
				
			    ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
			    
			    
			    
			    DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)
			    IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)
			    
			    GTYPE = top.initialize_gtype( ELN )
			    RESID = top.initialize_resid( ELN )
			    RESN = top.initialize_resn( ELN )
			    CHARN = top.initialize_charn( ELN )
			    
			    ASYMB = elements.eln_asymb(ELN)
			    AMASS = elements.eln_amass(ELN)
			
			    RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
			    RING_CONNECT  = top.find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB) 
	    
	
			    # Asign oplsaa atom types 
			    # ATYPE , CHARGES = atom_types.oplsaa( options, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
			    
	
			    ATYPE, CHARGES = atom_types.oplsaa(  options.ff_charges,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB )
			    
			    ATYPE , CHARGES = atom_types.biaryl_types( options.ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
			    #Refind inter ring types
			    ATYPE , CHARGES  = atom_types.interring_types(options.ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
			   
			    if( options.zero_charges ):
				for indx in range( len(CHARGES)):
				    CHARGES[indx] = 0.0 
				
			    # Set charge groups
					    
			    #
			    # Set charge groups
			    #
			    CG_SET = []
			    one = 1
			    for i in range( len(ELN) ):
				CG_SET.append(one)
				
			    #CHARN = top.set_chargegroups(CHARN,ATYPE,ASYMB,ELN,NBLIST,NBINDEX, RING_NUMB)
			    CHARN = top.set_chargegroups(options,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LV)
	
			    
			    debug = 0
			    if( debug):
				for i in range( len(ELN) ):
				    print i,ELN[i],ASYMB[i],ATYPE[i],CHARGES[i] # , CHARGES[i]
				    
				sys.exit( 'mk_ff')
				
			
	
			    # Identify total number of atom types for lammps output 
			    ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH)
	
			    # Check atom types to be sure each atom of the same type has the same number of neighbors 
			    ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)
			
			    # Modify dih potentials
			    ## FF_DIHTYPES = mod_dih( ATYPE_IND
			
			    # Print new itp file with only used atom types and interactions
			    if( options.ff_software == "gromacs"):
				##AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(options,ASYMB,ATYPE,BONDS,ANGLES,DIH,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
				new_itp = "ff-new.itp"
				AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(new_itp,options.norm_dihparam,ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
	
			    # Print output
			    
			    out_xyz =  job_name + "-ff.xyz"
			    xmol.write_xyz(ASYMB,R,out_xyz)
			    xmol.print_cgxyz(CHARN,R)
			    job_gro = job_name +  ".gro"
			    gromacs.print_gro(job_gro,GTYPE,RESID,RESN,R,LV)
			    
			    ATYPE_EP, ATYPE_SIG = top.atom_parameters(options.itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
			    BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(options.itp_file,BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
			    ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(options.itp_file,ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
			    DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(options.itp_file, options.norm_dihparam, DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )
			    
			    IMPTYPE_F  = top.imp_parameters(options.itp_file)
			    
			    if( options.ff_software == "gromacs"):
				top_file =  job_name +  ".top"
				const = []
				angle = []
				gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
					       ,const,angle
					       ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
					       ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
				
			    elif( options.ff_software == "lammps" ):
				data_file = "out.data" 
				lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
				      BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
				      ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
				      DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
				      RESN,ATYPE_IND,CHARGES,R , ATYPE,
				      BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV)
				#sys.exit("lammps test")
		    
			    # Add fixed dihedrals to fix list
			    DIH_CONST_ANGLE = [] #float( cent_angle )
			    DIH_CONST_ATOMS =[] # DIH_ATOMS			
    
			    for indx_c in range(len(DIH_ID)):
				if( DIH_TAG[indx_c].strip() == "fix" ):
			    
				    dih_i = DIH_ATOMS[indx_c][0] + 1
				    dih_j = DIH_ATOMS[indx_c][1] + 1
				    dih_k = DIH_ATOMS[indx_c][2] + 1
				    dih_l = DIH_ATOMS[indx_c][3] + 1
				    
				    DIH_CONST_ATOMS.append( [dih_i,dih_j,dih_k,dih_l] )
				    DIH_CONST_ANGLE.append( DIH_VAL[indx_c] )
				    
			    # Add loop as last fixed dihedral 
			    for indx_c in range(len(DIH_ID)):
				if( DIH_TAG[indx_c].strip() == "loop" ):
			    
				    dih_i = DIH_ATOMS[indx_c][0] + 1
				    dih_j = DIH_ATOMS[indx_c][1] + 1
				    dih_k = DIH_ATOMS[indx_c][2] + 1
				    dih_l = DIH_ATOMS[indx_c][3] + 1
				    
				    DIH_CONST_ATOMS.append( [dih_i,dih_j,dih_k,dih_l] )
				    DIH_CONST_ANGLE.append( DIH_VAL[indx_c] )
				    
			    # find compents of dihedral not in connection
			    debug = 0
			    DIH_fit = []
			    DIHTYPE_F_fit = []
			    DTYPE_IND_fit  = []
			    
			    # Fitted
			    dih_comp_name = job_name + "_fit.list"
			    f = open(dih_comp_name,'w')
			    DIH_fitset = []
			    
			    for indx_dih in range( len(DIH) ):
				dih_i = DIH[indx_dih][0] + 1
				dih_j = DIH[indx_dih][1] + 1
				dih_k = DIH[indx_dih][2] + 1
				dih_l = DIH[indx_dih][3] + 1
				#dih_real = [dih_i,dih_j,dih_k,dih_l]
				add_dih = 1
				if( debug ): print " checking ",indx_dih,DIH[indx_dih]
				for indx_c in range(len(DIH_ID)):
				    if( DIH_TAG[indx_c].strip() == "loop" ):
					if( debug ): print "   against loop atoms ",indx_c,DIH_ATOMS[indx_c]
					#same_dih = check_dihedrals(DIH[indx_dih],DIH_ATOMS[indx_c] )
					#if( same_dih ):
					#	add_dih = 0
					
					fix_j = DIH_ATOMS[indx_c][1] + 1 
					fix_k = DIH_ATOMS[indx_c][2] + 1 
					
					print " checking DIH ",dih_j,dih_k
					print " against looperf ",fix_j,fix_k
					
					fix_comp = 0
					if( dih_j == fix_j and dih_k == fix_k ):
					    fix_comp = 1
					if( dih_j == fix_k and dih_k == fix_j ):
					    fix_comp = 1
					    
					if( fix_comp ):
					    DIH_fitset.append( [ dih_i-1,dih_j-1,dih_k-1,dih_l-1 ])
					    f.write( " dihedral %8d  %8d  %8d  %8d %s %s %s %s  \n" % ( dih_i,dih_j,dih_k,dih_l,ATYPE[dih_i-1],ATYPE[dih_j-1],ATYPE[dih_k-1],ATYPE[dih_l-1]))
					    add_dih = 0
    
					    print "      dih ",dih_i,dih_j,dih_k,dih_l," part of fitting ",ATYPE[dih_i-1],ATYPE[dih_j-1],ATYPE[dih_k-1],ATYPE[dih_l-1]
    
    
				if( add_dih ):
				    if( debug ): print "   adding non loop atoms to DIH_FIT ",DIH[indx_dih]
				    DIH_fit.append( DIH[indx_dih] )
				    DTYPE_IND_fit.append( DTYPE_IND[indx_dih] )
				    
    
			    f.close()
    
	    
			    for dih_indx in range( len(dih_id_list) ):
				dih_id = dih_id_list[dih_indx]
				cent_min = cent_min_list[dih_indx]
				cent_max = cent_max_list[dih_indx]
				cent_step = cent_step_list[dih_indx]
				a_k = a_k_list[dih_indx]
				a_i = a_i_list[dih_indx]
				a_j = a_j_list[dih_indx]
				a_l = a_l_list[dih_indx]

				if( options.verbose ):
				    print "    Writing input files for ",options.ff_software ," torsion: ",cent_min,cent_max,cent_step
				#
				# Append torsional information 
				#
						
				ff_tor_data["cent_id"].append( dih_id )
				ff_tor_data["a_k"].append( a_k )
				ff_tor_data["a_i"].append( a_i )
				ff_tor_data["a_j"].append( a_j )
				ff_tor_data["a_l"].append( a_l )
				ff_tor_data["cent_min"].append( options.cent_min )
				ff_tor_data["cent_max"].append( options.cent_max+options.cent_step )
				ff_tor_data["cent_step"].append( options.cent_step )
				ff_tor_data["ff_type"].append( ff_type_id)
				
						
				
				
				# Loop over angels of central dihedrals
				cent_indx = 0
				for cent_angle in range(cent_min,cent_max,cent_step):
				    ff_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_ff' + options.ff_software + ff_type_id
				    qm_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_auxfix'
		
				    # Check to see if finished
					
				    fchk_file =   qm_id +"/"+qm_id+".fchk"
				    
				    if( options.verbose ):
					print "      Checking ",fchk_file
					
				    if( file_io.file_exists( fchk_file) ):
					if( options.verbose ):
					    print  '        Get results from ', fchk_file
					    
					# Get energy 
					NA, ELN, R, TOTAL_ENERGY, Q_ESP   = gaussian.parse_fchk( fchk_file )
					#		    
					ff_dir =   ff_id 
					if (not os.path.isdir(ff_dir)):
					    os.mkdir(ff_dir)
					#
					print 'making ff data'
					#
					os.chdir(ff_dir)
					#
					
					if( options.ff_software == "gromacs"):
					    #
					    # Print input structure from DFT optimization 
					    #
					    gro_file =  "full.gro"
					    gromacs.print_gro(gro_file,GTYPE,RESID,RESN,R,LV)
    
					    # make ff input
					    
					    #
					    # Print full topology 
					    #
					    top_file = 'full.top'
					    
					    const = []
					    angle = []
					    
					    gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
						   ,const,angle
						   ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
						   ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
						
					    #    
					    #    top_file = 'out_dih.top'
					    #    
					    #    gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
					    #	   ,const,angle
					    #	   ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
					    #	   ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
					    #    
				    
					    #
					    # Print full topology with  constrained dihedral 
					    #
					    DIH_CONST_ANGLE = float( cent_angle )
					    top_file = 'full_const.top'
					    
					    gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
						   ,DIH_CONST_ATOMS[0], DIH_CONST_ANGLE
						   ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
						   ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)

					    #
					    # Print topology without considered dihedral 
					    #					    
					    top_file = 'dih0.top'
					    gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH_fit , IMPS
						   ,const,angle
						   ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
						   ,DTYPE_IND_fit, DIHTYPE_F, IMPTYPE_F,LV)
					    
					    #
					    # Print topology without considered dihedral and constrained dihedral 
					    #					    
					    top_file = 'dih0_const.top'
					    gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH_fit , IMPS
						   ,DIH_CONST_ATOMS[0], DIH_CONST_ANGLE
						   ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
						   ,DTYPE_IND_fit, DIHTYPE_F, IMPTYPE_F,LV)
					    

					    #
					    # Print topology without bonded parameters for LJ and coulombic contributions 
					    #							    
					    top_file = 'out_zero.top'
					    const = []
					    angle = []
					    
					    BONDS_z = []
					    ANGLES_z  = []
					    DIH_z  = []
					    IMPS_z = []
					    
					    gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS_z , ANGLES_z , DIH_z , IMPS_z
						   ,const,angle
						   ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
						   ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
					    
					    #AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(options,ASYMB,ATYPE,BONDS,ANGLES,DIH,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
					    new_itp = "ff-new.itp"
					    AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(new_itp,options.norm_dihparam,ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
    
			
					elif( options.ff_software == "lammps" ):
					    
					    if( options.verbose ):
						print "   Writing lammps input files "
					    
					    last_index = len(DIH_CONST_ANGLE) - 1
					    
					    DIH_CONST_ANGLE[last_index] = float(cent_angle)
    

					    #
					    # Print full topology 
					    #
					    rest_file = 'full_sp.in'
					    data_file = "full.data"
					    lammps.print_sp(rest_file,data_file,options)

					    #
					    # Print full topology with  constrained dihedral 
					    #					    
					    rest_file = 'full_const.in'
					    data_file = "full.data"
					    lammps.print_rest2(rest_file,data_file,DIH_CONST_ANGLE,DIH_CONST_ATOMS,options)

					    #
					    # Print topology without considered dihedral 
					    #					    					    
					    rest_file = 'dih0_sp.in'
					    data_file = "dih0.data"
					    lammps.print_sp(rest_file,data_file,options)

					    #
					    # Print topology without considered dihedral and constrained dihedral 
					    #		
					    rest_file = 'dih0_const.in'
					    data_file = "dih0.data"
					    lammps.print_rest2(rest_file,data_file,DIH_CONST_ANGLE,DIH_CONST_ATOMS,options)
					    
					    # Print lammps structure file 
					    data_file = "dull.data"
					    lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
					       BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
					       ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
					       DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
					       RESN,ATYPE_IND,CHARGES,R , ATYPE,
					       BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV)
					    
					    # Print lammps structure file without considered dihedrals 
					    data_file = "dih0.data"
					    lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
					       BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
					       ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
					       DIH_fit,DTYPE_IND_fit,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
					       RESN,ATYPE_IND,CHARGES,R , ATYPE,
					       BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV)
				    	    
					os.chdir(work_dir)
					os.chdir(struct_dir)
					
				
			    json_file_loc = "%s%s" % ( job_name , ".json" )
			    f = open(json_file_loc	, 'w')
			    json.dump(json_data, f, indent=2)
			    f.close()
			    
			else:
			    print os.getcwd()
			    print " NEED ",dlist_name," which should have been generated by mk_qm.py "
			    sys.exit("no dih list file ")
								    
			os.chdir(work_dir)

if __name__=="__main__":
    main() 

