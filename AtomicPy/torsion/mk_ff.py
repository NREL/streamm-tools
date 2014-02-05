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
    
    # Cluster options
    parser.add_option("--cluster_host", dest="cluster_host",type="string",default="peregrine",help=" name of cluster ")


    parser.set_defaults(submit=False)
    parser.add_option("--submit", dest="submit",action="store_true",help=" submit calculations to the queue ")
    parser.set_defaults(localrun=False)
    parser.add_option("--localrun", dest="localrun",action="store_true",help=" Run calculations locally, ment to use if scripts included in submitted job")
 
    parser.set_defaults(submit_command="qsub")
    parser.add_option("--submit_command", dest="submit_command",type="string",help=" command used to submit script to the queue ")

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

    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")

    (options, args) = parser.parse_args()

    # Set options based on cluster 
    if( options.cluster_host == "peregrine" ):
        options.npros = options.nnodes* 24
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"
    elif( options.cluster_host == "redmesa" ):
        options.npros = options.nnodes *8
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/g09/C.01"
	    
	    
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
    
def write_esp_com(calc_id,ASYMB,R):
    
    esp_com = calc_id +  ".com"
    
    f = file(esp_com, "w")
    f.write(  '%chk='+str(calc_id)+'.chk' + '\n' )
    #lines_opt =  '#P B3LYP/6-31+g**  OPT nosym '
    lines_opt =  '#P HF/3-21g  SP  POP=MK'
    lines_opt = lines_opt + '\n' + ' '
    lines_opt = lines_opt + '\n' + str(calc_id) + " esp fit "
    lines_opt = lines_opt + '\n' + ' ' 
    lines_opt = lines_opt + '\n' + '0 1 ' + '\n' 
    f.write(lines_opt)

    for atom_i in range( len(ASYMB) ):
        f.write( " %5s %12.6f %12.6f %12.6f \n" % ( ASYMB[atom_i],R[atom_i][0], R[atom_i][1], R[atom_i][2] ))
        
    f.write( ' ' + '\n')
    f.close()
    
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
    import string  , numpy 
    from string import replace
    import file_io, gaussian, elements, gromacs, lammps , top, atom_types, xmol 
    
    options, args = get_options()
    
    #
    if( options.cluster_host == "peregrine" ):
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
    
    elif( options.cluster_host == "redmesa" ):
	load_gaussian = "module load gaussian/g09/C.01"
        user = 'twkempe'
	
	if( options.ff_software == "gromacs"):
	    md_min_template = "redmesa.md_min.template"
	    f = open(md_min_template,"r")
	    md_min_templ = f.read()
	    f.close()
	elif( options.ff_software == "lammps" ):	
	    lammps_min_template = "redmesa.lmp_min.template"
	    f = open(lammps_min_template,'r')
	    md_min_templ = f.read()
	    f.close()
	
	    
	f = open("redmesa.slurm.template",'r')
	pbs_templ = f.read()
	f.close()

    elif( options.cluster_host == "macbook" ):
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
    
    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
        for gen_line in Lines:
            col = gen_line.split()
            if( len(col) >= 4 and col[0] == "gen" ):
		
                mol_dir = col[1].strip()
                mol_id = col[2].strip()
                mol_repeat = int(col[3].strip() )
                mol_acc = col[4].strip()
                
		
                # File info
                struct_dir = mol_dir + "/" + mol_id + "/"
                job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)
		
		
		print "   Running job ",struct_dir,job_name

		# Get geometry from ZMAT
                calc_id = job_name + "-ZMAT"
		log_file = struct_dir +'/' + calc_id +"/"+calc_id+".log"
		fchk_file = struct_dir +'/' + calc_id +"/"+calc_id+".fchk"
		    
		zmat_finished = 0 
		esp_finished = 0 
		if( file_io.file_exists( fchk_file) ):
                    if( options.verbose ):
			print  '    Get results from ', fchk_file
			    
		    # Read in from zmatrix optimization 
		    NA, ELN, R, TOTAL_ENERGY, Q_ESP   = gaussian.parse_fchk( fchk_file )
		    ASYMB = elements.eln_asymb(ELN)
		    
		    # Get bonds from log file 
		    BONDS, ANGLES, DIH  = gaussian.read_optlog(log_file)
		    NBLIST,NBINDEX = build_nablist(ELN,BONDS)
		    zmat_finished = 1
    
		    # pars zmatirix file
		    comz_name = struct_dir+'/'+ job_name+'-ZMAT/'+ job_name +  "-ZMAT.com"
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
    
				
		else:
		    print "  could not find ",fchk_file

		    
		if( zmat_finished ):
		    # Run esp fit
		    calc_id_esp = job_name + "-ESP"
		    input_file = calc_id_esp + ".com"
		    log_file_esp = struct_dir +'/' + calc_id_esp +"/"+calc_id_esp+".log"
		    fchk_file_esp = struct_dir +'/' + calc_id_esp +"/"+calc_id_esp+".fchk"
		    
		    # Run esp fit if not finished
		    if( not file_io.file_exists( fchk_file_esp) ):
                        # Make esp input file
			if( options.verbose ):
			    print  '    Writing esp fit input file ',calc_id_esp," in ",struct_dir
			
					    
			os.chdir(struct_dir)
			write_esp_com(calc_id_esp,ASYMB,R)
    
                        if( options.submit ):
                            if( options.verbose ):
                                print "     Submitting ESP optimization to queue "
                            # Print pbs script
                            pbs_id = calc_id_esp+'.pbs'
                            pbs_name = pbs_id
                            pbs_dih = pbs_templ
                            pbs_dih = replace(pbs_dih,"<calc_id>",calc_id_esp)
                            pbs_dih = replace(pbs_dih,"<input_file>",input_file)
                            pbs_dih = replace(pbs_dih,"<nnodes>",str( options.nnodes) )
                            pbs_dih = replace(pbs_dih,"<pmem>",str( options.pmem) )
                            pbs_dih = replace(pbs_dih,"<npros>",str( options.npros) )
                            f = file(pbs_name, "w")
                            f.write(pbs_dih)
                            f.close()
                            # Submit job
                            submit_job( struct_dir, pbs_id ,options )
                            
                        elif( options.localrun ):
                            if( options.verbose ):
				print  '    Running ESP fit '
                            gaussian.run(options, calc_id_esp)
			    
                        else:
                            print " Please mark either qsub or runloc options as True to run qm"	
			
			os.chdir(work_dir)
			    
		    # Check if finished
		    #   basicaly if local run is done or should just skip over if submitted to queue
		    if(  file_io.file_exists( fchk_file_esp) ):
			if( options.verbose ):
			    print "    Getting charges from esp fit ",fchk_file_esp
			NA, ELN, R, TOTAL_ENERGY , CHARGES  = gaussian.parse_fchk( fchk_file_esp )
			    
			esp_finished = 1
			    

		dlist_exists = 0 
		if( zmat_finished and  esp_finished ):
		    

		    os.chdir(struct_dir)

                    dlist_name = job_name + "_dih.list"

                    dlist_exists = file_io.file_exists( dlist_name )
                    if ( dlist_exists ):                    
                        if( options.verbose ):
                            print "       Writing input files for all loop dihedrals "
                            print "         Nodes ",options.nnodes
                            print "         Processors  ",options.npros
                        DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS = read_dihlist(dlist_name)
    			    
		    else:
			print os.getcwd()
			print " NEED ",dlist_name," which should have been generated by mk_qm.py "
			sys.exit("no dih list file ")
			
		    os.chdir(work_dir)
		    
		os.chdir(work_dir)
			
		if ( dlist_exists ):    
		    ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
		    DIH = top.nblist_dih(NA,NBLIST, NBINDEX)
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
		    os.chdir(struct_dir)
		    
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
		
		
                    dlist_name = job_name + "_dih.list"

                    dlist_exists = file_io.file_exists( dlist_name )
                    if ( dlist_exists ):                    
                        if( options.verbose ):
                            print "       Writing input files for all loop dihedrals "
                            print "         Nodes ",options.nnodes
                            print "         Processors  ",options.npros
                        DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS = read_dihlist(dlist_name)
    
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

		    os.chdir(work_dir)

		    # Loop of dihedrals calculated at the qm level
		    for qm_dih_line in Lines:
			qm_dih_col = qm_dih_line.split()
			
			
			if( len(qm_dih_col) >= 4 and qm_dih_col[0] == "qm_dih" and qm_dih_col[2] == mol_id.strip() ):
			    
			    dih_id = qm_dih_col[5].strip() 
			    a_k  = int(  qm_dih_col[6].strip() )
			    a_i  = int(  qm_dih_col[7].strip() )
			    a_j  = int(  qm_dih_col[8].strip() )
			    a_l  = int(  qm_dih_col[9].strip() )
			    cent_min  = int(  qm_dih_col[10].strip() )
			    cent_max = int(  qm_dih_col[11].strip() )
			    cent_step = int(  qm_dih_col[12].strip() )
			    
			    if( options.verbose ):
				print "    Writing input files for ",options.ff_software ," torsion: ",cent_min,cent_max,cent_step

			    # Print calculation information     
			    file_info = open( indx_file,'a')
			    file_info.write( "\n ff_dih  %s %s  %s %s %s %8d %8d %8d %8d %d %d %d " %( struct_dir,job_name,options.ff_software, ff_type_id,dih_id, a_k,a_i, a_j ,a_l,cent_min,cent_max,cent_step))
			    file_info.close()				

			    # Loop over angels of central dihedrals
			    cent_indx = 0
			    for cent_angle in range(cent_min,cent_max,cent_step):
			        ff_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_ff' + options.ff_software + ff_type_id
			        qm_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_auxfix'
	    
				# Check to see if finished
				    
				fchk_file = struct_dir +  qm_id +"/"+qm_id+".fchk"
				
				if( options.verbose ):
				    print "      Checking ",fchk_file
				    
				if( file_io.file_exists( fchk_file) ):
				    if( options.verbose ):
					print  '        Get results from ', fchk_file
					
				    # Get energy 
				    NA, ELN, R, TOTAL_ENERGY, Q_ESP   = gaussian.parse_fchk( fchk_file )
							
				    ff_dir =  struct_dir + ff_id 
				    if (not os.path.isdir(ff_dir)):
					os.mkdir(ff_dir)
				    
				    print 'making ff data'
				    
				    os.chdir(ff_dir)

				    if( options.ff_software == "gromacs"):

					gro_file =  "out.gro"
					gromacs.print_gro(gro_file,GTYPE,RESID,RESN,R,LV)

					# make ff input
					top_file = 'out.top'
					
					const = []
					angle = []
					
					gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
					       ,const,angle
					       ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
					       ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
					
					top_file = 'out_dih.top'
					
					gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
					       ,const,angle
					       ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
					       ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
					
					DIH_CONST_ANGLE = float( cent_angle )
			
					top_file = 'out_const.top'
					gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
					       ,DIH_ATOMS[cent_indx], DIH_CONST_ANGLE
					       ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
					       ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
					
					top_file = 'out_fit.top'
					gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH_fit , IMPS
					       ,DIH_ATOMS[cent_indx], DIH_CONST_ANGLE
					       ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
					       ,DTYPE_IND_fit, DIHTYPE_F, IMPTYPE_F,LV)
					
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
					
					last_index = len(DIH_CONST_ANGLE) - 1
					
					DIH_CONST_ANGLE[last_index] = float(cent_angle)

					# Print lammps input file with dihedral restraint 
					rest_file = 'full_sp.in'
					data_file = "out.data"
					lammps.print_sp(rest_file,data_file,options)
					
					rest_file = 'd0_sp.in'
					data_file = "fit.data"
					lammps.print_sp(rest_file,data_file,options)
					
					rest_file = 'rest.in'
					data_file = "out.data"
					lammps.print_rest2(rest_file,data_file,DIH_CONST_ANGLE,DIH_CONST_ATOMS,options)
					
					rest_file = 'fit.in'
					data_file = "fit.data"
					lammps.print_rest2(rest_file,data_file,DIH_CONST_ANGLE,DIH_CONST_ATOMS,options)
					
					# Print lammps structure file 
					data_file = "out.data"
					lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
					   BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
					   ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
					   DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
					   RESN,ATYPE_IND,CHARGES,R , ATYPE,
					   BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV)
					
					# Print lammps structure file without considered dihedrals 
					data_file = "fit.data"
					lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
					   BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
					   ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
					   DIH_fit,DTYPE_IND_fit,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
					   RESN,ATYPE_IND,CHARGES,R , ATYPE,
					   BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV)
					
				    os.chdir(work_dir)
				
		os.chdir(work_dir)

if __name__=="__main__":
    main() 

