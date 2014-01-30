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
    

    parser.add_option("-v","--verbose", dest="verbose", default=True,action="store_true", help="Verbose output ")
    
    # Cluster options
    parser.add_option("--cluster_host", dest="cluster_host",type="string",default="peregrine",help=" name of cluster ")


    parser.set_defaults(submit=False)
    parser.add_option("--submit", dest="submit",action="store_true",help=" submit calculations to the queue ")
    parser.set_defaults(localrun=False)
    parser.add_option("--localrun", dest="localrun",action="store_true",help=" Run calculations locally, ment to use if scripts included in submitted job")
 
    parser.set_defaults(submit_command="qsub")
    parser.add_option("--submit_command", dest="submit_command",type="string",help=" command used to submit script to the queue ")

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
    
    parser.set_defaults(norm_dihparam=False)
    parser.add_option("--norm_dihparam", dest="norm_dihparam",action="store_true",help=" divide dihedral parameters in ff itp file by 4 ")

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
    
    if( options.cluster_host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	
	
	if( options.ff_software == "gromacs"):
	    md_min_template = "peregrine.md_min.template"
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
    
    LAT_CONST = numpy.zeros( (3,3) )
    
    LAT_CONST[0][0] = 50.0
    LAT_CONST[1][1] = 50.0
    LAT_CONST[2][2] = 50.0

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
                            gaussian.run1(options, calc_id_esp)
			    
                        else:
                            print " Please mark either qsub or runloc options as True to run qm"	
			
			os.chdir(work_dir)
			    
		    # Check if finished
		    #   basicaly if local run is done or should just skip over if submitted to queue
		    if(  file_io.file_exists( fchk_file_esp) ):
			if( options.verbose ):
			    print "    Getting charges from esp fit ",fchk_file_esp
			NA, ELN, R, TOTAL_ENERGY , CHARGES  = gaussian.parse_fchk( fchk_file_esp )
			
			debug = 0
			if( debug):
			    print len(ELN)
			    print ""
			    for i in range( len(ELN)):
				print ELN[i],R[i][0],R[i][1],R[i][2]
			    
			    sys.exit("debug pol_ff R")
			    
			    
			esp_finished = 1
			
		    if(  esp_finished ):
			    
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
			CHARN = top.set_chargegroups(options,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LAT_CONST)
    
			
			debug = 0
			if( debug):
			    for i in range( len(ELN) ):
				print i,ELN[i],ASYMB[i],ATYPE[i],CHARGES[i] # , CHARGES[i]
				
			    sys.exit( 'mk_ff')
			    
			# Read in units
			
			calc_id_esp = job_name + "-ESP"
			unit_file = struct_dir +'/' +job_name+".unit"
			F = open(unit_file , 'r' )
			Lines = F.readlines()
			F.close()
			
					    
			# Read in unit file
			#
			# unit = unit #
			# tag1 = ?
			# tag2
			#       U - central unitl
			#       T - Atom attached to terminal / connector 
			#       X - Atom to be replaced when another monomer is added 
			#       R - Functinal group atom
			#       F - Atom attached to Functinal group atom
			tag1 = []
			tag2 = []
			unit = []
			line_cnt = 0
			for line in Lines :
			    col = line.split()
			    if( len(col) > 3  ):
				#print col[4]
				u = int( col[4] ) 
				unit.append( u )
				tag1.append( "" )
				tag2.append( col[6] )
				
				
			
			unit_max = max( unit )
			print "  unit_max ",unit_max
			unit_q = numpy.zeros( unit_max )
			unit_cnt = numpy.zeros( unit_max )
			unit_q_exs = numpy.zeros( unit_max )
			#unit_cnt_exs = numpy.zeros( unit_max )
			
			
			
			for u in range(unit_max):
			    for i in range( len(ELN) ):
				if( u ==  unit[i] -1  ):
				    unit_q[u] += CHARGES[i]
				    
				    if(  tag2[i].strip() == 'U' ):
					unit_cnt[u] += 1
				    else:
					unit_q_exs[u] += CHARGES[i]
					
				    #print i,unit[i] ,u,CHARGES[i],unit_q[u]

			#sys.exit( 'mk_ff 1 - 2 ' )
				 
			n_chrgrp =  max( CHARN )
			q_chrgrp = numpy.zeros( n_chrgrp )
			print "   n_chrgrp ",n_chrgrp
			
			for i in range( len(ELN) ):
			    print i,ELN[i],ASYMB[i],ATYPE[i],CHARGES[i] ,CHARN[i], unit[i],tag1[i],tag2[i]
			    q_n = CHARN[i] - 1
			    q_chrgrp[q_n] += CHARGES[i]
			
			q_exs_u = numpy.zeros( unit_max )
			for u in range(unit_max):
			    q_exs_u[u] = unit_q[u] - unit_q_exs[u]
			    print u , " ",unit_q[u] , unit_q_exs[u]
			    print "      excess per atom ",q_exs_u[u]/unit_cnt[u] 
			
			
			#
			# Sum terminal hydrogen charges into carbon atoms 
			#
			CHARGES = top.zero_termq(ELN,ATYPE,CHARGES,tag2,NBINDEX,NBLIST,options.verbose)
			
			#
			# Set tag1 for cply output
			#
			if( options.verbose): print "    setting cply tags "
			for atom_i in range( len(ELN) ):
			    i_n = atom_i + 1
			    # Set terminal attached carbons 
			    if( tag2[atom_i] == "T" and ELN[atom_i] == 6 ):
				tag1[atom_i] = "term_C(" + str(i_n) + ")"
				#
				term_con_cnt = 0 # Terminal connections count
				#
				# Find terminal hydrogen
				#		    
				N_o = NBINDEX[atom_i]
				N_f = NBINDEX[atom_i+1] - 1
				for j_indx in range( N_o,N_f+1):
				    atom_j = NBLIST[j_indx]
				    j_n = atom_j + 1
				    if( tag2[atom_j] == "X" and ELN[atom_j] == 1 ):
					term_con_cnt += 1
					tag1[atom_j] = "termcap_H(" + str(j_n) + ")_on_C("+str(i_n)+")"
					
				if( term_con_cnt > 1 ):
				    print " Number of terminal atoms attached to atom ",i_n," greater than 1 "
				    sys.exit(" Error in terminal connections ")

			
			    # Set functional attached carbons 
			    if( tag2[atom_i] == "F" and ELN[atom_i] == 6 ):
				tag1[atom_i] = "f_C(" + str(i_n) + ")"
				#
				term_con_cnt = 0 # Terminal connections count
				#
				# Find function hydrogen
				#
				
				
				
				N_o = NBINDEX[atom_i]
				N_f = NBINDEX[atom_i+1] - 1
				for j_indx in range( N_o,N_f+1):
				    atom_j = NBLIST[j_indx]
				    j_n = atom_j + 1
				    print " ",atom_j,j_n,tag2[atom_j],ELN[atom_j]
				    if( tag2[atom_j] == "R" and ELN[atom_j] == 1 ):
					term_con_cnt += 1
					tag1[atom_j] = "fcap_H(" + str(j_n) + ")_on_C("+str(i_n)+")"
					
					print " func found ", tag1[atom_j] 
				
				if( term_con_cnt != 1 ):
				    print " Number of functional atoms attached to atom ",i_n," not equal to 1 "
				    sys.exit(" Error in functional connections ")

			
			
    
			# Identify total number of atom types for lammps output 
			ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH)
    
			# Check atom types to be sure each atom of the same type has the same number of neighbors 
			ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)
		    
			# Modify dih potentials
			## FF_DIHTYPES = mod_dih( ATYPE_IND
			
			
			# Print new itp file with only used atom types and interactions
			if( options.ff_software == "gromacs"):
			    #AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(options,ASYMB,ATYPE,BONDS,ANGLES,DIH,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
			    AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)

			# Print output
			os.chdir(struct_dir)
			
			out_xyz =  job_name + "-ff.xyz"
			xmol.write_xyz(ASYMB,R,out_xyz)
			xmol.print_cgxyz(CHARN,R)
			job_gro = job_name +  ".gro"
			gromacs.print_gro(job_gro,GTYPE,RESID,RESN,R,LAT_CONST)
			
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
					   ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LAT_CONST)
    
			    AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(options,ASYMB,ATYPE,BONDS,ANGLES,DIH,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
			
					    
			    
			elif( options.ff_software == "lammps" ):
			    
			    data_file = "out.data"
			    
			    lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
				  BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
				  ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
				  DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
				  RESN,ATYPE_IND,CHARGES,R , ATYPE,
				  BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LAT_CONST)
			    
			#
			# Print new building block with charges 
			#
			bb_dir = "BuildingBlock_local"
			if (not os.path.isdir(bb_dir)):
			    os.mkdir(bb_dir)
			    
			for bb_subdir in ( "acceptors","donors","functional_groups","spacers","terminals" ):
			    mk_dir  = bb_dir +"/"+bb_subdir
			    if (not os.path.isdir(mk_dir)):
			        os.mkdir(mk_dir)
			
			    
			bb_dir = "BuildingBlock_local/donors"
			    
			os.chdir(bb_dir)

			bb_file = job_name + ".cply"
			F = open(bb_file,'w')
			F.write('D(f1)')
			for atom_i in range( len(ELN) ):
			    print ASYMB[atom_i],R[atom_i][0],R[atom_i][1],R[atom_i][2],CHARGES[atom_i],tag1[atom_i],tag2[atom_i]
			    F.write('\n %s %16.6f %16.6f %16.6f %16.6f %s %s ' % (ASYMB[atom_i],R[atom_i][0],R[atom_i][1],R[atom_i][2],CHARGES[atom_i],str( tag2[atom_i]),str( tag1[atom_i]) ) )
			
			F.close()
			

		os.chdir(work_dir)

if __name__=="__main__":
    main() 

