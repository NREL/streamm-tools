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
    
    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    

    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")

    (options, args) = parser.parse_args()

    # Set options based on cluster 
    if( options.cluster_host == "peregrine" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"
    elif( options.cluster_host == "redmesa" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/g09/C.01"
	    
    return options, args


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


    
def main():
    import sys, os 
    import string  , numpy  , json 
    from string import replace
    import file_io, gaussian, elements, gromacs, lammps , top, atom_types, xmol 
    
    options, args = get_options()
    
    if( options.cluster_host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	
    elif( options.cluster_host == "redmesa" ):
	load_gaussian = "module load gaussian/g09/C.01"
        user = 'twkempe'
	

    elif( options.cluster_host == "macbook" ):
	load_gaussian = ""
        user = 'tkemper'


	
    # sufix for force field calcs to run multiple ff types or excultions, such as q(i) = 0.00 
    ff_type_id = "_fit"
    
    LAT_CONST = numpy.zeros( (3,3) )
    
    LAT_CONST[0][0] = 50.0
    LAT_CONST[1][1] = 50.0
    LAT_CONST[2][2] = 50.0

    # Store working dir  
    work_dir = os.getcwd()
    
    # Read index files from args
    for indx_file in args:
        # Get lines of record file   
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

		# Get geometry from optimized structure 
                calc_id = job_name
		log_file = struct_dir +'/' + calc_id +"/"+calc_id+".log"
		fchk_file = struct_dir +'/' + calc_id +"/"+calc_id+".fchk"
		
		opt_finished = 0 
		esp_finished = 0 
		if( file_io.file_exists( fchk_file) ):
                    if( options.verbose ):
			print  '    Get results from ', fchk_file
			    
		    # Read in from zmatrix optimization 
		    NA, ELN, R, TOTAL_ENERGY, Q_ESP   = gaussian.parse_fchk( fchk_file )
		    ASYMB = elements.eln_asymb(ELN)
		    
		    #   Build covalent nieghbor list for bonded information 
		    NBLIST, NBINDEX = top.build_covnablist(ELN,R)
		    opt_finished = 1
    
		    # pars zmatirix file
		    #comz_name = struct_dir+'/'+ job_name+'-ZMAT/'+ job_name +  "-ZMAT.com"
		    #f = open(comz_name,'r')
		    #zmatrix_lines = f.readlines()
		    #f.close()
		    
		    # need to change to reading out of connections 
		    #DIH_ID,DIH_ATOMS = get_dih_id( RING_CONNECT,RING_NUMB, zmatrix_lines )
		    # Find rings		    
		    RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
		    
		    RING_CONNECT  = top.find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB) 

		    #DIH_ID,DIH_ATOMS,zmatrix = gaussian.get_dih_id( RING_CONNECT,RING_NUMB, zmatrix_lines )
		    
		    #zmatrix = gaussian.com_zmatrix(comz_name)    
		    #DIH_ID, DIH_VAL, DIH_ATOMS = gaussian.get_dih_id( zmatrix)
    
				
		else:
		    print "  could not find ",fchk_file

		    
		if( opt_finished ):
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
			
                        else:
			    print " esp fit results do  not exist "
                            print " Please mark either qsub or runloc options as True to run qm"	
			    sys.exit(" Need esp fit ")
			    
			os.chdir(work_dir)
			    
		    # Check if finished
		    #   basicaly if local run is done or should just skip over if submitted to queue
		    if(  file_io.file_exists( fchk_file_esp) ):
			if( options.verbose ):
			    print "    Getting charges from esp fit ",fchk_file_esp
			NA, ELN, R, TOTAL_ENERGY , CHARGES  = gaussian.parse_fchk( fchk_file_esp )
			ASYMB = elements.eln_asymb(ELN)
			
			debug = 0
			if( debug):
			    print len(ELN)
			    print ""
			    for i in range( len(ELN)):
				print ELN[i],R[i][0],R[i][1],R[i][2]
			    
			    sys.exit("debug pol_ff R")
			    
			    
			esp_finished = 1
			
		    if(  esp_finished ):
			    
			# Read data from json file
			
				    
			json_name = struct_dir +"/" + job_name +".json"
			if(  file_io.file_exists( json_name) ):
			    if( options.verbose ):
				print "    Getting atomic data from  ",json_name
						    
			    f = open(json_name, 'r')
			    json_data = json.load(f)
			    f.close()
			    
			    
			    
			    
			    print json_data['metadata']["atomic"]["atype"] 
			    
			    sys.exit("gen_struct json read debug  ")
			

	    		# Update atomic information	
			CHARN = top.set_chargegroups(options.options,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LAT_CONST)
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
			# This step has be trasfered to donoracceptorsystems py 
			#CHARGES = top.zero_termq(ELN,ATYPE,CHARGES,tag2,NBINDEX,NBLIST,options.verbose)
			
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
				tag1[atom_i] = "func_C(" + str(i_n) + ")"
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
					tag1[atom_j] = "funccap_H(" + str(j_n) + ")_on_C("+str(i_n)+")"
					
					print " func found ", tag1[atom_j] 
				
				if( term_con_cnt != 1 ):
				    print " Number of functional atoms attached to atom ",i_n," not equal to 1 "
				    sys.exit(" Error in functional connections ")

			
						
			    
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
			F.write('D(R1)')
			for atom_i in range( len(ELN) ):
			    F.write('\n %s %16.6f %16.6f %16.6f %s ' % (ASYMB[atom_i],R[atom_i][0],R[atom_i][1],R[atom_i][2],str( tag1[atom_i]) ) )			
			F.close()
			

			    
		os.chdir(work_dir)

if __name__=="__main__":
    main() 

