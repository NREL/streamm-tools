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
    
    parser.add_option("-r","--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")

    
    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")
    parser.add_option("--pbs_template", dest="pbs_template",type="string",default="",help=" Template for job submission  ")

    # How to run the needed calculations 
    parser.add_option("--submit", dest="submit",action="store_true", default=False,help=" submit calculations to the queue ")
    parser.add_option("--localrun", dest="localrun",action="store_true", default=False,help=" Run calculations locally")
    parser.add_option("--submit_command", dest="submit_command",type="string", default="qsub",help=" command used to submit script to the queue ")

    
    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    parser.add_option("--qm_load", dest="qm_load",type="string",help=" string to load qm software module  ")

    parser.add_option("--qm_method", dest="qm_method", type="string",default="B3LYP", help="Method of QM calculation ")
    parser.add_option("--qm_basis", dest="qm_basis", type="string",default="6-31G**", help="Basis set of QM calculation ")
    parser.add_option("--qm_kywd", dest="qm_kywd", type="string",default="", help="Key words for QM calculation ")
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    parser.add_option("--qm_charge", type="int",action="append", default="0",help="Input gaussain log file ")
    parser.add_option("--qm_mult", dest="qm_mult", type="int",default="0", help=" Shift in default spin multiplicity ( singlet,doublet) QM calculation, allows for triplets ")
    
    parser.add_option("--pmem",dest="pmem", type="int", default="1700",help=" Memory per processor ")
    parser.add_option("--npros", dest="npros", type="int", default="4",help=" Number of processors ")
    parser.add_option("--nnodes", dest="nnodes", type="int", default="1",help=" Number of nodes ")

    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")
    parser.add_option("-s","--check_dihsym", dest="check_dihsym", default=True,action="store_true", help=" Prin inter-ring dihedral values")
    parser.add_option("-m","--mult_sym", dest="mult_sym", default=1,type=float,help=" value to multiply symetric dihedrals ")
    parser.add_option("-z","--sym_zmat", dest="sym_zmat", default=False,action="store_true", help=" Print symetric  zmatirix in input file")

    (options, args) = parser.parse_args()

    
    if( options.host == "peregrine" and options.submit ):
	options.pbs_template = "peregrine.pbs.template"
	print "  module load gaussian/.g09_C.01"
	
    
    if( options.host == "dale" and options.submit ):
	options.pbs_template = "dale.pbs.template"
	print "  module load gaussian/.g09_C.01"
	options.npros = 8
	
	    
    return options, args


def main():
    import string, os , sys 
    # atomicpy
    import jsonapy
    import gaussian, elements, xmol , file_io , cluster , top 
    from string import replace
    
    #
    # Set some defaults 
    #
    default_method = 'b3lyp'
    default_basis = '6-31G**'
    
    options, args = get_options()
    
    if(options.submit ):
	f = open(options.pbs_template,'r') # redmesa.slurm.template
	pbs_templ = f.read()
	f.close()

    work_dir = os.getcwd()

	
    json_files = options.json.split(',')
    print json_files
    if( len(json_files) > 0 ):
	# Read index files from args
	for json_file in json_files:
	    # Get lines of index file
		
	    
	    # Verbose output
	    if( options.verbose ):
		print "The molecules specified in json file ",options.json," will be read in "
    
	    json_data,json_success = jsonapy.read_jsondata(json_file)
	    if(  json_success ):
		#
		mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found ):
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
		    calc_id = "%s/%s%s" % (struct_dir, job_name , "-ESP" )
		    #
		    # 
		    #
		    fchk_file = "%s/%s%s" % ( calc_id ,job_name,"-ESP.fchk" )
		    print " Checking for complete esp fit ",fchk_file
                    calc_finished = file_io.file_exists( fchk_file )
		    
		    if( not calc_finished or options.recalc ):
			#
			# Asign method and basis 
			#      	
			if( len( options.qm_method) == 0 ):
			    options.qm_method = method		    
			if( len( options.qm_basis) == 0 ):
			    options.qm_basis = basis
			#
			# Get atomic data 
			#      
			if( options.verbose ):
			    print "     Getting atomic data from  ",json_file
	    
			ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,json_atomicdata  = jsonapy.read_atomic(json_data)
			
			if( not json_atomicdata ):
			    print "   json file ",json_file," exist, but does not contain any atomic data . "
			#
			# If get optimized atomic data from fchk file 
			#
			print " checking fchk files in ",struct_dir, " job name ",job_name+ '-ZMATOPT/' 
			
			fchk_file = struct_dir + job_name + '-ZMATOPT/' + job_name + "-ZMATOPT.fchk"
			
			fchk_atomicdata = 0
			if( file_io.file_exists(fchk_file) ):
			    if( options.verbose ):
				print " Reading atomic data from ",fchk_file
			    NA, ELN, R, TOTAL_ENERGY , Q_FCHK = gaussian.parse_fchk(fchk_file)
			    ASYMB = elements.eln_asymb(ELN)
			    #
			    # 
			    #
			    CHARGES = []				
			    ATYPE = []
			    ELECTRONS_i = 0
			    for atom_i in range( len(ELN) ):
				ATYPE.append( ASYMB[atom_i] )
				ELECTRONS_i += ELN[atom_i] 
				CHARGES.append( -100.0 )
			    
		    
			    # Optimize z-matrix to get bonding information
			    qm_kywd_o = options.qm_kywd 
			    options.qm_kywd = qm_kywd_o + " OPT=(tight) POP=MK FREQ "
			    #
			    # Print calculation information 
			    #
			    if( options.verbose ):
				print " Atomic data has been found creating zmatrix input files "
				print "   Method   ",options.qm_method
				print "   Basis    ",options.qm_basis
				print "   Keywords ",options.qm_kywd
	    
			    if( options.out_xyz ):
				if( options.verbose ):
				    print "      Writing xyz file of fchk geometry ",options.out_xyz
				xmol.write_xyz(ASYMB,R,options.out_xyz)
		    
			    os.chdir(struct_dir)
			    job_id =  "%s%s" % (job_name,"-ESP" )
			    com_name = job_id +".com"
			    
			    # Print com
			    if( not file_io.file_exists(com_name) ):
				if( options.verbose ):
				    print "       Printing temp zmatrix "
				    
				calc_id_temp = job_id + "-temp"
				gaussian.print_com( calc_id_temp, ASYMB,R,ATYPE,CHARGES,ELECTRONS_i,options.qm_method,options.qm_basis,options.qm_kywd,options.qm_charge,options.qm_mult)
				#  geometry z-matrix opt input files
				gaussian.com2zmat(calc_id_temp,job_id,options)
				
			    # Find dihedrals to make sure molecule is symetric
			    if( options.check_dihsym ):
							
				# Parse log file 
				NBLIST,NBINDEX = top.build_covnablist(ELN,R)
				BONDS  = top.nblist_bonds(NA,NBLIST, NBINDEX)
				
				# Find rings 
				RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
				RING_CONNECT  = top.find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB) 
			    
				N_CONECT = len(RING_CONNECT)
				if( options.verbose):
				    print "   Molecule ",job_name, " has ",max(RING_NUMB)," rings with ",N_CONECT," inter ring connections "
			    
				zmatrix = gaussian.com_zmatrix(com_name)    
				DIH_ID, DIH_VAL, DIH_ATOMS = gaussian.get_dih_id( zmatrix)
				    
				conect_indx = -1
				half_nconect = int(N_CONECT/2)
				
				#mod_connect =  N_CONECT%2
				#print " mod_connect ",mod_connect
				#if( mod_connect == 0 ):
				#    if( options.verbose): " An even # of connections "
				#   mult_sym = -1
				#else:
				#   mult_sym = 1
				#   if( options.verbose): " An odd  # of connections "
					
				CONECT_DIH = []
				CONECT_IND = []
				for dih_indx in  range(len(DIH_ID)):
				
				    zmat_l = DIH_ATOMS[dih_indx][0]
				    zmat_i = DIH_ATOMS[dih_indx][1]
				    zmat_j = DIH_ATOMS[dih_indx][2]
				    zmat_k = DIH_ATOMS[dih_indx][3]
				    
				    # Make sure inter ring connect not improper 
				    if (  RING_NUMB[zmat_i] != RING_NUMB[zmat_j] and RING_NUMB[zmat_l] != RING_NUMB[zmat_k] ):
					conect_indx += 1
					CONECT_DIH.append(DIH_VAL[dih_indx] )
					CONECT_IND.append(dih_indx )
					print "   Dihderal ",dih_indx,zmat_i, zmat_j,"  between ring ",RING_NUMB[zmat_i] , RING_NUMB[zmat_j] ," with unit #'s ", UNITNUMB[zmat_i], UNITNUMB[zmat_j]," = " ,DIH_VAL[dih_indx]
					
				    
				if( options.sym_zmat ):
    
				    for conect_indx in range(half_nconect):
					sym_connect = N_CONECT - conect_indx - 1
					print conect_indx,sym_connect
					print "   connection ",conect_indx,CONECT_DIH[conect_indx]," with sym ",sym_connect,CONECT_DIH[sym_connect]
					dih_ave = ( CONECT_DIH[conect_indx]  + options.mult_sym*CONECT_DIH[sym_connect])/2.0
					CONECT_DIH[conect_indx] = dih_ave
					CONECT_DIH[sym_connect] = options.mult_sym*dih_ave
					print "    Average dihdral ",dih_ave
				    
				    # Update dih values
				    for conect_indx in range(N_CONECT):
					dih_indx = CONECT_IND[conect_indx]
					DIH_VAL[dih_indx] = CONECT_DIH[conect_indx] 
					print "   connection ",conect_indx," is dih ",dih_indx," in zmat with new value ",DIH_VAL[dih_indx]
				    # Update zmatrix
				    am_opt = '%chk='+str(job_id)+'.chk' 
				    am_opt = am_opt + '\n%nproc='+str(options.npros)
				    #am_opt = am_opt + '\n' + '#P AM1/3-21g  popt=Zmat nosym '
				    keywd = str(options.qm_method) +"/"+ str(options.qm_basis) +" "+ str(options.qm_kywd)
				    am_opt = am_opt + "\n #P " + keywd
				    #f.write( "\n# P %s/%s  %s" % (opt:qions.qm_method,options.qm_basis,options.qm_kywd))
				    am_opt = am_opt + '\n' + ' '
				    am_opt = am_opt + '\n' + job_id
				    am_opt = am_opt + '\n' 
				    am_opt = am_opt + '\n 0 1 '
				    
				    # Prin all non constranied elments of zmatrix 
				    for line in iter(zmatrix.splitlines()) :
					print_z = 1
					colvar = line.split('=')
					var_id = colvar[0].strip()
					for dih_indx in  range(len(DIH_ID)):
					    if( DIH_ID[dih_indx] == var_id ):
						line = " %s = %f " % (var_id,DIH_VAL[dih_indx])
						break
					    
					if( print_z ):
					    am_opt =  am_opt+"\n" + line
									
				    am_opt = am_opt + ' \n'
				    am_opt = am_opt + ' \n'
				    
				    f = file(com_name, "w")
				    f.write(am_opt)
				    f.close()
			    # 
			    # Run optimization
			    # 
			    if( options.submit ):
				input_file =  "%s%s" % ( job_id,".com" )
				if( options.verbose ):
				    print "     Submitting Z-matrix optimization to queue "
				    print "         nodes        ",options.nnodes
				    print "         memory/node  ",options.pmem
				# Print pbs script
				pbs_id = cluster.write_pbs(pbs_templ,job_id,input_file,options)
				cluster.submit_job( struct_dir, pbs_id ,options )
				
			    elif( options.localrun ):
				if( options.verbose ):
				    print "       Running Z-matrix optimization  "
				gaussian.run(options, job_id)
				
			    else:
				print " Please mark either qsub or runloc options as True to run qm"
				
			    os.chdir(work_dir)

			    options.qm_kywd = qm_kywd_o

			else:
			    print " no atomic data found for ",job_name
			    			

		    else:
			print "  ESP fit has finished for ",job_name
			    						    
		
		    
    else:
	print " No json files specified "
        
    
if __name__=="__main__":
    main() 



#
#
#
#
#
#
#
#
#
#def submit_job( struct_dir, pbs_id ,options):
#    import sys, os
#
#    print " submitting job " , pbs_id
#    #submit = "sbatch " + pbs_id
#    submit = options.submit_command +" "+ pbs_id
#
#    #print " sumitting ",submit
#    #sys.exit(' submit_job ')
#    
#    os.system(submit)
# 
# 
#def main():
#    import string, os , sys 
#    # atomicpy
#    import gaussian, elements, xmol , file_io , cluster 
#    from string import replace
#    
#    options, args = get_options()
#    
#    # Set options
#    pmem = options.pmem
#
#    # Read in template files 
#    f = open(options.com_temp,'r')
#    fix_templ = f.read()
#    f.close()
#
#    # 
#    if( options.cluster_host == "peregrine" ):
#	f = open("peregrine.pbs.template",'r') # redmesa.slurm.template
#    elif( options.cluster_host == "redmesa" ):
#	f = open("redmesa.slurm.template",'r')
#    else:
#        # HACK !!!
#        f = open("peregrine.pbs.template",'r') # redmesa.slurm.template
#	
#    pbs_templ = f.read()
#    f.close()
#
#    # Read in index file produced by opv_generator or writen by hand
#    #   Format:
#    #     # - comment
#    #     gen ; main dir ; mol dir ; accuracy ; n units
#    #  Example:
#    #     entry:
#    #       gen  mol   D1_R2R200_A2_R3_  acc1   1
#    #    will read log file:
#    #      mol/D1_R2R200_A2_R3_/acc1_D1_R2R200_A2_R3__n1.fchk
#    #    and make new files in
#    #      mol/D1_R2R200_A2_R3_/
#    
#    
#    # Store working dir  
#    work_dir = os.getcwd()
#    
#    # Read index files from args
#    for indx_file in args:
#        # Get lines of index file   
#        f = open(indx_file,'r')
#        Lines = f.readlines()
#        f.close()
#        for line in Lines:
#            col = line.split()
#            if( len(col) >= 4 and col[0] == "gen" ):
#		
#		
#                mol_dir = col[1].strip()
#                mol_id = col[2].strip()
#                mol_repeat = int(col[3].strip() )
#                mol_acc = col[4].strip()
#                
#                # File info
#                struct_dir = mol_dir + "/" + mol_id + "/"
#                job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)
#                fchk_file = struct_dir + job_name + "/" + job_name + ".fchk"
#		
#                if( file_io.file_exists( fchk_file) ):
#		    if( options.verbose ):
#			
#			print "    Optimization finished  ",fchk_file
#			
#		    NA, ELN, R, TOTAL_ENERGY, Q_ESP    = gaussian.parse_fchk( fchk_file )
#		    ASYMB = elements.eln_asymb(ELN)
#		    
#
#		    # Run esp fit
#		    calc_id_esp = job_name + "-ESP"
#		    input_file = calc_id_esp + ".com"
#		    log_file_esp = struct_dir +'/' + calc_id_esp +"/"+calc_id_esp+".log"
#		    fchk_file_esp = struct_dir +'/' + calc_id_esp +"/"+calc_id_esp+".fchk"
#		    
#		    # Run esp fit if not finished
#		    if( not file_io.file_exists( fchk_file_esp) ):
#                        # Make esp input file
#			if( options.verbose ):
#			    print  '    Writing esp fit input file ',calc_id_esp," in ",struct_dir
#			
#					    
#			os.chdir(struct_dir)
#			gaussian.write_esp_com(calc_id_esp,ASYMB,R)
#    
#                        if( options.submit ):
#                            if( options.verbose ):
#                                print "     Submitting ESP optimization to queue "
#                            # Print pbs script
#                            pbs_id = calc_id_esp+'.pbs'
#                            pbs_name = pbs_id
#                            pbs_dih = pbs_templ
#                            pbs_dih = replace(pbs_dih,"<calc_id>",calc_id_esp)
#                            pbs_dih = replace(pbs_dih,"<input_file>",input_file)
#                            pbs_dih = replace(pbs_dih,"<nnodes>",str( options.nnodes) )
#                            pbs_dih = replace(pbs_dih,"<pmem>",str( options.pmem) )
#                            pbs_dih = replace(pbs_dih,"<npros>",str( options.npros) )
#                            f = file(pbs_name, "w")
#                            f.write(pbs_dih)
#                            f.close()
#                            # Submit job
#                            submit_job( struct_dir, pbs_id ,options )
#                            
#                        elif( options.localrun ):
#                            if( options.verbose ):
#				print  '    Running ESP fit '
#                            gaussian.run(options, calc_id_esp)
#			    
#                        else:
#                            print " Please mark either qsub or runloc options as True to run qm"	
#			
#			os.chdir(work_dir)
#			
#		
#    
#if __name__=="__main__":
#    main() 
#
