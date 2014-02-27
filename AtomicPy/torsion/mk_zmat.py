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
    parser.add_option("--pbs_template", dest="pbs_template",type="string",default="",help=" Template for job submission  ")

    # How to run the needed calculations 
    parser.add_option("--submit", dest="submit",action="store_true", default=False,help=" submit calculations to the queue ")
    parser.add_option("--localrun", dest="localrun",action="store_true", default=False,help=" Run calculations locally")
    parser.add_option("--submit_command", dest="submit_command",type="string", default="qsub",help=" command used to submit script to the queue ")

    parser.add_option("--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")
    
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
    import gaussian, elements, xmol , file_io , cluster 
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
		    calc_id = "%s/%s%s" % (struct_dir, job_name , "-ZMATOPT" )
		    #
		    # 
		    #
		    zmat_fchk = "%s/%s%s" % ( calc_id ,job_name,"-ZMATOPT.fchk" )
		    print " Checking for complete zmatrix optimiztion ",zmat_fchk
                    zmat_finished = file_io.file_exists( zmat_fchk )
		    
		    if( not zmat_finished or options.recalc ):
			
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
		    
			print " checking fchk files in ",struct_dir, " job name ",job_name
			
			fchk_file = struct_dir + job_name + '/' + job_name + ".fchk"
			
			fchk_atomicdata = 0
			if( file_io.file_exists(fchk_file) ):
			    if( options.verbose ):
				print " Reading atomic data from ",fchk_file
			    NA, ELN, R, TOTAL_ENERGY , Q_FCHK = gaussian.parse_fchk(fchk_file)
			    fchk_atomicdata = 1
			    
					
			    # Poppulate other atomic values                
			    ASYMB = elements.eln_asymb(ELN)
			    CHARGES = []
			    for atom_i in range(NA):
				CHARGES.append( -100.0 )
				
			if( json_atomicdata or fchk_atomicdata ):
			    
			    ATYPE = []
			    ELECTRONS_i = 0
			    for atom_i in range( len(ELN) ):
				ATYPE.append( ASYMB[atom_i] )
				ELECTRONS_i += ELN[atom_i] 
			    
		    
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
			    job_id =  "%s%s" % (job_name,"-ZMATOPT" )
		    
			    # Print com
			    calc_id_temp = job_id + "-temp"
			    gaussian.print_com( calc_id_temp, ASYMB,R,ATYPE,CHARGES,ELECTRONS_i,options.qm_method,options.qm_basis,options.qm_kywd,options.qm_charge,options.qm_mult)
			    #  geometry z-matrix opt input files
			    gaussian.com2zmat(calc_id_temp,job_id,options)
			    
			    # Run optimization
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
	print " No json files specified "
        
    
if __name__=="__main__":
    main() 

