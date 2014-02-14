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
    default_basis = '6-31G'
    
    options, args = get_options()
	
	
    json_files = options.json.split(',')
    print json_files
    if( len(json_files) > 0 ):
	# Read index files from args
	for json_file in json_files:
	    # Get lines of index file
		
	    
	    #
	    # Split json file to get struct_dir and job_name
	    #    This assumes 
	    #
	    json_dir_col = json_file.split('/')
	    struct_dir = ""
	    for dir_indx in range( len(json_dir_col) - 1 ):
		dir = json_dir_col[dir_indx]
		struct_dir = struct_dir + dir +'/'
		
	    json_file_col  = json_dir_col[len(json_dir_col)-1].split('.')
	    job_name = json_file_col[0]
	    
	    
	    
	    # Verbose output
	    if( options.verbose ):
		print "The molecules specified in json file ",options.json," will be read in "
    
	
    
	    json_atomicdata = 0
	    json_data,json_success = jsonapy.read_jsondata(json_file)
	    if(  json_success ):
		
		tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
		

		
		if( options.verbose ):
		    print "     Getting atomic data from  ",json_file
    
		ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,success  = jsonapy.read_atomic(json_data)
		if( success ):
		    json_atomicdata = 1
		else:
		    
		    print "   json file ",json_file," exist, but does not contain any atomic data . "
		    
		     
		
	    else:
		print "   json file ",json_file," does not exist. "
		
	    #
	    # Asign method and basis 
	    #
	    if( metadata_found ):
		if( options.verbose ):
		    print " Meta data found will use specified method and basis unless others are specified in the options "

		if( len( options.qm_method) == 0 ):
		    options.qm_method = method		    
		if( len( options.qm_basis) == 0 ):
		    options.qm_basis = basis
	    else:
		# If no meta data is found set accordingly 
		if( len( options.qm_method) == 0 ):
		    options.qm_method = default_method		    
		if( len( options.qm_basis) == 0 ):
		    options.qm_basis = default_basis

	    if( not  json_atomicdata ):
		#
		# If no atomic data try getting it from fchk file 
		#
	    
		print " checking fchk files in ",struct_dir, " job name ",job_name
		
		fchk_file = struct_dir + job_name + '/' + job_name + ".fchk"
		
		fchk_atomicdata = 0
		if( file_io.file_exists(fchk_file) ):
		    if( options.verbose ):
			print " Reading atomic data from ",fchk_file
		    NA, ELN, R, TOTAL_ENERGY , Q_ESP = gaussian.parse_fchk(fchk_file)
		    fchk_atomicdata = 1
		    
		    
	    if( json_atomicdata or fchk_atomicdata ):
		
		# Optimize z-matrix to get bonding information
		qm_kywd_o = options.qm_kywd 
		options.qm_kywd = qm_kywd_o + " OPT"
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
	
		# Print com
		calc_id_temp = 
		gaussian.print_com( calc_id_temp, ASYMB,R,ATYPE,CHARGES,ELECTRONS_i,qm_charge,qm_mult)
		#  geometry z-matrix opt input files
		gaussian.com2zmat(calc_id_temp,calc_id,options)
		
		# Run optimization
		if( options.submit ):
		    if( options.verbose ):
			print "     Submitting Z-matrix optimization to queue "
			print "         nodes        ",options.nnodes
			print "         memory/node  ",options.pmem
		    # Print pbs script
		    pbs_id = cluster.write_pbs(pbs_templ,calc_id,input_file,options)
		    cluster.submit_job( struct_dir, pbs_id ,options )
		    
		elif( options.localrun ):
		    if( options.verbose ):
			print "       Running Z-matrix optimization  "
		    gaussian.run(options, calc_id)
		else:
		    print " Please mark either qsub or runloc options as True to run qm"
					
		options.qm_kywd = qm_kywd_o
		
		    
		
		    
    else:
	print " No json files specified "
        
    
if __name__=="__main__":
    main() 

