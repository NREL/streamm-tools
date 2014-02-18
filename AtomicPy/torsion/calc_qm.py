#! /usr/bin/env python
# run ab initio torsional potential 

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

    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    parser.add_option("--qm_load", dest="qm_load",type="string",help=" string to load qm software module  ")
    
    (options, args) = parser.parse_args()

    # Set options based on cluster 
    if( options.host == "peregrine" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"
    elif( options.host == "redmesa" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/g09/C.01"

    return options, args

def main():
    import sys, os , string
    import jsonapy
    import file_io, gaussian, cluster

    options, args = get_options()
    
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
		    
		    for dih_indx in range( len(dih_id_list) ):
			dih_id = dih_id_list[dih_indx]
			cent_min = cent_min_list[dih_indx]
			cent_max = cent_max_list[dih_indx]
			cent_step = cent_step_list[dih_indx]
			
			print "   Running job ",struct_dir,job_name,dih_id 
			# Loop over angels of central dihedrals
			for cent_angle in range(cent_min,cent_max,cent_step):
			    calc_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_auxfix'
					
			    # Check to see if finished
			    #  log_file = struct_dir +'/' + calc_id +"/"+calc_id+".log"
			    #  run_qm = gaussian.check_log(log_file)
			    fchk_file = struct_dir +'/' + calc_id +"/"+calc_id+".fchk"
			    run_qm = gaussian.check_fchk(fchk_file)
	
			    qm_finished = 0
			    if( run_qm ):
				input_file =  calc_id + '.com'
				com_fdir =  struct_dir +'/' + calc_id + '.com'
		
				if ( file_io.file_exists(com_fdir) ):
				    os.chdir(struct_dir)
				    if( options.submit  ):
					pbs_id = calc_id+'.pbs'
					cluster.submit_job( struct_dir, pbs_id ,options )
					
				    elif(options.localrun):
					gaussian.run(options, calc_id)
				    os.chdir(work_dir)
				else:
				    print ' error in input file ',input_file
				    sys.exit(' no input file ')
				    

if __name__=="__main__":
    main() 

                