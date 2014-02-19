# Get data from remote cluster 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")
    parser.add_option("-d","--json_dir", dest="json_dir", default="/home/${USER}/",type="string",help=" location of json file on host ")    
    parser.add_option("--update_info", dest="update_info",type="string",default="True",help=" download reference file from cluster ")

    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")

    # should be reference file 
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    
    (options, args) = parser.parse_args()
    
    return options, args
    
def get_info(cluster_id,user_id,json_dir,json_file):
    import sys, os
    
    get_dat = "scp "+ user_id+"@"+cluster_id + json_dir +"/"+ json_file +' ./ \n'
    print get_dat
    os.system(get_dat)
    
def main():
    import sys, os 
    import string
    import file_io, jsonapy
    
    options, args = get_options()

    if( options.host == "peregrine" ):
        cluster_id = 'peregrine-login1.nrel.gov:'
        user_id = 'tkemper'
        
    else:
        print " unknow host id "
        sys.exit(" option error ")
    
        
    print " Downloading data files for ",user_id,"@",cluster_id


    
    json_files = options.json.split(',')
    print json_files
    if( len(json_files) > 0 ):
	# Read index files from args
	for json_file in json_files:
	    
                    
            if( options.update_info == "True" ):
                get_info(cluster_id,user_id,options.json_dir,options.json)

	    # Verbose output
	    if( options.verbose ):
		print "The molecules specified in json file ",options.json," will be read in "
    
	    json_data,json_success = jsonapy.read_jsondata(json_file)
	    if(  json_success ):
		
		mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found  ):
		    
		    if( options.verbose ):
			print " Meta data found will use specified method and basis unless others are specified in the options "
		    #
		    # Construct file names 
		    #
		    #short_name = "acc%d_%s_n%d" % (accuracy, tag, number )
		    job_name = "acc%d_%s_n%d" % (accuracy, tag, n_units )
		    struct_dir = "%s/%s/" % (mol_dir, tag )
	
                    if( options.verbose ):
                        print "   Getting ",struct_dir,job_name
            
                    qm_dih_id_list ,qm_cent_min_list ,qm_cent_max_list ,qm_cent_step_list,qm_a_k_list, qm_a_i_list, qm_a_j_list, qm_a_l_list,qmtor_found = jsonapy.read_qm_tor(json_data)
                    ff_dih_id_list ,ff_cent_min_list ,ff_cent_max_list ,ff_cent_step_list,ff_a_k_list, ff_a_i_list, ff_a_j_list, ff_a_l_list,ff_type_list,fftor_found = jsonapy.read_ff_tor(json_data)
                    if( qmtor_found or fftor_found ):
                    
                        mk_dir = 'mkdir -p '+struct_dir+'/'
                        os.system(mk_dir)
                        mv_json = 'mv '+json_file + " "+struct_dir+'/'
                        os.system(mv_json)
                        
                        
                        # Get data files 
                        get_dat = "scp "+ user_id+"@"+cluster_id + options.json_dir +"/"+ job_name +'*.dat ' + struct_dir+' \n'
                        print "  Getting data files ",get_dat
                        os.system(get_dat)

                        # Get xmol files 
                        get_dat = "scp "+ user_id+"@"+cluster_id + options.json_dir +"/"+ job_name +'*.xmol ' + struct_dir+' \n' 
                        print "  Getting xmol files ",get_dat
                        os.system(get_dat)
                        
                        # Append
                        json_rec = job_name + '.rec'
                        json_rec_lines = user_id+" "+cluster_id + " "+ options.json_dir  + " " + struct_dir +" "+ job_name
                        
                        F = open( json_rec , 'a' )
                        F.write("%s \n " % (json_rec_lines) )
                        F.close()
                        


if __name__=="__main__":
    main()

