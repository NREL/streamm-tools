# Plot torsional data based on reference file 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    # should be reference file 
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    
    (options, args) = parser.parse_args()
    
    return options, args


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
            DIH_VAL.append(  float(col[4]) )
            a_l = int( col[5] ) - 1
            a_i = int( col[6] ) - 1
            a_j = int( col[7] ) - 1
            a_k = int( col[8] ) - 1
            DIH_ATOMS.append( [ a_l,a_i,a_j,a_k ] )
        
  
    return (  DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS )
    	    
def main():
    import sys, os , numpy 
    import string  
    import file_io, jsonapy
    from string import replace
    
    options, args = get_options()
    
    lwidth =  5
    n_max  = 14

    # Initialize lists
    style_lines = []
    ang_plot_lines = []
    centang_plot_lines = []
    calc_i = 0 
    
	
    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
	
	
        for line in Lines:
            col = line.split()
            if( len(col) >= 4 and col[0] != "#" ):
		
		
		struct_dir =  col[3]
		job_name = col[4]
		
		json_file =  struct_dir + "/" + job_name +".json"
		json_data,json_success = jsonapy.read_jsondata(json_file)
		if(  json_success ):
		    
		    mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)		
		    #
		    # Need meta data to proceed 
		    #      		    
		    if( metadata_found  ):
	    
			print "   Plotting ",struct_dir,job_name
			# Open output file
				
			qm_dih_id_list ,qm_cent_min_list ,qm_cent_max_list ,qm_cent_step_list,qm_a_k_list, qm_a_i_list, qm_a_j_list, qm_a_l_list,qmtor_found = jsonapy.read_qm_tor(json_data)
			if( qmtor_found ):
			    
			    # Check to see if list of dihedrals to loop over has been creaeted
			    dlist_name = struct_dir +'/' + job_name + "_dih.list"
			    dlist_exists = file_io.file_exists( dlist_name )
			    if ( dlist_exists ):
			        if( options.verbose ):
				    print "       Reading in dihedral list from   ",dlist_name
				DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS = read_dihlist(dlist_name)
				
				
				    
				# Set dihedral count
				bar_cnt = -1 # numpy.zeros( [n_max*2])
	
				optangle_data = struct_dir +'/' + job_name + '_optangle.dat'
				optangle_data_cent = struct_dir +'/' + tag + '_cent_optangle.dat'
				
				if ( not file_io.file_exists(optangle_data) ):
				    
				    ang_dat = open(optangle_data ,"w")			
				    ang_dat.write("# n_units ; d_indx , optimized angle (deg) " )
				    ang_dat.close()
    
				    ang_line =  " \'"+optangle_data+"\' " + ' us 2:4  '+' with linespoints ls '+str(int(n_units)+1)+'  title '+ "\' " + job_name + " \'  , \\" + "\n"
				    ang_plot_lines.append( ang_line )
				    
				if ( not file_io.file_exists(optangle_data_cent) ):
				    
				    centang_dat = open(optangle_data_cent ,"w")			
				    centang_dat.write("# n_units ; d_indx , optimized angle (deg) " )
				    centang_dat.close()
    
				    centang_line =  " \'"+optangle_data_cent+"\' " + ' us 1:4  '+' with linespoints ls '+str(int(n_units)+1)+'  title '+ "\' " + tag + " \'  , \\" + "\n"
				    centang_plot_lines.append( centang_line )
				    
				
				for dih_indx in range(len(DIH_ID) ):
				    if( DIH_TAG[dih_indx] == "loop" or DIH_TAG[dih_indx] == "fix"  ):
					    
					bar_cnt += 1 
					d_o = 1*(n_units-2)
					dih_cnt =  int( bar_cnt )
					d_indx = n_max + d_o + dih_cnt*2
					d_pos =   d_o - dih_cnt*2
					
					if( d_pos >= 0):
						
					    print n_units , d_pos, DIH_VAL[dih_indx]
    
					    ang_dat = open(optangle_data ,"a")			
					    ang_dat.write("\n %d %d %f  %f " % (n_units , d_pos, DIH_VAL[dih_indx], numpy.absolute( DIH_VAL[dih_indx]) ))
					    ang_dat.close()
				    
					    if( d_pos  == 0 or d_pos  == 1  ):
						print " Center angle ",d_pos, DIH_VAL[dih_indx]
							
						centang_dat = open(optangle_data_cent ,"a")			
						centang_dat.write("\n %d %d %f  %f " % (n_units , d_pos, DIH_VAL[dih_indx], numpy.absolute( DIH_VAL[dih_indx]) ))
						centang_dat.close()


			    else:
				print "       Dihedral list   ",dlist_name," does not exist"


			    print " Optimized angles will be writen to ",optangle_data
			    print " Center chain optimized angles will be writen to ",optangle_data_cent
			    
			    print ""
			    data_files = struct_dir +'/' + "*" + tag + "*" +'_optangle.dat'
			    print " rm ",data_files
			    
		else:
		    print  " json file read not sucessful "
		    
    # Print gnuplot lines
    
    for plot_l in ang_plot_lines:
	print plot_l
	
    for plot_l in centang_plot_lines:
	print plot_l
	
    
if __name__=="__main__":
    main()

