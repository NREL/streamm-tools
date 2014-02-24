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

def qm_analysis(  dih_qm ):
    
    success = 0 
    
    qm_min =  1e16 
    qm_max = -1e16
    ang_min = "nan"
    ang_max = "nan"
    
    #ang_min 
    
    f_qm = open(dih_qm ,"r")
    f_qm_Lines = f_qm.readlines()
    f_qm.close()
    
    
    for f_qm_line in f_qm_Lines:
	f_qm_col = f_qm_line.split()
    
	if( len(f_qm_col) >= 3 and f_qm_col[0] != "#" ):
	    qm_angle = float( f_qm_col[1] )
	    qm_en = float( f_qm_col[2] )
	    if( qm_en < qm_min ):
		qm_min = qm_en
		ang_min = qm_angle
	    if( qm_en > qm_max ):
		qm_max = qm_en
		ang_max = qm_angle
		
	    success = 1 
	    
    return (success, qm_min,ang_min,qm_max,ang_max )


	    
def main():
    import sys, os , numpy 
    import string  
    import file_io, jsonapy
    from string import replace
    
    options, args = get_options()
    
    lwidth =  10
    n_max  = 14
    bar_cnt = numpy.zeros( [n_max*2])
    
    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
	
	
	# Initialize lists
	style_lines = []
	bar_plot_lines = []
	centbar_plot_lines = []
	calc_i = 0 
	
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
    
			    minmax_data = struct_dir +'/' + job_name + '_bar.dat'
			    
			    b_dat = open(minmax_data ,"w")			
			    b_dat.write("# n_units ; d_indx , qm_bar (eV) " )
			    b_dat.close()
				
			    bplot_l =  " \'"+minmax_data+"\' " + ' us 2:3  '+' with linespoints ls '+str(int(n_units)+1)+'  title '+ "\' " + tag+"n="+str(n_units) +" \'  , \\" + "\n"
			    bar_plot_lines.append( bplot_l )
				
			    centbar_data = struct_dir +'/' + tag + options.qm_sufix + '_centerbar.dat'
			    if ( not file_io.file_exists(centbar_data) ):
				b_dat = open(centbar_data ,"w")			
				b_dat.write("# n_units ; d_indx , qm_bar " )
				b_dat.close()
				
				centbar_line =  " \'"+centbar_data+"\' " + ' us 2:3  '+' with linespoints ls '+str(int(n_units)+1)+'  title '+ "\' " + tag+"n="+str(n_units) +" \'  , \\" + "\n"
				centbar_plot_lines.append( bplot_l )
					    
			    for dih_indx in range( len(qm_dih_id_list) ):
				dih_id = qm_dih_id_list[dih_indx]
				cent_min = qm_cent_min_list[dih_indx]
				cent_max = qm_cent_max_list[dih_indx]
				cent_step = qm_cent_step_list[dih_indx]
				
				dih_qm = struct_dir +'/' +job_name + '-' + dih_id + options.qm_sufix + '.dat'
				
				success, qm_min, ang_min, qm_max, ang_max = qm_analysis(  dih_qm )
		    
		    
				if( success ):
				    calc_i += 1 
				    qm_bar = qm_max - qm_min
		    
				    d_o = -1*(n_units-2)
				    
				    dih_cnt =  int( bar_cnt[n_units] )
				    d_indx = n_max + d_o + dih_cnt*2
				    d_pos = d_indx - n_max
				    if( options.verbose ):
					print "   grid point ",n_units, dih_cnt, d_o, d_pos ,d_indx, qm_bar
				    
				    
				    style_l =  'set style line ' + str(calc_i+1) + ' lt  ' +  str(calc_i) + ' lw '+str(lwidth) + ' lc ' +  str(calc_i) + "\n"
				    style_lines.append(  style_l )
		    
				    b_dat = open(minmax_data ,"a")			
				    b_dat.write("\n %d %d %f" % ( n_units , d_pos, qm_bar ) )
				    b_dat.close()
				    
				    # Write
				    if( d_pos  == 0 or d_pos  == 1  ):
					centbar_data = struct_dir +'/' + tag + options.qm_sufix + '_centerbar.dat'
					b_dat = open(centbar_data ,"a")	
					b_dat.write("\n %d %d %f" % ( n_units , d_pos, qm_bar ) )
					b_dat.close()
					
				    bar_cnt[n_units] += 1 
		
        if( options.verbose ):
            print "  Plotting data from " ,indx_file

        style_i = ""
        for s_line in style_lines:
            style_i += s_line
            
        plot_i = ""
        for p_line in bplot_l :
            plot_i += p_line
	    
        plot_i2 = ""
        for p_line in centbar_plot_lines :
            plot_i2 += p_line
    
	plt_template = 'dihbar_plt.template'
	f = open( plt_template , 'r')
	plt_file = f.read()
	f.close()
	
        plt_file = replace(plt_file,'<style_lines>',style_i)
        plt_file = replace(plt_file,'<plot_lines>',plot_i)
        plt_file = replace(plt_file,'<plot_lines2>',plot_i2)
        plt_file = replace(plt_file,'<structure_name>',tag)
	    
	    
	plt_dir_file =  replace(indx_file,'.rec','_bar.plt')
	
        f = open( plt_dir_file , 'w')
        f.write(plt_file)
        f.close()
	

if __name__=="__main__":
    main()

