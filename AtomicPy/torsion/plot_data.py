# Plot torsional data based on reference file 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    # should be reference file 
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    parser.add_option("--tor_id", dest="tor_id",type="string",default="dih",help=" plot id for postscript output  ")
    
    parser.add_option("--ff_software", dest="ff_software",type="string",default="lammps",help=" FF software ")
    
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


def ff_analysis(  dih_file ):
    
    success = 0 
    
    en_min =  1e16 
    en_max = -1e16
    ang_min = "nan"
    ang_max = "nan"
    
    #ang_min 
    
    f_en = open(dih_file ,"r")
    f_en_Lines = f_en.readlines()
    f_en.close()
    
    
    for f_en_line in f_en_Lines:
	f_en_col = f_en_line.split()
    
	if( len(f_en_col) >= 3 and f_en_col[0] != "#" ):
	    en_angle = float( f_en_col[1] )
	    en_en = float( f_en_col[2] )
	    if( en_en < en_min ):
		en_min = en_en
		ang_min = en_angle
	    if( en_en > en_max ):
		en_max = en_en
		ang_max = en_angle
		
	    success = 1 
	    
    return (success, en_min,ang_min,en_max,ang_max )


def find_barriers():
    
    qm_bar = qm_max - qm_min

    d_o = -1*(n_units-2)
    
    dih_cnt = bar_cnt[mol_i,n_units]
    d_indx = n_max + d_o + dih_cnt*2
    
    bar_grid[n_units,d_indx ] = qm_bar
    
    if( options.verbose ):
	print "   grid point ",n_units, d_o + dih_cnt*2 ,d_indx
    
    minmax_data = struct_dir +'/' + mol_acc + "_" + mol_id +"-"+ str(n_units)+ options.qm_sufix + '_bar.dat'
    if( dih_cnt == 0 ):
	b_dat = open(minmax_data ,"w")			
	b_dat.write("# n_units ; d_indx , qm_bar " )
	b_dat.close()
	
	bplot_l =  " \'"+minmax_data+"\' " + ' us 2:3  '+' with linespoints ls '+str(int(n_units)+1)+'  title '+ "\' " + mol_id+"n="+str(n_units) +" \'  , \\" + "\n"
	bar_plot_lines[mol_i].append( bplot_l )
	
    b_dat = open(minmax_data ,"a")			
    b_dat.write("\n %d %d %f" % ( n_units , d_indx - n_max, qm_bar ) )
    b_dat.close()
    
    # Write
    if( (d_indx - n_max) == 0 ):
	zerobar_data = struct_dir +'/' + mol_acc + "_" + mol_id + options.qm_sufix + '_bar.dat'
	b_dat = open(zerobar_data ,"a")	
	b_dat.write("\n %d %d %f" % ( n_units , d_indx - n_max, qm_bar ) )
	b_dat.close()
	
	#z_plot_line[mol_i] += " \'"+zerobar_data+"\' " + ' us 1:3  with linespoints  '
	
	#b_dat = open(minmax_data ,"w")			
	#b_dat.write("# n_units ; d_indx , qm_bar " )
	#b_dat.close()

    # Write
    if( (d_indx - n_max) == 1 ):
	zerobar_data = struct_dir +'/' + mol_acc + "_" + mol_id + options.qm_sufix + 'odd_bar.dat'
	b_dat = open(zerobar_data ,"a")	
	b_dat.write("\n %d %d %f" % ( n_units , d_indx - n_max, qm_bar ) )
	b_dat.close()
	
	#z_plot_line[mol_i] += " \'"+zerobar_data+"\' " + ' us 1:3  with linespoints  '
			
    bar_cnt[mol_i,n_units] += 1
    
def main():
    import sys, os , numpy 
    import string  
    import file_io, jsonapy
    from string import replace
    
    options, args = get_options()
    
    lwidth =  10
    
    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
	
	
	# Initialize lists
	style_lines = []
	plot_lines = []
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
	    
			    for dih_indx in range( len(qm_dih_id_list) ):
				dih_id = qm_dih_id_list[dih_indx]
				cent_min = qm_cent_min_list[dih_indx]
				cent_max = qm_cent_max_list[dih_indx]
				cent_step = qm_cent_step_list[dih_indx]
				
				dih_qm = struct_dir +'/' +job_name + '-' + dih_id + options.qm_sufix + '.dat'
				
				success, qm_min, ang_min, qm_max, ang_max = qm_analysis(  dih_qm )
		    
		    
				if( success ):
				    
				    print " qm min ",qm_min
				    calc_i += 1 
				    
				    
				    style_l =  'set style line ' + str(calc_i+1) + ' lt  ' +  str(calc_i) + ' lw '+str(lwidth) + ' lc ' +  str(calc_i) + "\n"
				    style_lines.append(  style_l )

				    plot_l =  " \'"+dih_qm+"\' " + ' us 2:($3- ' + str(qm_min) + ')*EVKC  axis x1y1 w l ls '+str(calc_i+1)+'  title '+ "\' " + mol_dir +" n="+str(n_units)+" dihedral "+dih_id+" \' smooth unique  , \\" + "\n"
				    plot_lines.append(  plot_l )
			    
				    plot_l =  " \'"+dih_qm+"\' " + ' us 2:($3- ' + str(qm_min) + ')*1000  axis x1y2 w l ls '+str(calc_i+1)+" notitle  smooth unique  , \\" + "\n"
				    plot_lines.append(  plot_l )
				    
				    # xmol_qm = job_name+'-'+dih_id+"_qm.xmol"
				    
				    
				    #if( options.plot_barrier ):
							
				    
			
	
			ff_dih_id_list ,ff_cent_min_list ,ff_cent_max_list ,ff_cent_step_list,ff_a_k_list, ff_a_i_list, ff_a_j_list, ff_a_l_list,ff_type_list,fftor_found = jsonapy.read_ff_tor(json_data)
    
			if( fftor_found ):
	    
			    for dih_indx in range( len(ff_dih_id_list) ):
				dih_id = ff_dih_id_list[dih_indx]
				cent_min = ff_cent_min_list[dih_indx]
				cent_max = ff_cent_max_list[dih_indx]
				cent_step = ff_cent_step_list[dih_indx]
				ff_type_id = ff_type_list[dih_indx]

				dih_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + options.ff_software + ff_type_id +".dat"
				
				success, ff_min, ang_min, ff_max, ang_max = ff_analysis(  dih_ff )
		    

				if( success ):
				    
				    print " ff min ",ff_min
				    calc_i += 1 
				    
				    style_l =  'set style line ' + str(calc_i+1) + ' lt  ' +  str(calc_i) + ' lw '+str(lwidth) + ' lc ' +  str(calc_i) + "\n"
				    plot_l =  " \'"+dih_ff+"\' " + ' us 2:($3- ' + str(ff_min) + ')*EVKC w l ls '+str(calc_i+1)+'  title '+ "\' " + mol_dir+" n="+str(n_units)+" dihedral "+dih_id+"FF \' smooth unique  , \\" + "\n"
			    
				    style_lines.append(  style_l )
				    plot_lines.append(  plot_l )
				
				
	# Use name of index file as plot names
	plt_id =  replace(indx_file,'rec','')
	
        if( options.verbose ):
            print "  Plotting data from " ,indx_file

        style_i = ""
        for s_line in style_lines:
            style_i += s_line
            
        plot_i = ""
        for p_line in plot_lines :
            plot_i += p_line
    
	plt_template = 'dih_plt.template'
	f = open( plt_template , 'r')
	plt_file = f.read()
	f.close()
	
        plt_file = replace(plt_file,'<style_lines>',style_i)
        plt_file = replace(plt_file,'<plot_lines>',plot_i)
        plt_file = replace(plt_file,'<structure_name>',plt_id)
	    
	    
	plt_dir_file =  plt_id + 'plt'
	
        f = open( plt_dir_file , 'w')
        f.write(plt_file)
        f.close()
	

if __name__=="__main__":
    main()

