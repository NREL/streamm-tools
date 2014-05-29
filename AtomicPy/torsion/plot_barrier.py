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

def d_en_2p(en_1,en_3,h):
    # en_1  eV
    # en_2  eV 
    
    EVKC=23.0605
    
    d_en = (en_3 - en_1)/float(h)/2.0  # eV/deg
    # change to meV
    #r_delta_en = int( 1000*d_en )      # meV/deg
    r_delta_en = 1000*d_en 	    # meV/deg
    
    d_en_p = 0
    if( r_delta_en > 0 ): d_en_p = 1
    if( r_delta_en < 0 ): d_en_p =  -1
	
    return d_en_p

def d2_en_3p(en_1,en_2,en_3,h):
    # en_1  eV
    # en_2  eV 
    # en_3  eV 
    
    d2_en = ( en_3 - 2.0*en_2  + en_1 )/ float(h)**2
    
    return d2_en

def qm_analysis_2(  dih_qm ):
    import numpy, sys  
    
    debug = 1
    
    EVKC=23.0605
    
    success = 0 
    
    qm_min =  1e16 
    qm_max = -1e16
    ang_min = "nan"
    ang_max = "nan"
    
    
    f_qm = open(dih_qm ,"r")
    f_qm_Lines = f_qm.readlines()
    f_qm.close()
    
    # Put energies in arrays 
    tor_en = [] #numpy.zeros()
    tor_angle = [] #numpy.zeros()
    
    for f_qm_line in f_qm_Lines:
	f_qm_col = f_qm_line.split()
    
	if( len(f_qm_col) >= 3 and f_qm_col[0] != "#" ):
	    tor_angle.append( float( f_qm_col[1] ) )
	    tor_en.append( float( f_qm_col[2] ) )
	    
	    success = 1
	    
    qm_min = min( tor_en)
    
    print "  Mim ",qm_min
    
    # Find inversion points in energy
    min_indx = []
    max_indx = []
    trans_list = []
    
    indx_m = len(tor_en) - 1
    en_p = tor_en[ indx_m ]

        
    if(debug):    print "  Initializing last energy as ",len(tor_en),en_p
    
    
    # if( en_p > tor_en[indx_m - 2 ]    ): d_en_p = 1
    # if( en_p < tor_en[indx_m - 2 ]   ): d_en_p =  -1
    
    d_en_p =  d_en(tor_en[indx_m - 2 ] ,en_p) 
	
    if(debug):  print "   en-1",tor_en[indx_m - 2 ] - qm_min ," en ",en_p- qm_min  ," angle ",tor_angle[ indx_m ],"   dE ", d_en_p
    
    # Test first point
    calc_minmaxtrans = 0 
    calc_maxmintrans = 0 
    indx = 0
    en_i = tor_en[indx]
    if( en_i  > en_p ): d_en_i =   1
    if( en_i  < en_p ): d_en_i =  -1
    
    d_en_i  = d_en(en_p ,en_i )
    
    if( d_en_i != 0 ):
	    
	if( d_en_i > d_en_p ):
	    min_indx.append(indx_m)
	    min_en_i = en_i
	    calc_minmaxtrans = 1 
	if( d_en_i < d_en_p ):
	    max_indx.append(indx_m)
	    max_en_i = en_i
	    calc_maxmintrans = 1 
	
	if(debug):  print "   en-1",en_p - qm_min," en ",en_i - qm_min," angle ",tor_angle[ indx ],"   dE ", d_en_i
    
	en_p = en_i
	d_en_p = d_en_i	
	
    # loop over subsequent points 
    for indx in range(1,len(tor_en)):
	en_i = tor_en[indx]
	if( en_i  > en_p ): d_en_i = 1
	if( en_i  < en_p ): d_en_i =  -1
	
	d_en_i  = d_en(en_p ,en_i )
	    
	if( d_en_i != 0 ):
		    
	    if( d_en_i > d_en_p ):
		min_indx.append(indx - 1 )
		min_en_i = en_i
		calc_minmaxtrans = 1 		
		# Calculate transition 
		if( calc_maxmintrans ):
		    trans_ev = max_en_i - min_en_i 
		    trans_list.append(  trans_ev )
		    
		    print "   Found max min trans ",max_en_i- qm_min," -> ",min_en_i - qm_min
		    
		    calc_maxmintrans = 0 
		    
	    if( d_en_i < d_en_p ):
		max_indx.append(indx - 1 )
		max_en_i = en_i
		calc_maxmintrans = 1
		    
		# Calculate transition 
		if( calc_minmaxtrans ):
		    trans_ev = max_en_i - min_en_i
		    trans_list.append(  trans_ev )
		    calc_minmaxtrans = 0 
		    
		    print "   Found min max trans ",max_en_i - qm_min," -> ",min_en_i - qm_min
		    
		    
	
	    
	    if(debug):  print "   en-1",en_p - qm_min," en ",en_i - qm_min," angle ",tor_angle[ indx ],"   dE ", d_en_i
    
	    en_p = en_i
	    d_en_p = d_en_i
    
    # Step over initial point in case of transition
    indx = 1
    en_i = tor_en[indx]
    if( en_i  > en_p ): d_en_i = 1
    if( en_i  < en_p ): d_en_i =  -1
    
    d_en_i  = d_en(en_p ,en_i )
    
    if(debug):  print "   en-1",en_p - qm_min," en ",en_i - qm_min," angle ",tor_angle[ indx ],"   dE ", d_en_i
    
    if( d_en_i != 0 ):
			
	if( d_en_i > d_en_p ):
	    min_en_i = en_i
	    # Calculate transition 
	    if( calc_maxmintrans ):
		trans_ev = max_en_i - min_en_i
		trans_list.append(  trans_ev )
		
		print "   Found max min trans ",max_en_i - qm_min," -> ",min_en_i - qm_min
		
		calc_maxmintrans = 0 
		
	if( d_en_i < d_en_p ):
	    max_en_i = en_i
	    # Calculate transition 
	    if( calc_minmaxtrans ):
		trans_ev = max_en_i - min_en_i
		trans_list.append(  trans_ev )
		calc_minmaxtrans = 0 
		
		print "   Found min max trans ",max_en_i - qm_min," -> ",min_en_i - qm_min
		
	    
    for inv_indx in range( len(min_indx) ):
	indx = min_indx[inv_indx]
	if(debug): print "  Min ",inv_indx," found at ",tor_angle[indx]," w energy ",tor_en[indx]
        
    for inv_indx in range( len(max_indx) ):
	indx = max_indx[inv_indx]
	if(debug): print "  Max ",inv_indx," found at ",tor_angle[indx]," w energy ",tor_en[indx]
    
    for inv_indx in range( len(trans_list) ):
	bar = trans_list[inv_indx]
	print "  Transion ",inv_indx," =  ",bar*EVKC," kcal/mol"
        
	
    if(debug):    sys.exit(" inversion testing ")
	    
    return (success,  tor_angle,tor_en,min_indx,max_indx,trans_list )


def qm_analysis_d2(  dih_qm ):
    import file_io
    import numpy, sys  
    # Find max and mins based on the second derivative 
    
    debug = 0
    
    success = 0 
    
    EVKC=23.0605
    
    success = 0 
    
    qm_min =  1e16 
    qm_max = -1e16
    ang_min = "nan"
    ang_max = "nan"
        

    if( file_io.file_exists(dih_qm) ):
		
		    
	print " reading ",dih_qm
	
	f_qm = open(dih_qm ,"r")
	f_qm_Lines = f_qm.readlines()
	f_qm.close()
	
	# Put energies in arrays 
	tor_en = [] #numpy.zeros()
	tor_angle = [] #numpy.zeros()
	
	for f_qm_line in f_qm_Lines:
	    f_qm_col = f_qm_line.split()
	
	    if( len(f_qm_col) >= 3 and f_qm_col[0] != "#" ):
		tor_angle.append( float( f_qm_col[1] ) )
		tor_en.append( float( f_qm_col[2] ) )
		
		success = 1
		
	qm_min = min( tor_en)
	h=5.0
	
	if( debug): print "  Mim ",qm_min
	
	# Find inversion points in energy
	min_indx = []
	max_indx = []
	trans_list = []
	
	# Test first point
	calc_minmaxtrans = 0 
	calc_maxmintrans = 0
	print " length ",len(tor_en)  
	# loop over subsequent points 
	for indx_i in range(len(tor_en) ):
	    
	    
	    indx_m = indx_i - 1  # _m minus 
	    indx_m_m = indx_m - 1 
	    indx_p = indx_i + 1  # _p plus 
	    indx_p_p = indx_p + 1 
	    
	    # apply boundry conditions 
	    if( indx_m < 0 ): indx_m  =  -1*indx_m 
	    if( indx_m_m < 0 ): indx_m_m  = -1*indx_m_m 
		    
	    if(debug): print "  m i p ",indx_m_m,indx_m,indx_i,indx_p,indx_p_p
    
	    if( indx_p > len(tor_en)  -1  ): indx_p  = indx_p - ( indx_p -  len(tor_en) + 1)*2  # as 0.0 == 180.0 
	    if( indx_p_p > len(tor_en)  -1  ): indx_p_p  = indx_p_p -  ( indx_p_p -  len(tor_en) + 1)*2  #  as 0.0 == 180.0 
	    
	    d_en_m =d_en_2p(tor_en[indx_m] ,tor_en[indx_i],h) 
	    d_en_m_m =d_en_2p(tor_en[indx_m_m] ,tor_en[indx_m],h) 
	    d_en_p =d_en_2p(tor_en[indx_i] ,tor_en[indx_p],h)
	    d_en_p_p =d_en_2p(tor_en[indx_p] ,tor_en[indx_p_p],h)
	    
	    if(debug):
		print "  m i p ",indx_m_m,indx_m,indx_i,indx_p,indx_p_p
		print "  m i p ",tor_en[indx_m_m] - qm_min,tor_en[indx_m] - qm_min,tor_en[indx_i] - qm_min,tor_en[indx_p] - qm_min,tor_en[indx_p_p] - qm_min
		print "     dm dp ",d_en_m_m,d_en_m,d_en_p,d_en_p_p
		#print "     dm dp ",d_en_m,d_en_p
	    
	    
	    if( d_en_m_m < 0 and d_en_m < 0  and d_en_p > 0  and d_en_p_p > 0 ):
	    #if(  d_en_m < 0  and d_en_p > 0  ):
		min_indx.append( indx_i )
		min_en_i = tor_en[indx_i]
		calc_minmaxtrans = 1
		
		if(debug): print " min found ",tor_en[indx_i] - qm_min
		
		if( calc_maxmintrans ):
		    trans_ev = max_en_i - min_en_i
		    trans_list.append(  trans_ev )
		    calc_maxmintrans = 0 
			    
		    if(debug): print "   Found  max min trans ",max_en_i - qm_min," -> ",min_en_i - qm_min," trans = ",trans_ev
			    
		
		
	    if( d_en_m_m > 0 and d_en_m > 0  and d_en_p < 0  and d_en_p_p < 0 ):
	    #if( d_en_m > 0  and d_en_p < 0   ):
		max_indx.append( indx_i )
		max_en_i = tor_en[indx_i]
		calc_maxmintrans = 1
    
		if(debug): print " max found ",tor_en[indx_i] - qm_min
		
		if( calc_minmaxtrans ):
		    trans_ev = max_en_i - min_en_i
		    trans_list.append(  trans_ev )
		    calc_minmaxtrans = 0 
			    
			    
			    
		    if(debug): print "   Found min max trans ",min_en_i - qm_min," -> ",max_en_i - qm_min," trans = ",trans_ev
			    
    
		
	for inv_indx in range( len(min_indx) ):
	    indx = min_indx[inv_indx]
	    if(debug): print "  Min ",inv_indx," found at ",tor_angle[indx]," w energy ",tor_en[indx]
	    
	for inv_indx in range( len(max_indx) ):
	    indx = max_indx[inv_indx]
	    if(debug): print "  Max ",inv_indx," found at ",tor_angle[indx]," w energy ",tor_en[indx]
	
	for inv_indx in range( len(trans_list) ):
	    bar = trans_list[inv_indx]
	    #ang = tor_angle[inv_indx]
	    print "  Transion ",inv_indx," =  ",bar*EVKC," kcal/mol",bar," eV " #at ",ang," degrees "
	    
	    
	if(debug):    sys.exit(" inversion testing ")
    else:
	print " file not found ",dih_qm
	print " return success = 0"
	tor_angle = 0
	tor_en = 0
	min_indx = 0
	max_indx = 0
	trans_list = 0
	 
    return (success,  tor_angle,tor_en,min_indx,max_indx,trans_list )


	    
def main():
    import sys, os , numpy 
    import string  
    import file_io, jsonapy
    from string import replace
    
    options, args = get_options()
    
    lwidth =  5
    n_max  = 14
    
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
			    
			    # Set dihedral count
			    bar_cnt = [-1]*n_max*2 # numpy.zeros( [n_max*2])
    
			    minmax_data = struct_dir +'/' + job_name + '_bar.dat'
			    
			    b_dat = open(minmax_data ,"w")			
			    b_dat.write("# n_units ; d_indx , qm_bar (eV) " )
			    b_dat.close()

			    minmax_data2 = struct_dir +'/' + job_name + '_bar2.dat'

			    b2_dat = open(minmax_data2 ,"w")			
			    b2_dat.write("\n # n_units ; d_pos ; min angle/val N; max angle/val N ; barriers N" ) #% ( len(min_indx),len( max_indx ) , len( trans_list) ) )
			    #b2_dat.write("\n # n_units ; d_pos ; min angle/val %d ; max angle/val %d ; barriers %d " ) #% ( len(min_indx),len( max_indx ) , len( trans_list) ) )
			    b2_dat.close()

				
			    bplot_l =  " \'"+minmax_data+"\' " + ' us 2:3  '+' with linespoints ls '+str(int(n_units)+1)+'  title '+ "\' " + tag+"n="+str(n_units) +" \'  , \\" + "\n"
			    
			    print bplot_l 
			    
			    bar_plot_lines.append( bplot_l )
				
			    calc_i += 1 
				    
			    style_l =  'set style line ' + str(calc_i+1) + ' lt  ' +  str(calc_i) + ' lw '+str(lwidth) + ' lc ' +  str(calc_i) + "\n"
			    style_lines.append(  style_l )
			    
			    # center barriers 
			    minmax_data_cent = struct_dir +'/' + tag + '_bar.dat'
			    
			    centbar_data = struct_dir +'/' + tag + options.qm_sufix + '_centerbar.dat'
			    if ( not file_io.file_exists(centbar_data) ):
				
				b_dat = open(centbar_data ,"w")			
				b_dat.write("# n_units ; d_indx , qm_bar " )
				b_dat.close()
				
				centbar_line =  " \'"+centbar_data+"\' " + ' us 2:3  '+' with linespoints ls '+str(int(n_units)+1)+'  title '+ "\' " + tag+"n="+str(n_units) +" \'  , \\" + "\n"
				centbar_plot_lines.append( centbar_line )
				
				print centbar_line
					    
			    for dih_indx in range( len(qm_dih_id_list) ):
				dih_id = qm_dih_id_list[dih_indx]
				cent_min = qm_cent_min_list[dih_indx]
				cent_max = qm_cent_max_list[dih_indx]
				cent_step = qm_cent_step_list[dih_indx]
				
				dih_qm = struct_dir +'/' +job_name + '-' + dih_id + options.qm_sufix + '.dat'
				
				# success, qm_min, ang_min, qm_max, ang_max = qm_analysis_2(  dih_qm )
				
				success,  tor_angle,tor_en,min_indx,max_indx,trans_list = qm_analysis_d2(  dih_qm )
				
				bar_cnt[n_units] += 1 
		    
				if( success ):		    
				    d_o = 1*(n_units-2)
				    
				    dih_cnt =  int( bar_cnt[n_units] )
				    d_indx = n_max + d_o + dih_cnt*2
				    d_pos =   d_o - dih_cnt*2 
				    
				    min_angles = ""
				    for indx in range( len(min_indx) ):
					dih_indx = min_indx[indx]
					angle = tor_angle[dih_indx]
					min_angles = min_angles +"  " +str( angle )
				    
				    min_enval = ""
				    for indx in range( len(min_indx) ):
					dih_indx = min_indx[indx]
					en_val = tor_en[dih_indx]
					min_enval = min_enval +"  " +str( en_val )
					
				    
				    max_angles = ""
				    for indx in range( len(max_indx) ):
					dih_indx = max_indx[indx]
					angle = tor_angle[dih_indx]
					max_angles = max_angles +"  " +str( angle )
				    
				    max_enval = ""
				    for indx in range( len(max_indx) ):
					dih_indx = max_indx[indx]
					en_val = tor_en[dih_indx]
					max_enval = max_enval +"  " +str( en_val )
					
				    qm_bar = ""
				    for indx in range( len(trans_list) ):
					qm_bar = qm_bar +"  " + str( trans_list[indx] )


				    glb_min_en = 1e16
				    glb_min_angle = 1e16
				    glb_max_en = -1e16
				    glb_max_angle = -1e16
				    
				    
				    
#    for inv_indx in range( len(trans_list) ):
#	bar = trans_list[inv_indx]
#	ang = tor_angle[inv_indx]
#	print "  Transion ",inv_indx," =  ",bar*EVKC," kcal/mol",bar," eV at ",ang," degrees "
        
	
				    for indx in range( len(tor_en) ):
					if(  tor_en[indx] < glb_min_en ):
					    glb_min_en =  tor_en[indx]
					    glb_min_angle = tor_angle[dih_indx]
					

				    for indx in range( len(tor_en) ):					
					if(  tor_en[indx] > glb_max_en ):
					    glb_max_en =  tor_en[indx]
					    glb_max_angle = tor_angle[dih_indx]
				    
                                    delat_maxmin = glb_max_en - glb_min_en
				    

				    if( options.verbose ):
					print "   units ",n_units, " dih id ",dih_id,"dih pos",d_pos," dih cnt ",dih_cnt," do ", d_o, d_pos ,d_indx, qm_bar
				    
				    b2_dat = open(minmax_data2 ,"a")			
				    b2_dat.write("\n  %s %d %d %s %s %s  %s %s " % (dih_id, n_units , d_pos, min_angles,min_enval, max_angles, max_enval , qm_bar ) )
				    b2_dat.close()

				    b_dat = open(minmax_data ,"a")			
				    b_dat.write("\n  %s  %d %d %f %f %f %f %f " % ( dih_id,n_units , d_pos, glb_min_angle,glb_min_en, glb_max_angle, glb_max_en, delat_maxmin ) )
				    b_dat.close()
	    				    
				    # Write
				    if( d_pos  == 0 or d_pos  == 1  ):
					centbar_data = struct_dir +'/' + tag + options.qm_sufix + '_centerbar.dat'
					b_dat = open(centbar_data ,"a")	
					b_dat.write("\n  %s %d %d %s" % ( dih_id, n_units , d_pos, delat_maxmin ) )
					b_dat.close()
					
					print " Center relative barriers ",minmax_data_cent
					
					b_datcent = open(minmax_data_cent ,"a")			
					b_datcent.write("\n  %s %d %d %s %s %s  %s %s " % (dih_id, n_units , d_pos, min_angles,min_enval, max_angles, max_enval , qm_bar ) )
					b_datcent.close()
	    
    if( options.verbose ):
	print "  Plotting data from " ,indx_file

    style_i = ""
    for s_line in style_lines:
	style_i += s_line
	
    plot_i = ""
    for p_line in bar_plot_lines :
	plot_i += p_line
	
    plot_cent = ""
    for p_line in centbar_plot_lines :
	plot_cent += p_line

    plt_template = 'dihbar_plt.template'
    f = open( plt_template , 'r')
    plt_file = f.read()
    f.close()
    
    plt_posbar = replace(plt_file,'<style_lines>',style_i)
    plt_posbar = replace(plt_posbar,'<plot_lines>',plot_i)
    #plt_posbar = replace(plt_posbar,'<plot_lines2>',plot_i2)
    plt_posbar = replace(plt_posbar,'<structure_name>',tag)
	
	
    plt_dir_file =  replace(indx_file,'.rec','_bar.plt')
    
    f = open( plt_dir_file , 'w')
    f.write(plt_posbar)
    f.close()
    
    # centbar_plot_lines.append( centbar_line )
    
    

    print " plot_cent ",plot_cent

    plt_posbar = replace(plt_file,'<style_lines>',style_i)
    plt_posbar = replace(plt_posbar,'<plot_lines>',plot_cent)
    plt_posbar = replace(plt_posbar,'<structure_name>',tag)
	
	
    plt_dir_file =  replace(indx_file,'.rec','_centbar.plt')
    
    f = open( plt_dir_file , 'w')
    f.write(plt_posbar)
    f.close()
    

if __name__=="__main__":
    main()

