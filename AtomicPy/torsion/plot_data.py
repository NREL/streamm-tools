# Plot torsional data based on reference file 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=True, help="Verbose output ")

    # should be reference file 
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    
    (options, args) = parser.parse_args()
    
    return options, args

def qm_analysis(  dih_qm ):
    
    qm_min =  1e16 
    qm_max = -1e16
    
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
	    
    return ( qm_min,ang_min,qm_max,ang_max )
	    
def main():
    import sys, os , numpy 
    import string  
    import file_io
    from string import replace
    
    options, args = get_options()


    mol_list = []
    mdir_list = []
    mid_list = []
    n_list = []
    acc_list = []
    
    style_lines = []
    plot_lines = []
    bar_plot_lines = []
    qm_min_list  = []
    z_plot_line = []
    
    n_max  = 14
    bar_cnt = numpy.zeros( [n_max,n_max*2])
    bar_grid = numpy.zeros( [n_max,n_max*2])



    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
        for line in Lines:
            col = line.split()
            if( len(col) >= 4 and col[0] == "gen" ):
		
			    
		mol_dir = col[1].strip()
		mol_id = col[2].strip()
		mol_repeat = int(col[3])
		mol_acc = col[4].strip()
		
		# File info
		struct_dir = mol_dir + "/" + mol_id + "/"
		job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)
		
		
                new_mol = 1
                for m_indx in range(len(mol_list)):
                    if( mol_id == mol_list[m_indx] ):
                        new_mol = 0
                        mdir_list[m_indx].append( mol_dir )
                        mid_list[m_indx].append( mol_id )
                        n_list[m_indx].append( mol_repeat )
                        acc_list[m_indx].append( mol_acc )
                        
                if( new_mol ):
                    mol_list.append( mol_id)
                    m_indx = len(mol_list)  - 1
                    mdir_list.append( [] )
                    mdir_list[m_indx].append( mol_dir )
                    mid_list.append( [] )
                    mid_list[m_indx].append( mol_id )
                    n_list.append( [] )
                    n_list[m_indx].append( mol_repeat )
                    acc_list.append( [] )
                    acc_list[m_indx].append( mol_acc )
                    
                    style_lines.append( [] )
                    plot_lines.append( [] )
                    qm_min_list.append( [] )
		    bar_plot_lines.append( [] )
		    z_plot_line.append( [] )
		    
		    
		    # Initialize max min data file 
		    minmax_data = struct_dir +'/' + mol_acc + "_" + mol_id + options.qm_sufix + '_bar.dat'
		    b_dat = open(minmax_data ,"w")			
		    b_dat.write("# mol_i , mol_repeat , dih_id, qm_min, ang_min, qm_max, ang_max , qm_bar " )
		    b_dat.close()
		    
		    zerobar_data = struct_dir +'/' + mol_acc + "_" + mol_id + options.qm_sufix + '_bar.dat'
		    b_dat = open(zerobar_data ,"w")	
    		    b_dat.write("# mol_repeat , d_indx - n_max, qm_bar " )
		    b_dat.close()		    
		    z_plot_line[m_indx] = " \'"+zerobar_data+"\' " + ' us 1:3  with linespoints  title \' even n \' '
		    
		    
		    zerobar_data = struct_dir +'/' + mol_acc + "_" + mol_id + options.qm_sufix + 'odd_bar.dat'		    
		    b_dat = open(zerobar_data ,"w")	
    		    b_dat.write("# mol_repeat , d_indx - n_max, qm_bar " )
		    b_dat.close()
		    z_plot_line[m_indx] += ", \'"+zerobar_data+"\' " + ' us 1:3  with linespoints title \' odd n  \'  '
		    
		    
		    
		    
	calc_i = 0
	qm_minmax = []
        for line in Lines:
            col = line.split()

            if( len(col) >= 4 and col[0] == "qm_dih" ):
	    
		mol_dir = col[1].strip()
		mol_id = col[2].strip()
		mol_repeat = int(col[3])
		mol_acc = col[4].strip()
		dih_id = col[5].strip()
		
		# File info
		struct_dir = mol_dir + "/" + mol_id + "/"
		job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)
		
		#a_k,a_i, a_j ,a_l,
		cent_min  = int(  col[10].strip() )
		cent_max = int(  col[11].strip() )
		cent_step = int(  col[12].strip() )
                    
                for m_indx in range(len(mol_list)):
                    if( mol_id == mol_list[m_indx] ):		
			mol_i = m_indx
			break
                        
                    
                    
                if( options.verbose ):
                    print "  getting plot lines " ,mol_id, mol_i,struct_dir,job_name,dih_id

		dih_qm = struct_dir +'/' +job_name + '-' + dih_id + options.qm_sufix + '.dat'

                # Find qm_min
		
		qm_min, ang_min, qm_max, ang_max = qm_analysis(  dih_qm )
		
                print " qm min ",qm_min
                calc_i += 1 
                
                
                style_l =  'set style line ' + str(calc_i+1) + ' lt  ' +  str(calc_i) + ' lw 1 '+ ' lc ' +  str(calc_i) + "\n"
                plot_l =  " \'"+dih_qm+"\' " + ' us 2:($3- ' + str(qm_min) + ') w l ls '+str(calc_i+1)+'  title '+ "\' " + mol_id+"n="+str(mol_repeat)+" dihedral "+dih_id+" \' smooth unique  , \\" + "\n"
        
                style_lines[mol_i].append(  style_l )
                plot_lines[mol_i].append(  plot_l )
                
			    
		qm_bar = qm_max - qm_min

		d_o = -1*(mol_repeat-2)
		
		dih_cnt = bar_cnt[mol_i,mol_repeat]
		d_indx = n_max + d_o + dih_cnt*2
		
		bar_grid[mol_repeat,d_indx ] = qm_bar
		
                if( options.verbose ):
		    print "   grid point ",mol_repeat, d_o + dih_cnt*2 ,d_indx
		
		minmax_data = struct_dir +'/' + mol_acc + "_" + mol_id +"-"+ str(mol_repeat)+ options.qm_sufix + '_bar.dat'
		if( dih_cnt == 0 ):
		    b_dat = open(minmax_data ,"w")			
		    b_dat.write("# mol_repeat ; d_indx , qm_bar " )
		    b_dat.close()
		    
		    bplot_l =  " \'"+minmax_data+"\' " + ' us 2:3  '+' with linespoints ls '+str(int(mol_repeat)+1)+'  title '+ "\' " + mol_id+"n="+str(mol_repeat) +" \'  , \\" + "\n"
		    bar_plot_lines[mol_i].append( bplot_l )
		    
		b_dat = open(minmax_data ,"a")			
		b_dat.write("\n %d %d %f" % ( mol_repeat , d_indx - n_max, qm_bar ) )
		b_dat.close()
		
		# Write
		if( (d_indx - n_max) == 0 ):
		    zerobar_data = struct_dir +'/' + mol_acc + "_" + mol_id + options.qm_sufix + '_bar.dat'
		    b_dat = open(zerobar_data ,"a")	
    		    b_dat.write("\n %d %d %f" % ( mol_repeat , d_indx - n_max, qm_bar ) )
		    b_dat.close()
		    
		    #z_plot_line[mol_i] += " \'"+zerobar_data+"\' " + ' us 1:3  with linespoints  '
		    
		    #b_dat = open(minmax_data ,"w")			
		    #b_dat.write("# mol_repeat ; d_indx , qm_bar " )
		    #b_dat.close()

		# Write
		if( (d_indx - n_max) == 1 ):
		    zerobar_data = struct_dir +'/' + mol_acc + "_" + mol_id + options.qm_sufix + 'odd_bar.dat'
		    b_dat = open(zerobar_data ,"a")	
    		    b_dat.write("\n %d %d %f" % ( mol_repeat , d_indx - n_max, qm_bar ) )
		    b_dat.close()
		    
		    #z_plot_line[mol_i] += " \'"+zerobar_data+"\' " + ' us 1:3  with linespoints  '
		    		    
		bar_cnt[mol_i,mol_repeat] += 1
		
		

            if( len(col) >= 4 and col[0] == "ff_dih"  ):
                
                struct_dir = col[1].strip()
                job_name = col[2].strip()
                ff_software = col[3].strip()
                ff_type_id = col[4].strip()
    
                dih_id = col[5].strip() 
                a_k  = int(  col[6].strip() )
                a_i  = int(  col[7].strip() )
                a_j  = int(  col[8].strip() )
                a_l  = int(  col[9].strip() )
                cent_min  = int(  col[10].strip() )
                cent_max = int(  col[11].strip() )
                cent_step = int(  col[12].strip() )
                    
                for m_indx in range(len(mol_list)):
                    mol_id_ref = mol_list[m_indx]
                    
                    for calc_indx in range(len(method_list[m_indx])):
                        method_id = method_list[m_indx][calc_indx]
                        mol_dir = moldir_list[m_indx][calc_indx]
                        mol_repeat = repeat_list[m_indx][calc_indx]
                        
                        struct_dir_ref = method_id + "/" + mol_dir + "/"
                        job_name_ref = mol_id_ref + "_n" + str(mol_repeat)
                        
                        
                        if( struct_dir == struct_dir_ref and job_name == job_name_ref ):
                            mol_id = mol_id_ref
                            mol_i = m_indx
                            calc_i = calc_indx
                            repeat_n = mol_repeat
                            m_id = method_id
                            break
                        
                    
                    
                if( options.verbose ):
                    print "  getting plot lines " ,mol_id, mol_i,calc_i,struct_dir,job_name,dih_id

		#dih_qm = struct_dir +'/' +job_name + '-' + dih_id + options.qm_sufix + '.dat'
                dih_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + ff_software + ff_type_id +".dat"

                # Find qm_min
                ff_min = 100000.0 
                f_ff = open(dih_ff ,"r")
                f_ff_Lines = f_ff.readlines()
                f_ff.close()
                for f_ff_line in f_ff_Lines:
                    f_ff_col = f_ff_line.split()
                
                    if( len(f_ff_col) >= 5 and f_ff_col[0] != "#" ):
                        # constrained optimization energy 
                        ff_energy_c = float( f_ff_col[5] )
                        if( ff_energy_c < ff_min ): ff_min = ff_energy_c
                        
                        
                print " ff min ",ff_min
                
                
                
                style_l =  'set style line ' + str(calc_i+1) + ' lt  ' +  str(calc_i) + ' lw 5 '+ ' lc ' +  str(calc_i) + "\n"
                plot_l =  " \'"+dih_ff+"\' " + ' us 2:($6- ' + str(ff_min) + ')*KCEV w l ls '+str(calc_i+1)+'  title '+ "\' " + m_id+"n="+str(repeat_n)+" dihedral "+dih_id+" \' smooth unique  , \\" + "\n"
        
                style_lines[mol_i].append(  style_l )
                plot_lines[mol_i].append(  plot_l )
                

    for m_indx in range(len(mol_list)):
        mol_id = mol_list[m_indx]
        if( options.verbose ):
            print "  Plotting " ,mol_id

        style_i = ""
        for s_line in style_lines[m_indx] :
            style_i += s_line
            
        plot_i = ""
        for p_line in plot_lines[m_indx] :
            plot_i += p_line
    
	plt_template = 'dih_plt.template'
	f = open( plt_template , 'r')
	plt_file = f.read()
	f.close()
	
        plt_file = replace(plt_file,'<style_lines>',style_i)
        plt_file = replace(plt_file,'<plot_lines>',plot_i)
        plt_file = replace(plt_file,'<structure_name>',mol_id)
	    
        plt_dir_file = mol_id+'.plt'
        f = open( plt_dir_file , 'w')
        f.write(plt_file)
        f.close()
	
	plt_template = 'dih_plt.template'
	f = open( plt_template , 'r')
	plt_file = f.read()
	f.close()
	
	
	plt_template = 'dih_b_plt.template'
	f = open( plt_template , 'r')
	plt_file = f.read()
	f.close()
	
        plot_i = ""
        for p_line in bar_plot_lines[m_indx] :
            plot_i += p_line
	
        plt_file = replace(plt_file,'<style_lines>',style_i)
        plt_file = replace(plt_file,'<plot_lines>',plot_i)
        plt_file = replace(plt_file,'<structure_name>',mol_id)
        plt_file = replace(plt_file,'<z_plot_lines>',z_plot_line[m_indx] )
	
	
	
        plt_dir_file = mol_id+'_bar.plt'
        f = open( plt_dir_file , 'w')
        f.write(plt_file)
        f.close()
	
	# Write data file of barriers
        bar_data = mol_id+'_bar.grid'
        f = open( bar_data , 'w')
	f.write( "  repeat (n) ; dih postion ; barrier value (eV) " )
	
	for mol_repeat in range( n_max ):
	    f.write("\n")
	    for dih_indx in range( n_max*2 ):
		b_val = bar_grid[mol_repeat,dih_indx]
		f.write( "\n %d %d %f " %  ( mol_repeat, dih_indx, b_val) )
		
        f.close()

if __name__=="__main__":
    main()

