# Plot torsional data based on reference file 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=True, help="Verbose output ")

    parser.add_option("--update_img", dest="update_img",default=True,help=" Update images ")

    # should be reference file 
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    
    (options, args) = parser.parse_args()
    
    return options, args

def main():
    import sys, os 
    import string  
    import file_io
    from string import replace

    options, args = get_options()

    mol_list = []
    method_list = []
    moldir_list = []
    repeat_list = []
    
    style_lines = []
    plot_lines = []
    qm_min_list  = []

    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
        for line in Lines:
            col = line.split()
            if( len(col) >= 4 and col[0] == "gen" ):
                mol_id = col[1].strip()
                main_dir = col[2].strip()
                mol_dir = col[3].strip()
                mol_repeat = int(col[4])
                
                # File info
                struct_dir = main_dir + "/" + mol_dir + "/"
                job_name = mol_id + "_n" + str(mol_repeat)

                new_mol = 1
                for m_indx in range(len(mol_list)):
                    if( mol_id == mol_list[m_indx] ):
                        new_mol = 0
                        method_list[m_indx].append( main_dir )
                        moldir_list[m_indx].append( mol_dir )
                        repeat_list[m_indx].append( mol_repeat )
                        
                if( new_mol ):
                    mol_list.append( mol_id)
                    m_indx = len(mol_list)  - 1
                    method_list.append( [] )
                    method_list[m_indx].append( main_dir )
                    moldir_list.append( [] )
                    moldir_list[m_indx].append( mol_dir )
                    repeat_list.append( [] )
                    repeat_list[m_indx].append( mol_repeat )
                    
                    style_lines.append( [] )
                    plot_lines.append( [] )
                    qm_min_list.append( [] ) 
                
        for line in Lines:
            col = line.split()

            if( len(col) >= 4 and col[0] == "qm_dih" ):
                struct_dir = col[1].strip()
                job_name = col[2].strip()
		dih_id = col[3].strip()
		#a_k,a_i, a_j ,a_l,
                a_k  = int(  col[4].strip() )
                a_i  = int(  col[5].strip() )
                a_j  = int(  col[6].strip() )
                a_l  = int(  col[7].strip() )
		cent_min  = int(  col[8].strip() )
		cent_max = int(  col[9].strip() )
		cent_step = int(  col[10].strip() )
                    
                    
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
		str_file = struct_dir +'/' + job_name+'-'+dih_id+"_qm.xmol"
		str_type = 'xyz'
		snapshot_file = dih_id+'_qm'

		style_2_select = 'index '+ str(a_k) +' or index '+ str(a_i) +' or index '+ str(a_j)  +' or index '+ str(a_l)

		vmd_template = 'dih_vmd.template'
		f = open( vmd_template , 'r')
		vmd_file = f.read()
		f.close()

		vmd_file = replace(vmd_file,'<str_file>',str_file)
		vmd_file = replace(vmd_file,'<str_type>',str_type)
		vmd_file = replace(vmd_file,'<style_2_select>',style_2_select)
		vmd_file = replace(vmd_file,'<snapshot_file>',snapshot_file)
		
	    
		vmd_in = dih_id+'qm.vmd'
		f = open( vmd_in , 'w')
		f.write(vmd_file)
		f.close()
	    
		
		if( options.update_img ):
		    make_snapshot = "vmd_MACOSXX86 < " + vmd_in
		    os.system(make_snapshot)
		#if( open_files ):
		#    open_tga = " open "+ snapshot_file+'.pdf' 
		#    os.system(open_tga)
	    
			    

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


		str_file = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + ff_software + ff_type_id +".xmol"
		str_type = 'xyz'
		snapshot_file = dih_id+'_ff'

		style_2_select = 'index '+ str(a_k) +' or index '+ str(a_i) +' or index '+ str(a_j)  +' or index '+ str(a_l)

		vmd_template = 'dih_vmd.template'
		f = open( vmd_template , 'r')
		vmd_file = f.read()
		f.close()

		vmd_file = replace(vmd_file,'<str_file>',str_file)
		vmd_file = replace(vmd_file,'<str_type>',str_type)
		vmd_file = replace(vmd_file,'<style_2_select>',style_2_select)
		vmd_file = replace(vmd_file,'<snapshot_file>',snapshot_file)

		vmd_in = dih_id+'ff.vmd'
		f = open( vmd_in , 'w')
		f.write(vmd_file)
		f.close()
	    
		
		if( options.update_img ):
		    make_snapshot = "vmd_MACOSXX86 < " + vmd_in
		    os.system(make_snapshot)
		    
if __name__=="__main__":
    main()

