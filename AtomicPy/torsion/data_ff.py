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
    parser.add_option("--userid", dest="userid",type="string",default="tkemper",help=" name of user  ")

    # How to run the needed calculations 
    parser.add_option("--submit", dest="submit",action="store_true", default=False,help=" submit calculations to the queue ")
    parser.add_option("--localrun", dest="localrun",action="store_true", default=False,help=" Run calculations locally")
    parser.add_option("--submit_command", dest="submit_command",type="string", default="qsub",help=" command used to submit script to the queue ")

    parser.add_option("--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")


    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    
    # Force field generation options     
    parser.set_defaults(ff_software="gromacs")
    parser.add_option("--ff_software", dest="ff_software",type="string",help=" what software to use for the ff calculations   ")

    # Gromacs related options 
    parser.add_option("--gromacs_dir", dest="gromacs_dir",type="string",default="",help=" Directory of gromacs run files   ")
    parser.add_option("--gromacs_sufix", dest="gromacs_sufix",type="string",default="",help=" sufix for gromacs such as _d or _mpi  ")
    
    # Lammps options 
    parser.set_defaults(lammp_dir="$HOME/Software/lammps/src/")
    parser.add_option("--lammp_dir", dest="lammp_dir",type="string",help=" Directory of lammps run files  ")

    (options, args) = parser.parse_args()


    return options, args


def main():
    import sys, os , string 
    from string import replace
    import jsonapy
    import file_io, elements, gaussian, gromacs , lammps , xmol 

    options, args = get_options()

	    
    # Store working dir  
    work_dir = os.getcwd()
    
    
    # Set options based on cluster 
    if( options.host == "peregrine" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"
    elif( options.host == "dale" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian "
	
	load_gromacs = 'module load gromacs/4.6.1'
	options.gromacs_sufix = "_mpi"
	
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
		dih_id_list ,cent_min_list ,cent_max_list ,cent_step_list,a_k_list, a_i_list, a_j_list, a_l_list,ff_type_list,fftor_found = jsonapy.read_ff_tor(json_data)
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found and fftor_found  ):
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
		    
		    rec_qm = job_name + ".rec"	
			
		    for dih_indx in range( len(dih_id_list) ):
			dih_id = dih_id_list[dih_indx]
			cent_min = cent_min_list[dih_indx]
			cent_max = cent_max_list[dih_indx]
			cent_step = cent_step_list[dih_indx]
			a_k = a_k_list[dih_indx]
			a_i = a_i_list[dih_indx]
			a_j = a_j_list[dih_indx]
			a_l = a_l_list[dih_indx]
			ff_type_id = ff_type_list[dih_indx]
			
			print "   Getting energies for ",struct_dir,job_name,dih_id
			
			# Open output file 
			dih_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + options.ff_software + ff_type_id +".dat"
			if( options.verbose ):
			    print "  Writing data file ",dih_ff
			    
			ff_out = open(dih_ff,'w')
			ff_out.write( '# cent_indx,cent_angle,ff_sp_full (eV) ,ff_sp_d0 (eV) ,ff_energy_f (eV) ,ff_energy_c (eV) ,n_angles,dih_angles ')
			
			# xmol 
			xmol_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + options.ff_software + ff_type_id +".xmol"
			if( options.verbose ):
			    print "  Writing xmol file ",xmol_ff
			if( file_io.file_exists( xmol_ff ) ): os.remove(xmol_ff)
			
			# Read in all components of fitted dihedral 
			dih_comp_name = struct_dir +'/' + job_name + "_fit.list"
			f_fit = open(dih_comp_name,'r')
			Lines_fit = f_fit.readlines()
			f_fit.close()
	
			DIH_fitset = []
			for line_fit in Lines_fit:
			    col_fit = line_fit.split()
			    if( len(col_fit) >= 4 and col_fit[0] == "dihedral" ):
				DIH_fitset.append( [int( col_fit[1])-1,int( col_fit[2])-1,int( col_fit[3])-1,int( col_fit[4])-1] )
				
			# Get elements form -ZMAT file, should get from out.dat
	
			# Get geometry from ZMAT
			zmat_fchk = "%s/%s%s/%s%s" % ( struct_dir, job_name , "-ZMATOPT"  ,job_name,"-ZMATOPT.fchk" )
			
			print " Checking for complete zmatrix optimiztion ",zmat_fchk
			zmat_finished = file_io.file_exists( zmat_fchk )
			
			
			if( zmat_finished  ):
			    
			    NA, ELN, R, TOTAL_ENERGY , Q_ESP   = gaussian.parse_fchk( zmat_fchk )
			    ASYMB = elements.eln_asymb(ELN)
				
			    # Loop over angels of central dihedrals
			    cent_indx = 0
			    print "     looping over ",cent_min,cent_max,cent_step
			    for cent_angle in range(cent_min,cent_max,cent_step):
				print cent_angle
			
			    calc_success = 0
			    for cent_angle in range(cent_min,cent_max,cent_step):
				ff_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_ff' + options.ff_software + ff_type_id
				ff_dir = struct_dir + ff_id
	    
				if( os.path.isdir(ff_dir) ):                    
				    # Get ff energy
				    os.chdir(ff_dir)
				    if( options.ff_software == "gromacs"):
			
					g_top = 'full.top'
					g_gro = 'full.gro'
					# input_correct = gromacs.check_input( g_gro,g_top,options )
					input_correct = gromacs.check_input( g_gro,g_top,load_gromacs,options.gromacs_sufix,options.gromacs_dir  )
					
					if( input_correct ):
					    g_file = "full_sp.log"
					    run_fsp = gromacs.check_g(g_file)
					    g_file = "dih0_sp.log"
					    run_d0sp = gromacs.check_g(g_file)                
					    
					    if( run_fsp == 0 and run_d0sp == 0 ):
						#
						ff_sp_full = gromacs.get_g_potenergy('full_sp',load_gromacs,options.gromacs_sufix,options.gromacs_dir )
						
						
						# ff_sp_d0 = gromacs.get_logenergy('dih0_sp')
						ff_sp_d0 = gromacs.get_g_potenergy('dih0_sp',load_gromacs,options.gromacs_sufix,options.gromacs_dir )
						
						#
						# ff_r = gromacs.gro_coord('sp_2')
						#
						
						#
						# copy full.gro to full_sp.gro 
						#
						cp_gro = " cp full.gro full_sp.gro "
						os.system(cp_gro)
						
						ff_r = gromacs.get_coord('full_sp',load_gromacs,options.gromacs_sufix,options.gromacs_dir )

						dih_angles = ""
						n_angles = 0
						for indx_fs in range( len( DIH_fitset )):
						    angle_indx = DIH_fitset[indx_fs]
						    angle = gromacs.get_dihangle('full_sp',angle_indx,load_gromacs,options.gromacs_sufix,options.gromacs_dir )
						    dih_angles += " "+str(angle)
						    n_angles += 1
						    
						#
						# Set relaxed energies to zero as the relaxations are not working properly yet
						#
						ff_rlx_full = 0.0 
						ff_rlx_d0 = 0.0
						
						os.chdir(work_dir)
						#ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f " % ( cent_indx,cent_angle, ff_energy_min_1,ff_energy_min_2,ff_energy_sp ))
						# ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f " % ( cent_indx,cent_angle,dih_angles[0],dih_angles[1],dih_angles[2],dih_angles[3],ff_energy_f,ff_energy_c ))
						ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f %16.8f %6d %s" % ( cent_indx,cent_angle,ff_sp_full,ff_sp_d0,ff_rlx_full,ff_rlx_d0,n_angles,dih_angles ))
						xmol.print_xmol(ASYMB,ff_r,xmol_ff)
			
					else:
					    print ' error in ff input files '
				    
				    elif( options.ff_software == "lammps" ):
					sp_file = "full_sp.log"
					d0_file = "dih0_sp.log"
					if( file_io.file_exists(sp_file) and file_io.file_exists(d0_file) ):  #and file_io.file_exists(sp_file) and file_io.file_exists(d0_file) ):
					    ff_sp_full = lammps.get_pe(sp_file)
					    ff_sp_d0 = lammps.get_pe(d0_file)
					    
					    #ff_energy_c = lammps.get_pe(rest_file)
					    #ff_energy_f = lammps.get_pe(fit_file)
					
					    #
					    # Set relaxed energies to zero as the relaxations are not working properly yet
					    #
					    ff_rlx_full = 0.0 
					    ff_rlx_d0 = 0.0
					    
					    ff_r = lammps.last_xmol('sp.xyz',options)
					    dih_angles = ""
					    n_angles = 0
					    for indx_fs in range( len( DIH_fitset )):
						angle_indx = DIH_fitset[indx_fs]
						angle = lammps.get_dihangle(ff_r,angle_indx,options)
						#dih_angles.append( angle )
						dih_angles += " "+str(angle)
						n_angles += 1
	    
					    os.chdir(work_dir)
					    #ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f " % ( cent_indx,cent_angle,dih_angles[0],dih_angles[1],dih_angles[2],dih_angles[3],ff_sp_full,ff_sp_d0,ff_energy_f,ff_energy_c ))
					    print   cent_indx,cent_angle,ff_sp_full,ff_sp_d0,ff_rlx_full,ff_rlx_d0,n_angles,dih_angles 
					    ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f %16.8f %6d %s" % ( cent_indx,cent_angle,ff_sp_full,ff_sp_d0,ff_rlx_full,ff_rlx_d0,n_angles,dih_angles ))
					    #ff_energy.append([cent_indx,cent_angle,dih_angles[0],dih_angles[1],dih_angles[2],dih_angles[3],ff_sp_full,ff_sp_d0,ff_energy_f,ff_energy_c ] )
					    xmol.print_xmol(ASYMB,ff_r,xmol_ff)
					    
					    calc_success = 1
					    
				    os.chdir(work_dir)
				
			    ff_out.close()
			    
			if( calc_success ):
			    
			    # make sure rec has not been previously recorded 
			    wrte_recline = 1
			    json_path = "%s/%s" % (work_dir, struct_dir)
			    
			    if( file_io.file_exists( rec_qm ) ):
				
				F = open(rec_qm,'r')
				Lines = F.readlines()
				F.close()
				# Check for a complete exicution
				for line in Lines:
				    col = line.split()
				    if( len(col) >= 5 and col[0] != "#" ):
					if( col[0].strip() == options.userid and col[1].strip() ==  options.host  and col[2].strip() == json_path  and col[3].strip() ==  struct_dir  and col[4].strip() ==  job_name ):
					    wrte_recline = 0 
				    
		    
			    if( wrte_recline ):
				    
				F = open(rec_qm,'a')
				F.write( " %s %s %s %s %s \n " % ( options.userid , options.host, json_path, struct_dir, job_name) )
				F.close()

			    
			    
			else:
			    print os.getcwd()
			    print " NEED ",zmat_finished," which should have been generated by mk_zmat.py "
			    
if __name__=="__main__":
    main() 

