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

    parser.add_option("--pmem",dest="pmem", type="int", default="1700",help=" Memory per processor ")
    parser.add_option("--npros", dest="npros", type="int", default="16",help=" Number of processors ")
    parser.add_option("--nnodes", dest="nnodes", type="int", default="1",help=" Number of nodes ")

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

    # Set options based on cluster 
    if( options.host == "peregrine" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"

    return options, args

def main():
    import sys, os , string 
    from string import replace
    import jsonapy
    import file_io, gaussian, cluster, gromacs , lammps 

    options, args = get_options()

    # Store working dir  
    work_dir = os.getcwd()
    
    if( options.host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	md_min_template = "peregrine.md_min.template"
	lammps_min_template = "peregrine.lmp_min.template"
	
	lammps_dir= options.lammp_dir 
	lammps_src="lmp_serial"

    if( options.host == "dale" ):
	load_gaussian = "module load gaussian"
        user = 'tkemper'
	
	
	md_min_template = "dale.gromacs_min.template"
	lammps_min_template = "dale.lmp_min.template"
	
	load_gromacs = 'module load gromacs/4.6.1'
	
	lammps_dir= options.lammp_dir 
	lammps_src="lmp_serial"
		
    if( options.host == "macbook" ):
        user = 'tkemper'
	lammps_dir= '/Users/'+user+'/Software/lammps-24Apr13/src/'
	lammps_src="lmp_mactk"
				
	    		
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
			
			print "   Running job ",struct_dir,job_name,dih_id
			
			if( options.submit and options.ff_software == "gromacs"):
			    
			    f = open(md_min_template,'r')
			    md_min_templ = f.read()
			    f.close()
				    
			if( options.submit and options.ff_software == "lammps" ):
			    f = open(lammps_min_template,'r')
			    md_min_templ = f.read()
			    f.close()
				    
					    
			# Loop over angels of central dihedrals
			cent_indx = 0
			for cent_angle in range(cent_min,cent_max,cent_step):
			    ff_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_ff' + options.ff_software + ff_type_id
			    ff_dir = struct_dir + ff_id 
			    
			    print ff_dir
		    
			    if( os.path.isdir(ff_dir) ):                    
				# Get ff energy
				os.chdir(ff_dir)
				if( options.ff_software == "gromacs"):
		    
	
				    g_top = 'full.top'
				    g_gro = 'full.gro'
				    input_correct = gromacs.check_input( g_gro,g_top,load_gromacs,options.gromacs_sufix,options.gromacs_dir  )
				    
				    if( input_correct ):
					g_file = "full_sp.log"
					run_fsp = gromacs.check_g(g_file)
					g_file = "dih0_sp.log"
					run_d0sp = gromacs.check_g(g_file)                
					
					if( run_fsp or run_d0sp ):
		    
					    if( options.submit ):
						# Print mdp files 
						g_mdp = 'min.mdp'
						gromacs.print_min(g_mdp)                        
						g_mdp = 'sp.mdp'
						gromacs.print_sp(g_mdp)
					    
						md_min_id =  ff_id+'md_min'
						pbs_id = md_min_id+'.pbs'
						pbs_name = pbs_id
						pbs_dih = md_min_templ
						pbs_dih = replace(pbs_dih,"<calc_id>",md_min_id)
						# pbs_dih = replace(pbs_dih,"<input_file>",input_file)
						pbs_dih = replace(pbs_dih,"<pmem>",str(options.pmem))
						pbs_dih = replace(pbs_dih,"<nnodes>",str(options.nnodes))
						f = file(pbs_name, "w")
						f.write(pbs_dih)
						f.close()
							
						mols_dir_ff = ""
						#submit_job( mols_dir_ff, pbs_id,options )
						cluster.submit_job( mols_dir_ff, pbs_id ,options )
							
						#submit_job( mols_dir_ff, pbs_id,options )
						
					    elif ( options.localrun):
						
						g_top = 'full.top'
						g_gro = 'full.gro'
						g_mdp = 'full_sp.mdp'
						s_suf = options.gromacs_sufix + '_d'
						gromacs.print_sp(g_mdp)
						ff_energy_fsp = gromacs.run_gromacs(g_gro,g_top,g_mdp,load_gromacs,options.gromacs_sufix,options.gromacs_dir )
		    
						g_top = 'dih0.top'
						g_gro = 'full.gro'
						g_mdp = 'dih0_sp.mdp'
						s_suf =  options.gromacs_sufix + '_d'
						gromacs.print_sp(g_mdp)
						ff_energy_d0sp = gromacs.run_gromacs(g_gro,g_top,g_mdp,load_gromacs,options.gromacs_sufix,options.gromacs_dir )
				    else:
					print ' error in ff input files '
					
				elif( options.ff_software == "lammps" ):
				    
				    print "  Running lammps ",ff_id
							
				    if( options.submit ):
					
	
					md_min_id =  ff_id+'md_min'
					pbs_id = md_min_id+'.pbs'
					pbs_name = pbs_id
					pbs_dih = md_min_templ
					pbs_dih = replace(pbs_dih,"<calc_id>",md_min_id)
					pbs_dih = replace(pbs_dih,"<pmem>",str(options.pmem))
					pbs_dih = replace(pbs_dih,"<npros>",str(options.npros))
					pbs_dih = replace(pbs_dih,"<run_dir>",ff_dir)
					f = file(pbs_name, "w")
					f.write(pbs_dih)
					f.close()
					mols_dir_ff = ""
					#submit_job( mols_dir_ff, pbs_id,options )
					cluster.submit_job( mols_dir_ff, pbs_id ,options )
					
				    elif ( options.localrun ):
					
		    
					lmp_in = "full_sp.in"
					lammps.run_lmp(lammps_dir,lammps_src,lmp_in)
					
					lmp_in = "d0_sp.in"
					lammps.run_lmp(lammps_dir,lammps_src,lmp_in)
					
					lmp_in = "rest.in"
					lammps.run_lmp(lammps_dir,lammps_src,lmp_in)
					
					lmp_in = "fit.in"
					lammps.run_lmp(lammps_dir,lammps_src,lmp_in)
					
			    os.chdir(work_dir)


if __name__=="__main__":
    main() 

