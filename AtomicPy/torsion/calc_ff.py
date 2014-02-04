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

    # Cluster options
    parser.add_option("--cluster_host", dest="cluster_host",type="string",default="peregrine",help=" name of cluster ")


    parser.add_option("--pmem",dest="pmem", type="int", default="1700",help=" Memory per processor ")
    parser.add_option("--npros", dest="npros", type="int", default="16",help=" Number of processors ")
    parser.add_option("--nnodes", dest="nnodes", type="int", default="1",help=" Number of nodes ")


    parser.set_defaults(submit=False)
    parser.add_option("--submit", dest="submit",action="store_true",help=" submit calculations to the queue ")
    parser.set_defaults(localrun=False)
    parser.add_option("--localrun", dest="localrun",action="store_true",help=" Run calculations locally, ment to use if scripts included in submitted job")
 
    parser.set_defaults(submit_command="qsub")
    parser.add_option("--submit_command", dest="submit_command",type="string",help=" command used to submit script to the queue ")


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
    if( options.cluster_host == "peregrine" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"

    return options, args

def main():
    import sys, os , string 
    from string import replace
    import file_io, gaussian, cluster, gromacs , lammps 

    options, args = get_options()

    # Store working dir  
    work_dir = os.getcwd()
    
    if( options.cluster_host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	md_min_template = "peregrine.md_min.template"
	lammps_min_template = "peregrine.lmp_min.template"
	
	lammps_dir= options.lammp_dir 
	lammps_src="lmp_serial"
	
    if( options.cluster_host == "macbook" ):
        user = 'tkemper'
	lammps_dir= '/Users/'+user+'/Software/lammps-24Apr13/src/'
	lammps_src="lmp_mactk"
				
	    		    
     
    # Read index files from args
    for indx_file in args:
	print " Reading ",indx_file
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
        for line in Lines:
            col = line.split()
            if( len(col) >= 4 and col[0] == "ff_dih" ):
		
		print " ff_dh  line ",line
		
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
		
			    
		if( ff_software == "gromacs"):
		    
		    f = open(md_min_template,'r')
		    md_min_templ = f.read()
		    f.close()
			    
		if( ff_software == "lammps" ):
		    f = open(lammps_min_template,'r')
		    md_min_templ = f.read()
		    f.close()
			    
			    
		# Loop over angels of central dihedrals
		cent_indx = 0
		for cent_angle in range(cent_min,cent_max,cent_step):
		    ff_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_ff' + ff_software + ff_type_id
		    ff_dir = struct_dir + ff_id 

	    
		    if( os.path.isdir(ff_dir) ):                    
			# Get ff energy
			os.chdir(ff_dir)
			if( ff_software == "gromacs"):
	    

			    g_top = 'out_const.top'
			    g_gro = 'out.gro'
			    input_correct = gromacs.check_input( g_gro,g_top,options )
			    
			    if( input_correct ):
				g_file = "min_f.log"
				run_min_f = gromacs.check_g(g_file)
				g_file = "min_c.log"
				run_min_c = gromacs.check_g(g_file)                
				
				if( run_min_f or run_min_c ):
	    
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
					pbs_dih = replace(pbs_dih,"<input_file>",input_file)
					pbs_dih = replace(pbs_dih,"<pmem>",str(options.pmem))
					pbs_dih = replace(pbs_dih,"<npros>",str(options.npros))
					f = file(pbs_name, "w")
					f.write(pbs_dih)
					f.close()
					
					#submit_job( mols_dir_ff, pbs_id,options )
					
				    elif ( options.localrun):
					
					g_top = 'out_fit.top'
					g_gro = 'out.gro'
					g_mdp = 'min_f.mdp'
					s_suf = options.gromacs_sufix + '_d'
					gromacs.print_min(g_mdp)
					ff_energy_f = gromacs.run_gromacs(g_gro,g_top,g_mdp,s_suf,options )
	    
					g_top = 'out_const.top'
					g_gro = 'min_f.gro'
					g_mdp = 'min_c.mdp'
					s_suf =  options.gromacs_sufix + '_d'
					gromacs.print_min(g_mdp)
					ff_energy_c = gromacs.run_gromacs(g_gro,g_top,g_mdp,s_suf,options )
			    else:
				print ' error in ff input files '
				
			elif( ff_software == "lammps" ):
						
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

