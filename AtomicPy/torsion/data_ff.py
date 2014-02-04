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
    

    parser.add_option("-v","--verbose", dest="verbose", default=True, help="Verbose output ")

    # Cluster options
    parser.add_option("--cluster_host", dest="cluster_host",type="string",default="peregrine",help=" name of cluster ")


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
    elif( options.cluster_host == "redmesa" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/g09/C.01"

    return options, args


def main():
    import sys, os , string 
    from string import replace
    import file_io, elements, gaussian, gromacs , lammps , xmol 

    options, args = get_options()
    
    # Store working dir  
    work_dir = os.getcwd()

    if( options.cluster_host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	md_min_template = "peregrine.md_min.template"
	lammps_min_template = "peregrine.lmp_min.template"
    elif( options.cluster_host == "redmesa" ):
	load_gaussian = "module load gaussian/g09/C.01"
        user = 'twkempe'
	md_min_template = "redmesa.md_min.template"
	lammps_min_template = "redmesa.lmp_min.template"

    if( options.ff_software == "gromacs"):	
	f = open(md_min_template,'r')
	md_min_templ = f.read()
	f.close()
    elif( options.ff_software == "lammps" ):	
	f = open(lammps_min_template,'r')
	md_min_templ = f.read()
	f.close()
    
    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
        for line in Lines:
            col = line.split()
            if( len(col) >= 4 and col[0] == "ff_dih" ):
		
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
		
		# Open output file 
		dih_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + ff_software + ff_type_id +".dat"
		if( options.verbose ):
		    print "  Writing data file ",dih_ff
		ff_out = open(dih_ff,'w')
		ff_out.write( '# cent_indx,cent_angle,ff_sp_full (eV) ,ff_sp_d0 (eV) ,ff_energy_f (eV) ,ff_energy_c (eV) ,n_angles,dih_angles ')
		
		# xmol 
		xmol_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + ff_software + ff_type_id +".xmol"
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
                zmat_id = job_name + "-ZMAT"
		fchk_file = struct_dir +'/' + zmat_id +"/"+zmat_id+".fchk"
		NA, ELN, R, TOTAL_ENERGY , Q_ESP   = gaussian.parse_fchk( fchk_file )
		ASYMB = elements.eln_asymb(ELN)
		    
		# Loop over angels of central dihedrals
		cent_indx = 0
		print "     looping over ",cent_min,cent_max,cent_step
		for cent_angle in range(cent_min,cent_max,cent_step):
		    print cent_angle
		    
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
				
				if( run_min_f == 0 and run_min_c == 1 ):
				    #
				    ff_energy_f = gromacs.get_logenergy('min_f')
				    ff_energy_c = gromacs.get_logenergy('min_c')
				    #
				    ##ff_r = gromacs.gro_coord('sp_2')
				    ff_r = gromacs.get_coord('min_f',options)
				    dih_angles = [] 
				    for indx_fs in range( len( DIH_fitset )):
					angle_indx = DIH_fitset[indx_fs]
					angle = gromacs.get_dihangle('min_f',angle_indx,options)
					dih_angles.append( angle )
					
				    os.chdir(work_dir)
				    #ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f " % ( cent_indx,cent_angle, ff_energy_min_1,ff_energy_min_2,ff_energy_sp ))
				    ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f " % ( cent_indx,cent_angle,dih_angles[0],dih_angles[1],dih_angles[2],dih_angles[3],ff_energy_f,ff_energy_c ))
				    xmol.print_xmol(ASYMB,ff_r,xmol_ff)
	    
			    else:
				print ' error in ff input files '
			
			elif( ff_software == "lammps" ):
			    rest_file = "rest.log"
			    fit_file = "fit.log"
			    sp_file = "full_sp.log"
			    d0_file = "d0_sp.log"
			    if( file_io.file_exists(rest_file) and file_io.file_exists(fit_file) and file_io.file_exists(sp_file) and file_io.file_exists(d0_file) ):
				ff_sp_full = lammps.get_pe(sp_file)
				ff_sp_d0 = lammps.get_pe(d0_file)
				ff_energy_c = lammps.get_pe(rest_file)
				ff_energy_f = lammps.get_pe(fit_file)
				ff_r = lammps.last_xmol('fit.xyz',options)
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
				print  cent_indx,cent_angle,ff_sp_full,ff_sp_d0,ff_energy_f,ff_energy_c,n_angles,dih_angles 
				ff_out.write( " \n %8d %8.4f %16.8f %16.8f %16.8f %16.8f %6d %s" % ( cent_indx,cent_angle,ff_sp_full,ff_sp_d0,ff_energy_f,ff_energy_c,n_angles,dih_angles ))
				#ff_energy.append([cent_indx,cent_angle,dih_angles[0],dih_angles[1],dih_angles[2],dih_angles[3],ff_sp_full,ff_sp_d0,ff_energy_f,ff_energy_c ] )
				xmol.print_xmol(ASYMB,ff_r,xmol_ff)
				
			os.chdir(work_dir)
	
		ff_out.close()

if __name__=="__main__":
    main() 

