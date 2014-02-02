#! /usr/bin/env python
# Make input files for torsional potential energy surface 

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
    parser.add_option("--qm_load", dest="qm_load",type="string",help=" string to load qm software module  ")

    parser.add_option("--qm_method", dest="qm_method", type="string",default="HF", help="Method of QM calculation ")
    parser.add_option("--qm_basis", dest="qm_basis", type="string",default="3-21G", help="Basis set of QM calculation ")
    parser.add_option("--qm_kywd", dest="qm_kywd", type="string",default="", help="Key words for QM calculation ")
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="",help=" sufix of qm data file  ")
    parser.add_option("--qm_charge", type="int",action="append", default="0",help="Input gaussain log file ")
    parser.add_option("--qm_mult", dest="qm_mult", type="int",default="0", help=" Shift in default spin multiplicity ( singlet,doublet) QM calculation, allows for triplets ")
    
    parser.add_option("--pmem",dest="pmem", type="int", default="1700",help=" Memory per processor ")
    parser.add_option("--npros", dest="npros", type="int", default="4",help=" Number of processors ")
    parser.add_option("--nnodes", dest="nnodes", type="int", default="1",help=" Number of nodes ")

    parser.add_option("--high_basis", dest="high_basis", type="string", default="cc-pVTZ",help=" Basis set for hihgh level energy calculations ")
    parser.add_option("--com_temp", dest="com_temp", type="string", default="esp.com.template",help=" Template for Links of dihedral calculation ")

    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")

    (options, args) = parser.parse_args()

    # Set options based on cluster 
    if( options.cluster_host == "peregrine" ):
        options.npros = options.nnodes* 24
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"
    elif( options.cluster_host == "redmesa" ):
        options.npros = options.nnodes *8
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/g09/C.01"
	    
	    
    return options, args

def main():
    import string, os , sys 
    # atomicpy
    import gaussian, elements, xmol , file_io , cluster 
    from string import replace
    
    options, args = get_options()
    
    # Set options
    pmem = options.pmem

    # Read in template files 
    f = open(options.dih_temp,'r')
    fix_templ = f.read()
    f.close()

    # 
    if( options.cluster_host == "peregrine" ):
	f = open("peregrine.pbs.template",'r') # redmesa.slurm.template
    elif( options.cluster_host == "redmesa" ):
	f = open("redmesa.slurm.template",'r')
    else:
        # HACK !!!
        f = open("peregrine.pbs.template",'r') # redmesa.slurm.template
	
    pbs_templ = f.read()
    f.close()

    # Read in index file produced by opv_generator or writen by hand
    #   Format:
    #     # - comment
    #     gen ; main dir ; mol dir ; accuracy ; n units
    #  Example:
    #     entry:
    #       gen  mol   D1_R2R200_A2_R3_  acc1   1
    #    will read log file:
    #      mol/D1_R2R200_A2_R3_/acc1_D1_R2R200_A2_R3__n1.fchk
    #    and make new files in
    #      mol/D1_R2R200_A2_R3_/
    
    
    # Store working dir  
    work_dir = os.getcwd()
    
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
                mol_repeat = int(col[3].strip() )
                mol_acc = col[4].strip()
                
                # File info
                struct_dir = mol_dir + "/" + mol_id + "/"
                job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)
                fchk_file = struct_dir + job_name + "/" + job_name + ".fchk"
		
                if( file_io.file_exists( fchk_file) ):
		    if( options.verbose ):
			
			print "    Optimization finished  ",fchk_file

		    # Run esp fit
		    calc_id_esp = job_name + "-ESP"
		    input_file = calc_id_esp + ".com"
		    log_file_esp = struct_dir +'/' + calc_id_esp +"/"+calc_id_esp+".log"
		    fchk_file_esp = struct_dir +'/' + calc_id_esp +"/"+calc_id_esp+".fchk"
		    
		    # Run esp fit if not finished
		    if( not file_io.file_exists( fchk_file_esp) ):
                        # Make esp input file
			if( options.verbose ):
			    print  '    Writing esp fit input file ',calc_id_esp," in ",struct_dir
			
					    
			os.chdir(struct_dir)
			gaussian.write_esp_com(calc_id_esp,ASYMB,R)
    
                        if( options.submit ):
                            if( options.verbose ):
                                print "     Submitting ESP optimization to queue "
                            # Print pbs script
                            pbs_id = calc_id_esp+'.pbs'
                            pbs_name = pbs_id
                            pbs_dih = pbs_templ
                            pbs_dih = replace(pbs_dih,"<calc_id>",calc_id_esp)
                            pbs_dih = replace(pbs_dih,"<input_file>",input_file)
                            pbs_dih = replace(pbs_dih,"<nnodes>",str( options.nnodes) )
                            pbs_dih = replace(pbs_dih,"<pmem>",str( options.pmem) )
                            pbs_dih = replace(pbs_dih,"<npros>",str( options.npros) )
                            f = file(pbs_name, "w")
                            f.write(pbs_dih)
                            f.close()
                            # Submit job
                            submit_job( struct_dir, pbs_id ,options )
                            
                        elif( options.localrun ):
                            if( options.verbose ):
				print  '    Running ESP fit '
                            gaussian.run(options, calc_id_esp)
			    
                        else:
                            print " Please mark either qsub or runloc options as True to run qm"	
			
			os.chdir(work_dir)
			
			
			
			
			
			
                read_fchk = 1 
                try:
                    with open(fchk_file) as f:
                        read_fchk = 1
                except IOError:
		    fchk_file = struct_dir + job_name + "-ZMAT/" + job_name +"-ZMAT"+ ".fchk"
                    if( options.verbose ):
                        print "    file  ",fchk_file," does not exist trying ZMAT file ",
			
		    try:
			with open(fchk_file) as f:
			    read_fchk = 1
		    except IOError:
			print " no fchk file found "
			sys.exit("no reference file ")
			
                    
                run_qm = 0
                if( read_fchk ):
                    NA, ELN, R, TOTAL_ENERGY, Q_ESP  = gaussian.parse_fchk( fchk_file )
                    run_qm = 1

                # Poppulate other atomic values                
                ASYMB = elements.eln_asymb(ELN)
                ATYPE = []
                CHARGES = []
                ELECTRONS_i = 0
                for atom_i in range(NA):
                    ATYPE.append( ASYMB[atom_i] )
                    CHARGES.append( 0.0 )
                    ELECTRONS_i += ELN[atom_i] 
                
                
                if( options.out_xyz ):
                    if( options.verbose ):
                        print "      Writing xyz file of fchk geometry ",options.out_xyz
                        xmol.write_xyz(ASYMB,R,options.out_xyz)


                # Run subsequent calculations if structure information was found 
                if( run_qm ):
                    os.chdir(struct_dir)
                    # Check for zmatrix file 
                    calc_id = job_name + "-ZMAT"
                    input_file =  calc_id+'.com'
                    calc_id_temp =  job_name+'_temp'
                    log_file =  calc_id+'.log'
                    zmat_fchk =  calc_id +"/" + calc_id +".fchk"
                    zmat_finished = file_io.file_exists( zmat_fchk )
                    # 
                    if( not zmat_finished):
			# Optimize z-matrix to get bonding information
			qm_kywd_o = options.qm_kywd 
			options.qm_kywd = qm_kywd_o + " OPT"
                        # Print com
                        gaussian.print_com( calc_id_temp, ASYMB,R,ATYPE,CHARGES,ELECTRONS_i,options)
                        #  geometry z-matrix opt input files
                        gaussian.com2zmat(calc_id_temp,calc_id,options)
                        
                        # Run optimization
                        if( options.submit ):
                            if( options.verbose ):
                                print "     Submitting Z-matrix optimization to queue "
                                print "         nodes        ",options.nnodes
                                print "         memory/node  ",options.pmem
                            # Print pbs script
                            pbs_id = cluster.write_pbs(pbs_templ,calc_id,input_file,options)
                            cluster.submit_job( struct_dir, pbs_id ,options )
                            
                        elif( options.localrun ):
                            if( options.verbose ):
                                print "       Running Z-matrix optimization  "
                            gaussian.run(options, calc_id)
                        else:
                            print " Please mark either qsub or runloc options as True to run qm"
                        			
			options.qm_kywd = qm_kywd_o

                # Check if finished
                #   basicaly if local run is done or should just skip over if submitted to queue
                zmat_finished = file_io.file_exists( zmat_fchk )
                # 
                if( zmat_finished):
                    if( options.verbose ):
                        print "       Parsing Z-matrix file to creat input files for run  "
                    RING_CONNECT, RING_NUMB, DIH_ID, DIH_VAL, DIH_ATOMS, zmatrix = pars_zmatrix(calc_id)
                        
                    # Check to see if list of dihedrals to loop over has been creaeted
                    dlist_name = job_name + "_dih.list"
                    dlist_exists = file_io.file_exists( dlist_name )
                    if ( dlist_exists ):
                        if( options.verbose ):
                            print "       Reading in dihedral list from   ",dlist_name
                        DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS = read_dihlist(dlist_name)
                        
                    else:
                        if( options.verbose ):
                            print "       Writing dihedral list ",dlist_name
                        DIH_TAG = tag_dih(RING_CONNECT, RING_NUMB,  DIH_ID, DIH_VAL, DIH_ATOMS)
                        write_dihlist(dlist_name, RING_NUMB, DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS )
    
                    dlist_exists = file_io.file_exists( dlist_name )
                    if ( dlist_exists ):                    
                        if( options.verbose ):
                            print "       Writing input files for all loop dihedrals "
                            print "         Nodes ",options.nnodes
                            print "         Processors  ",options.npros
			    
			qm_kywd_o = options.qm_kywd 
			options.qm_kywd = " popt=Zmat  nosym "
                        write_input(options, mol_dir,mol_id,mol_repeat,mol_acc, DIH_ID,DIH_TAG,DIH_VAL, DIH_ATOMS, zmatrix,pbs_templ,fix_templ,indx_file,work_dir)
			options.qm_kywd = qm_kywd_o
			
		
                os.chdir(work_dir)
        
    
if __name__=="__main__":
    main() 

