#! /usr/bin/env python
# collect data from ab initio torsional potential run

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
    
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    
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
    import file_io, gaussian, elements, xmol
    
    options, args = get_options()

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
            if( len(col) >= 4 and col[0] == "qm_dih" ):
                
                mol_dir = col[1].strip()
                mol_id = col[2].strip()
                mol_repeat = int(col[3])
                mol_acc = col[4].strip()
                
                # File info
                struct_dir = mol_dir + "/" + mol_id + "/"
                job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)
                
                
                
		dih_id = col[5].strip()
		#a_k,a_i, a_j ,a_l,
		cent_min  = int(  col[10].strip() )
		cent_max = int(  col[11].strip() )
		cent_step = int(  col[12].strip() )
                
                
                print struct_dir , job_name , dih_id , cent_min , cent_max , cent_step

		print "   Running job ",struct_dir,job_name
		# Open output file
		
		dih_qm = struct_dir +'/' +job_name + '-' + dih_id + options.qm_sufix + '.dat'
		qm_out = open(dih_qm,'w')
		qm_out.write( '# cent_indx, cent_angle, qm_energy (eV) ,  all cojugated dih angles')


		xmol_qm = job_name+'-'+dih_id+"_qm.xmol"
		xmol_dir_qm = struct_dir +'/' + xmol_qm
		if( file_io.file_exists( xmol_dir_qm ) ): os.remove(xmol_dir_qm) 

		
		# Loop over angels of central dihedrals
		cent_indx = 0
		for cent_angle in range(cent_min,cent_max,cent_step):
		    calc_id =   job_name  + '-' + dih_id + '_' + str(cent_angle) + '_auxfix'
		    cent_indx += 1
	    
		    # Check to see if finished
		    fchk_file = struct_dir +'/' + calc_id +"/"+calc_id+".fchk"
		    
		    if( file_io.file_exists( fchk_file) ):
                        if( options.verbose ):
			    print  '    Get results from ', fchk_file
			    
			# Get energy 
			NA, ELN, R, TOTAL_ENERGY   = gaussian.parse_fchk( fchk_file )
			ASYMB = elements.eln_asymb(ELN)

			# Place 0.0 value as place holder
			DIH_ID  = []
			DIH_VAL = [] 
			DIH_ID.append( dih_id )
			dih_val_str = ''
		        for dih_indx in  range(len(DIH_ID)):
			    DIH_VAL.append( 0.0 )

			for indx_str in range( len(DIH_ID) ):
			    dih_val_str = dih_val_str + str( " %8.4f " % DIH_VAL[indx_str] )
			    
			debug = 0
			if( debug):
			    print ' printing qm energies '
			    print cent_indx,cent_angle, TOTAL_ENERGY,dih_val_str
			    sys.exit('qm print ')

			qm_out.write( " \n %8d %8.4f %16.8f  %50s" % ( cent_indx,cent_angle, TOTAL_ENERGY,dih_val_str))
			xmol.print_xmol(ASYMB,R,xmol_dir_qm)

	qm_out.close()
	


if __name__=="__main__":
    main() 

