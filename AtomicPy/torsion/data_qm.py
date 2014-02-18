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
    

    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")
    
    
    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")

    # How to run the needed calculations 
    parser.add_option("--submit", dest="submit",action="store_true", default=False,help=" submit calculations to the queue ")
    parser.add_option("--localrun", dest="localrun",action="store_true", default=False,help=" Run calculations locally")
    parser.add_option("--submit_command", dest="submit_command",type="string", default="qsub",help=" command used to submit script to the queue ")

    parser.add_option("--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")

    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    parser.add_option("--qm_load", dest="qm_load",type="string",help=" string to load qm software module  ")
    
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    
    (options, args) = parser.parse_args()

    # Set options based on cluster 
    if( options.host == "peregrine" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/.g09_C.01"
    elif( options.host == "redmesa" ):
        if( options.qm_software == "gaussian" ):
            options.qm_load = "module load gaussian/g09/C.01"
	    
	    
    return options, args

def main():
    import sys, os , string
    import jsonapy
    import file_io, gaussian, elements, xmol
    
    options, args = get_options()

    # Store working dir  
    work_dir = os.getcwd()
    
    
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
		dih_id_list ,cent_min_list ,cent_max_list ,cent_step_list,a_k_list, a_i_list, a_j_list, a_l_list,qmtor_found = jsonapy.read_qm_tor(json_data)
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found and qmtor_found  ):
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
				NA, ELN, R, TOTAL_ENERGY , Q_ESP  = gaussian.parse_fchk( fchk_file )
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

