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
    

    parser.add_option("-v","--verbose", dest="verbose", default=True,action="store_true", help="Verbose output ")
    
    # Cluster options
    parser.add_option("--cluster_host", dest="cluster_host",type="string",default="peregrine",help=" name of cluster ")
    
    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    

    parser.set_defaults(ff_charges=False)
    parser.add_option("--ff_charges", dest="ff_charges",action="store_true",help=" Use ff charges ")
    
    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")

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
    import sys, os 
    import string, numpy, jsonapy , json 
    from string import replace
    import file_io, gaussian, elements, gromacs, lammps , top, atom_types, xmol 
    
    options, args = get_options()
    
    if( options.cluster_host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	
    elif( options.cluster_host == "redmesa" ):
	load_gaussian = "module load gaussian/g09/C.01"
        user = 'twkempe'
	

    elif( options.cluster_host == "macbook" ):
	load_gaussian = ""
        user = 'tkemper'
	
    # sufix for force field calcs to run multiple ff types or excultions, such as q(i) = 0.00 
    ff_type_id = "_fit"
    
    LAT_CONST = numpy.zeros( (3,3) )
    
    LAT_CONST[0][0] = 50.0
    LAT_CONST[1][1] = 50.0
    LAT_CONST[2][2] = 50.0

    # Store working dir  
    work_dir = os.getcwd()
    
    # Read index files from args
    for indx_file in args:
        # Get lines of record file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
        for gen_line in Lines:
            col = gen_line.split()
            if( len(col) >= 4 and col[0] == "gen" ):
		
                mol_dir = col[1].strip()
                mol_id = col[2].strip()
                mol_repeat = int(col[3].strip() )
                mol_acc = col[4].strip()
                
                # File info
                struct_dir = mol_dir + "/" + mol_id + "/"
                job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)

		# Read data from json file
		json_name = struct_dir +"/" + job_name +".json"
		json_atomicdata = 0
		json_data,json_success = jsonapy.read_jsondata(json_name)
		if(  json_success ):
		    if( options.verbose ):
			print "     Getting atomic data from  ",json_name
    
		    ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,success  = jsonapy.read_atomic(json_data)
		    if( success ):
			json_atomicdata = 1
		    else:
			 print "   json file ",json_name," exist, but does not contain any atomic data . "
		    
		else:
		    print "   json file ",json_name," does not exist. "
		#
		# CHeck for optimized geometry 
		#
                calc_id = job_name
		log_file = struct_dir +'/' + calc_id +"/"+calc_id+".log"
		fchk_file = struct_dir +'/' + calc_id +"/"+calc_id+".fchk"
		opt_finished = 0 
		if( file_io.file_exists( fchk_file) ):
                    if( options.verbose ):
			print  '    Get optimized geometry from ', fchk_file
			    
		    # Read in from zmatrix optimization 
		    NA, ELN, R, TOTAL_ENERGY, Q_ESP    = gaussian.parse_fchk( fchk_file )
		    
		    opt_finished = 1
		    
				
		else:
		    print "  could not find optimized geometry ",fchk_file

		
		# Get ESP chages 
		calc_id_esp = job_name + "-ESP"
		fchk_file_esp = struct_dir +'/' + calc_id_esp +"/"+calc_id_esp+".fchk"
		esp_finished = 0 
		if( file_io.file_exists( fchk_file_esp) ):
		    esp_finished = 1
		    if( options.verbose ):
			print "    Getting charges from esp fit ",fchk_file_esp
			print "        these charges will be used unless the ff_charges option is set to true "
		    NA, ELN, R, TOTAL_ENERGY , Q_ESP  = gaussian.parse_fchk( fchk_file_esp )
		    CHARGES = Q_ESP 
		    
		else:
		    print "  could not find esp calc ",fchk_file_esp
		    
		#
		# If atomic information is found make new files 
		#
		if( opt_finished or esp_finished or json_atomicdata ):
		    #
		    #   Build covalent nieghbor list for bonded information
		    #
		    NBLIST, NBINDEX = top.build_covnablist(ELN,R)
		    
		    if( not esp_finished ):
			if( options.verbose ):
			    print " Initialize charges to zero  since no ESP charges found "
			#
			# Initialize charges to zero 
			#
			CHARGES = top.initialize_charges( ELN )

		    #  Get atomic symbols and masses 
		    ASYMB = elements.eln_asymb(ELN)
		    AMASS = elements.eln_amass(ELN)
		    #   Initialize charge group #
		    CHARN = top.initialize_charn( ELN )

		    RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
		    RING_CONNECT  = top.find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB) 

		    # Asign oplsaa atom types
		    ATYPE, CHARGES = atom_types.oplsaa(  options.ff_charges,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB )
		    # Asign biaryl_types
		    ATYPE , CHARGES = atom_types.biaryl_types( options.ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
		    #Refind inter ring types
		    ATYPE , CHARGES  = atom_types.interring_types(options.ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )

		    #
		    # Set charge groups
		    #
		    CG_SET = []
		    one = 1
		    for i in range( len(ELN) ):
			CG_SET.append(one)			    
		    CHARN = top.set_chargegroups(options.verbose,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LAT_CONST)
		    
					    
		    if( not json_atomicdata ):
			
			if( options.verbose ):
			    print " Initialize ctype and "
			    
			#
			# Set unknowns to some default value 
			#
			CTYPE = []
			UNITNUMB = []
			UNITTYPE = []
			one = 1
			for i in range( len(ELN) ):
			    CTYPE.append("UNKNOWN") 
			    UNITNUMB.append(one)
			    UNITTYPE.append("UNKNOWN")
			    
			
		    # Write new atomic properties to json file

		    json_data = jsonapy.append_atomic(json_data,ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R)
		    
#                   Print new info files 
		    #
		    # Set cply_tag for cply output
		    #
		    cply_tag = top.set_cply_tags(  options.verbose, ELN, CTYPE,UNITNUMB ,NBLIST, NBINDEX )
		    # !!!! need to update CTYPES to show connections made !!!!
		    #
		    # Zero unit charges need to figure out though 
		    #
		    zero_term = 1
		    zero_func = 1 
		    CHARGES = top.zero_unitq(ELN,ATYPE,CHARGES,CTYPE,NBINDEX,NBLIST,options.verbose,zero_term,zero_func) 
		    #
		    # Print new building block with optimized geometry   and charges "qply" file type 
		    #
		    if( options.verbose): print "    Creating BuildingBlock_local directories "
		    bb_dir = struct_dir + "BuildingBlock_local"
		    if (not os.path.isdir(bb_dir)):
			os.mkdir(bb_dir)
			
		    for bb_subdir in ( "acceptors","donors","functional_groups","spacers","terminals" ):
			mk_dir  = bb_dir +"/"+bb_subdir
			if (not os.path.isdir(mk_dir)):
			    os.mkdir(mk_dir)
		    
		    
		    bb_dir_donor =  struct_dir + "BuildingBlock_local/donors"
		    bb_file = bb_dir_donor +"/" + mol_id + "n" + str(mol_repeat) + ".cply"
		    if( options.verbose): print "    Creating BuildingBlock_local file ",bb_file
		    F = open(bb_file,'w')
		    F.write('D()')
		    for atom_i in range( len(ELN) ):
			F.write('\n %s %16.6f %16.6f %16.6f %s ' % (ASYMB[atom_i],R[atom_i][0],R[atom_i][1],R[atom_i][2],str( cply_tag[atom_i]) ) )			
		    F.close()
		    #
		    # Prin qply file for building large oligomers for MD simulations
		    #
		    bb_file = bb_dir_donor +"/" + mol_id + "n" + str(mol_repeat) + ".cply"
		    if( options.verbose): print "    Creating BuildingBlock_local file ",bb_file
		    F = open(bb_file,'w')
		    F.write('D()')
		    for atom_i in range( len(ELN) ):
			#                                                              0           1            2             3              4              5                6                   7   
			F.write('\n %s %16.6f %16.6f %16.6f %16.6f %d %s %s ' % (ASYMB[atom_i],R[atom_i][0],R[atom_i][1],R[atom_i][2],CHARGES[atom_i],UNITNUMB[atom_i],UNITTYPE[atom_i],str( cply_tag[atom_i]) ) )			
		    F.close()
		    #
		    # Print new json with updated atomic charges 
		    #
		    # Need to merge with opv-project to get general json file write_meta
		    #frag.write_meta(json_data, xyz_name)
		    json_file = bb_dir_donor +"/" + mol_id + "n" + str(mol_repeat) + ".json"
		    
		    if( options.verbose): print "    Creating new json file ",json_file 
		    
		    f = open(json_file, 'w')
		    json.dump(json_data, f, indent=2)
		    f.close()
		    
		else:
		    print " No atomic information found for ",job_name
		    
		    
		    
		    
		
		

if __name__=="__main__":
    main() 

