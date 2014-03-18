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
        

    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    parser.add_option("-r","--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")

    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")

    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    

    parser.set_defaults(ff_charges=False)
    parser.add_option("--ff_charges", dest="ff_charges",action="store_true",help=" Use ff charges ")
    
    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")

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
    import sys, os 
    import string, numpy, jsonapy , json 
    from string import replace
    import file_io, gaussian, elements, gromacs, lammps , top, atom_types, xmol 
    
    options, args = get_options()
    
    if( options.host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
	
    elif( options.host == "dale" ):
	load_gaussian = "module load gaussian" #/g09/C.01"
        user = 'tkemper'	

    elif( options.host == "macbook" ):
	load_gaussian = ""
        user = 'tkemper'
	
    # sufix for force field calcs to run multiple ff types or excultions, such as q(i) = 0.00 
    ff_type_id = "_fit"
    
    LV = numpy.zeros( (3,3) )
    
    LV[0][0] = 50.0
    LV[1][1] = 50.0
    LV[2][2] = 50.0

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
			
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found ):
		    if( options.verbose ):
			print " Meta data found will use specified method and basis unless others are specified in the options "
		    #
		    # Construct file names 
		    #
		    #short_name = "acc%d_%s_n%d" % (accuracy, tag, number )
		    job_name = "acc%d_%s_n%d" % (accuracy, tag, n_units )
		    struct_dir = "%s/%s/" % (mol_dir, tag )
		    calc_id = "%s%s" % ( job_name , "-ESP" )
		    
		    # json_atomicdata = 0
		    if( options.verbose ):
			print "     Getting atomic data from  ",json_file
    
		    ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,success  = jsonapy.read_atomic(json_data)
		    
		    if( success ):
			fchk_file = struct_dir +'/' + calc_id +"/"+calc_id+".fchk"
			
			if( file_io.file_exists( fchk_file) ):
			    if( options.verbose ):
				print  '    Get optimized geometry from ', fchk_file
			    
			    # Read in from zmatrix optimization 
			    NA, ELN, R, TOTAL_ENERGY, Q_ESP    = gaussian.parse_fchk( fchk_file )
			    
			    #  Get atomic symbols and masses 
			    ASYMB = elements.eln_asymb(ELN)
			    AMASS = elements.eln_amass(ELN)
			    
			    CHARGES = Q_ESP 
		    
			    #
			    #   Build covalent nieghbor list for bonded information
			    #
			    NBLIST, NBINDEX = top.build_covnablist(ELN,R)

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
				
			    CHARN = top.set_chargegroups(options.verbose,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LV)
			    
			    #
			    #    
			    #    if( not json_atomicdata ):
			    #	
			    #	if( options.verbose ):
			    #	    print " Initialize ctype and "
			    #	#
			    #	# Set unknowns to some default value 
			    #	#
			    #	CTYPE = []
			    #	UNITNUMB = []
			    #	UNITTYPE = []
			    #	one = 1
			    #	for i in range( len(ELN) ):
			    #	    CTYPE.append("UNKNOWN") 
			    #	    UNITNUMB.append(one)
			    #	    UNITTYPE.append("UNKNOWN")
			    #	    
			    #
			    #
			    # Write new atomic properties to json file
	
			    json_data = jsonapy.append_atomic(json_data,ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R)
				
			    out_xyz = struct_dir +'/' + calc_id + '.xyz'
			    
			    if( options.verbose ):
				print "  Writing xyz file for new atomic order reference ",out_xyz
			    xmol.write_xyz(ASYMB,R,out_xyz)

			    debug = 0 
			    if( debug ):
				for i in range( len(ELN) ):
				    print " ",i+1,ELN[i]
						      
				sys.exit("  atomic order issues ")
			    
			    #
			    # Print new info files 
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
			    bb_file = bb_dir_donor +"/" + tag + "n" + str(n_units) + ".cply"
			    if( options.verbose): print "    Creating BuildingBlock_local file ",bb_file
			    F = open(bb_file,'w')
			    F.write('D()')
			    for atom_i in range( len(ELN) ):
				F.write('\n %s %16.6f %16.6f %16.6f %s ' % (ASYMB[atom_i],R[atom_i][0],R[atom_i][1],R[atom_i][2],str( cply_tag[atom_i]) ) )			
			    F.close()
			    #
			    # Prin qply file for building large oligomers for MD simulations
			    #
			    bb_file = bb_dir_donor +"/" + tag + "n" + str(n_units) + ".cply"
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
			    json_file = bb_dir_donor +"/" + tag + "n" + str(n_units) + ".json"
			    
			    if( options.verbose): print "    Creating new json file ",json_file 
			    
			    f = open(json_file, 'w')
			    json.dump(json_data, f, indent=2)
			    f.close()
	
				
			else:
			    print "  could not find ZMATOPT calc ",fchk_file
			
			
		    else:
			 print "   json file ",json_file," exist, but does not contain any atomic data . "
		    
		else:
		    print "   no meta data in ",json_file,"  "
		 
		 
		    
		
		

if __name__=="__main__":
    main() 

