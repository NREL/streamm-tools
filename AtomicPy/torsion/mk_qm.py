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
    
    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")
    
    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")

    # How to run the needed calculations 
    parser.add_option("--submit", dest="submit",action="store_true", default=False,help=" submit calculations to the queue ")
    parser.add_option("--localrun", dest="localrun",action="store_true", default=False,help=" Run calculations locally")
    parser.add_option("--submit_command", dest="submit_command",type="string", default="qsub",help=" command used to submit script to the queue ")

    parser.add_option("--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")

    # Torsion
    parser.add_option("--cent_min", dest="cent_min", type="int", default="0",help=" Initial torsional angle ")
    parser.add_option("--cent_max", dest="cent_max", type="int", default="180",help=" Final torsional angle ")
    parser.add_option("--cent_step", dest="cent_step", type="int", default="5",help=" Step size torsional angle ")

    parser.add_option("--dih_temp", dest="dih_temp", type="string", default="mp2_dih.com.template",help=" Template for Links of dihedral calculation ")

    # QM calculation options 
    parser.set_defaults(qm_software="gaussian")
    parser.add_option("--qm_software", dest="qm_software",type="string",help=" what software to use for the qm calculations   ")
    parser.add_option("--qm_load", dest="qm_load",type="string",help=" string to load qm software module  ")

    parser.add_option("--qm_kywd", dest="qm_kywd", type="string",default="", help="Key words for QM calculation ")
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    parser.add_option("--qm_charge", type="int",action="append", default="0",help="Input gaussain log file ")
    parser.add_option("--qm_mult", dest="qm_mult", type="int",default="0", help=" Shift in default spin multiplicity ( singlet,doublet) QM calculation, allows for triplets ")

    parser.add_option("--qm_method", dest="qm_method", type="string",default="B3LYP", help="Method of QM calculation ")
    parser.add_option("--qm_basis", dest="qm_basis", type="string",default="6-31G**", help="Basis set of QM calculation ")

    parser.add_option("--sp_meth", dest="sp_meth", type="string", default="MP2",help=" Method for set for hihgh level single point energy calculations ")
    parser.add_option("--sp_basis", dest="sp_basis", type="string", default="cc-pVTZ",help=" Basis for set for hihgh level single point energy calculations")

    parser.add_option("--pmem",dest="pmem", type="int", default="1700",help=" Memory per processor ")
    parser.add_option("--npros", dest="npros", type="int", default="4",help=" Number of processors ")
    parser.add_option("--nnodes", dest="nnodes", type="int", default="1",help=" Number of nodes ")

    # Output options 
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help="Output single frame xyz file in xmol format ")

    (options, args) = parser.parse_args()
	    
    return options, args

def read_dihlist(dlist_name):
    
    
    DIH_ID = []
    DIH_VAL = []
    DIH_ATOMS = []
    DIH_TAG = []
    
    # get lines 
    f = open(dlist_name,'r')
    Lines = f.readlines()
    f.close()
    
    
    for line in Lines: 
        col = line.split()
        if( col[0] == "dih" and len(col) >= 11 ):
            dih_indx = col[1]
            DIH_TAG.append( col[2] )
            DIH_ID.append( col[3] )
            DIH_VAL.append(  col[4] )
            a_l = int( col[5] ) - 1
            a_i = int( col[6] ) - 1
            a_j = int( col[7] ) - 1
            a_k = int( col[8] ) - 1
            DIH_ATOMS.append( [ a_l,a_i,a_j,a_k ] )
        
  
    return (  DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS )
    
def write_dihlist(dlist_name, RING_NUMB, DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS ):
    
    flist = open(dlist_name,"w")
    
    for dih_indx in  range(len(DIH_ID)):
        a_l = DIH_ATOMS[dih_indx][0]
        a_i = DIH_ATOMS[dih_indx][1]
        a_j = DIH_ATOMS[dih_indx][2]
        a_k = DIH_ATOMS[dih_indx][3]
        
        flist.write( " dih %d %s %s %f %d %d %d %d %d %d %d %d \n"  % (dih_indx,DIH_TAG[dih_indx],DIH_ID[dih_indx], DIH_VAL[dih_indx],a_l+1,a_i+1,a_j+1,a_k+1 ,RING_NUMB[a_l],RING_NUMB[a_i],RING_NUMB[a_j],RING_NUMB[a_k] ) )
        
    flist.close()

def print_dih_com( options,job_name,struct_dir,fix_templ,cent_indx, cent_angle ,DIH_ID,DIH_TAG,DIH_VAL, zmatrix ):
    from string import replace

    cent_name =   job_name  + '-' + DIH_ID[cent_indx].strip() + '_' + str(cent_angle) + '_auxfix'
    cent_tor_cord = 'tor_cord '+str(cent_angle) + '_auxfix'
    calc_id =  cent_name
    cent_id = DIH_ID[cent_indx]
       
    input_file =  calc_id + '.com'
    com_name = struct_dir +'/'+ input_file
    am_opt = '%chk='+str(calc_id)+'.chk' 
    am_opt = am_opt + '\n%nproc='+str(options.npros)
    #am_opt = am_opt + '\n' + '#P AM1/3-21g  popt=Zmat nosym '
    keywd = str(options.qm_method) +"/"+ str(options.qm_basis) +" "+ str(options.qm_kywd)
    am_opt = am_opt + "\n #P " + keywd
    #f.write( "\n# P %s/%s  %s" % (opt:qions.qm_method,options.qm_basis,options.qm_kywd))
    am_opt = am_opt + '\n' + ' '
    am_opt = am_opt + '\n' + 'tor_cord  ' + str(cent_angle) #str(job_name)
    am_opt = am_opt + '\n' 
    am_opt = am_opt + '\n 0 1 '
    
    # Prin all non constranied elments of zmatrix 
    for line in iter(zmatrix.splitlines()) :
        print_z = 1
	colvar = line.split('=')
        var_id = colvar[0].strip()
        for dih_indx in  range(len(DIH_ID)):
            if( DIH_ID[dih_indx] == var_id ):
                if(  DIH_TAG[dih_indx] != "relax"  ):
                    print_z = 0
                    break
            
        if( print_z ):
            am_opt =  am_opt+"\n" + line 
            
            
    am_opt = am_opt + '\n Constants: '
    # Print dihedral id's and angles as fixed
    for  dih_indx in  range(len(DIH_ID)):
	if (  DIH_ID[dih_indx] !=  DIH_ID[cent_indx]):
            if(  DIH_TAG[dih_indx] != "relax"  ):
                am_opt = am_opt +"\n" + str(DIH_ID[dih_indx]) + " " + str(DIH_VAL[dih_indx]) + '  F  '
	else:
	    am_opt = am_opt +"\n" + str(DIH_ID[dih_indx]) + " " + str(cent_angle) + '.0  F  '
	    
		
    am_opt = am_opt + ' \n'
    am_opt = am_opt + ' \n'
    temp_fix = fix_templ 
    temp_fix = replace(temp_fix,"<calc_id>",calc_id)
    temp_fix = replace(temp_fix,"<structure_name>",job_name)
    temp_fix = replace(temp_fix,"<dih_id>",cent_id)
    temp_fix = replace(temp_fix,"<dih_angle>",str(cent_angle)+'.0')
    temp_fix = replace(temp_fix,"<sp_meth>",options.sp_meth)
    temp_fix = replace(temp_fix,"<sp_basis>",options.sp_basis)
    temp_fix = replace(temp_fix,"<npros>",str( options.npros) )
    
    #temp_fix = replace(temp_fix,"<qm_method>",str( options.qm_method) )
    #temp_fix = replace(temp_fix,"<qm_basis>",str( options.qm_basis) )
    
    com_dih = am_opt + temp_fix
    
    return com_dih
	    

def build_nablist(ELN,BONDS):
    import sys,numpy

    debug = 0
    NNAB = 0

    maxnnab = len(BONDS)*2 + 1
    NBLIST = numpy.empty( maxnnab,  dtype=int )
    NBINDEX = numpy.empty( maxnnab,  dtype=int )

    if(debug ):
	for b in range( len(BONDS) ) :
	    bnd_i = BONDS[b][0]
	    bnd_j = BONDS[b][1]
	    print b, bnd_i,bnd_j
	    
    for i in range(len(ELN)):
	NBINDEX[i] = NNAB + 1
	for b in range( len(BONDS) ) :
	    bnd_i = BONDS[b][0]
	    bnd_j = BONDS[b][1]
	    if ( i == bnd_i ):
		NNAB = NNAB + 1
		NBLIST[NNAB] =  bnd_j
		if(debug): print " adding bnd_i",i,b,bnd_i,bnd_j,len(NBLIST),NBINDEX[i]
	
	    if ( i == bnd_j ):
		NNAB = NNAB + 1
		NBLIST[NNAB] =  bnd_i
		if(debug): print " adding bnd_j",i,b,bnd_i,bnd_j,len(NBLIST),NBINDEX[i]
		    
    #if(debug): sys.exit('debug')
    
    # Account for final atom position
    NBINDEX[i+1] =  NNAB + 1

    debug = 0
    if ( debug ):
       print ' total nbs ',NNAB,' total bonds ',len(BONDS)
       for i in range(len(ELN)):
	    N_o = NBINDEX[ i  ]
	    N_f = NBINDEX[ i + 1 ] - 1
	    NNAB = N_f - N_o + 1
	    print ' atom ', i,' ',ELN[i],NNAB,N_o
    
	    #
	    # Find number of elements
	    #
	    #ELCNT = numpy.zeros(120)
	    for indx in range( N_o,N_f+1):
		j = NBLIST[indx]
		print ELN[j],j
		#    el_j = ELN[j]
		#    ELCNT[j] = ELCNT[j] + 1

       sys.exit('debug')


    return (NBLIST,NBINDEX)

def pars_zmatrix( calc_id, job_name ):
    from string import replace
    # atomicpy
    import gaussian, top
    
    log_name = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.log" )
    com_name = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.com" )
    # Get lines of log file     
    f = open(log_name,'r')
    Lines = f.readlines()
    f.close()

    #Parse fchk
    fchk_file = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.fchk" )

    print " fchk_file" , fchk_file
    NA, ELN, R, TOTAL_ENERGY, Q_ESP  = gaussian.parse_fchk( fchk_file )
    
    # Parse log file 
    BONDS   , ANGLES, DIH = gaussian.read_optlog(log_name)
    NBLIST,NBINDEX = build_nablist(ELN,BONDS)
    
    # Find rings 
    RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
    RING_CONNECT  = top.find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB) 

    zmatrix = gaussian.com_zmatrix(com_name)    
    
    DIH_ID, DIH_VAL, DIH_ATOMS = gaussian.get_dih_id( zmatrix)
    
    return ( RING_CONNECT, RING_NUMB, DIH_ID, DIH_VAL, DIH_ATOMS, zmatrix )


def tag_dih(options,RING_CONNECT, RING_NUMB, DIH_ID, DIH_VAL, DIH_ATOMS):
    
    DIH_TAG = []
    
    # Initialize all dihedrals to relax 
    for dih_indx in  range(len(DIH_ID)):
        DIH_TAG.append("relax")

    
    for dih_indx in  range(len(DIH_ID)):
    
        zmat_l = DIH_ATOMS[dih_indx][0]
        zmat_i = DIH_ATOMS[dih_indx][1]
        zmat_j = DIH_ATOMS[dih_indx][2]
        zmat_k = DIH_ATOMS[dih_indx][3]
        
        # Make sure inter ring connect not improper 
        if (  RING_NUMB[zmat_i] != RING_NUMB[zmat_j] and RING_NUMB[zmat_l] != RING_NUMB[zmat_k] ):
            DIH_TAG[dih_indx] = "fix"
	    
	    if ( options.verbose ):
		print "  Fixing dih ",DIH_ID[dih_indx], DIH_VAL[dih_indx], DIH_ATOMS[dih_indx]
	    
            if( RING_NUMB[zmat_i] > 0 and  RING_NUMB[zmat_j] > 0   ):
                DIH_TAG[dih_indx] = "loop"
    
		if ( options.verbose ):
		    print "  Looping dih ",DIH_ID[dih_indx], DIH_VAL[dih_indx], DIH_ATOMS[dih_indx]
		
                        
                        
    return DIH_TAG

def write_dihlist(dlist_name, RING_NUMB, DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS ):
    
    flist = open(dlist_name,"w")
    
    for dih_indx in  range(len(DIH_ID)):
        a_l = DIH_ATOMS[dih_indx][0]
        a_i = DIH_ATOMS[dih_indx][1]
        a_j = DIH_ATOMS[dih_indx][2]
        a_k = DIH_ATOMS[dih_indx][3]
        
        flist.write( " dih %d %s %s %f %d %d %d %d %d %d %d %d \n"  % (dih_indx,DIH_TAG[dih_indx],DIH_ID[dih_indx], DIH_VAL[dih_indx],a_l+1,a_i+1,a_j+1,a_k+1 ,RING_NUMB[a_l],RING_NUMB[a_i],RING_NUMB[a_j],RING_NUMB[a_k] ) )
        
    flist.close()

def write_input( options,  json_data, struct_dir ,job_name , DIH_ID,DIH_TAG,DIH_VAL, DIH_ATOMS, zmatrix,pbs_templ,fix_templ ,work_dir):
    import cluster, json 
    
    json_file_loc = "%s%s" % ( job_name , ".json" )

    
    #append_qm_tor_json
    qm_tor_data = {}
    json_data['metadata']["qm_tor_data"] = qm_tor_data
    

    qm_tor_data["cent_id"] = []
    qm_tor_data["a_k"] =  []
    qm_tor_data["a_i"] =  []
    qm_tor_data["a_j"] =  []
    qm_tor_data["a_l"] =  []
    qm_tor_data["cent_min"] =  []
    qm_tor_data["cent_max"] =  []
    qm_tor_data["cent_step"] = []
			
			    
    # Loop over central dihedrals 
    for cent_indx in range(len(DIH_ID)):
        print DIH_TAG[cent_indx].strip() 
        if( DIH_TAG[cent_indx].strip()  == "loop" ):
                
            a_l = DIH_ATOMS[cent_indx][0]
            a_i = DIH_ATOMS[cent_indx][1]
            a_j = DIH_ATOMS[cent_indx][2]
            a_k = DIH_ATOMS[cent_indx][3]
                
        
            cent_id = DIH_ID[cent_indx].strip()
            print " The connection between ", a_k,a_i, a_j ,a_l, " is ", cent_id,  DIH_ATOMS[cent_indx]
    
	    #
	    # Append torsional information 
	    #
	    
	    qm_tor_data["cent_id"].append( cent_id )
	    qm_tor_data["a_k"].append( a_k )
	    qm_tor_data["a_i"].append( a_i )
	    qm_tor_data["a_j"].append( a_j )
	    qm_tor_data["a_l"].append( a_l )
	    qm_tor_data["cent_min"].append( options.cent_min )
	    qm_tor_data["cent_max"].append( options.cent_max+options.cent_step )
	    qm_tor_data["cent_step"].append( options.cent_step )
			    
				    
            # Loop over angels of central dihedrals
            for cent_angle in range(options.cent_min,options.cent_max+options.cent_step,options.cent_step):
                cent_name =   job_name  + '-' + cent_id + '_' + str(cent_angle) + '_auxfix'
                cent_tor_cord = 'tor_cord '+str(cent_angle) + '_auxfix'
                calc_id =  cent_name
                com_dih = print_dih_com(options,job_name,struct_dir,fix_templ,cent_indx, cent_angle ,DIH_ID,DIH_TAG,DIH_VAL, zmatrix )
                        
                input_file =  calc_id + '.com'
                com_name = input_file
    
                f = file(com_name, "w")
                f.write(com_dih)
                f.close()
                
		if( options.host == "peregrine" ):
		    pbs_id = cluster.write_pbs(pbs_templ,calc_id,input_file,options)


    f = open(json_file_loc	, 'w')
    json.dump(json_data, f, indent=2)
    f.close()
    
    
def main():
    import string, os , sys 
    # atomicpy
    import jsonapy 
    import gaussian, elements, xmol , file_io , cluster 
    from string import replace
    
    options, args = get_options()
    
    # Verbose output
    if( options.verbose ):
        print "   Calculations (sufix): "
        print "     Z-matrix optimization with ESP charges  (-ZMATOPT) "
        if( options.submit):
            print "  Calculations "
        print "  Please load module ",options.qm_load

    # Read in template files 
    f = open(options.dih_temp,'r')
    fix_templ = f.read()
    f.close()

    # 
    if( options.host == "peregrine" ):
	f = open("peregrine.pbs.template",'r') # redmesa.slurm.template
	pbs_templ = f.read()
	f.close()
	    
	# Verbose output
	if( options.verbose ):   print "  Please load module ",options.qm_load
	
    elif( options.host == "macbook" ):
	pbs_templ= ""


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
		    calc_id = "%s/%s%s" % (struct_dir, job_name , "-ZMATOPT" )
		    
		    zmat_fchk = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.fchk" )
		    print " Checking for complete zmatrix optimiztion ",zmat_fchk
		    
		    zmat_finished = 0
		    if( gaussian.check_fchk( zmat_fchk ) == 0 ): zmat_finished =1 
                    #zmat_finished =  file_io.file_exists( zmat_fchk )
		    
		    if( zmat_finished  ):

			NA, ELN, R, TOTAL_ENERGY, Q_ESP  = gaussian.parse_fchk( zmat_fchk )
			RING_CONNECT, RING_NUMB, DIH_ID, DIH_VAL, DIH_ATOMS, zmatrix = pars_zmatrix(calc_id,job_name)
	
	
			# Poppulate other atomic values                
			ASYMB = elements.eln_asymb(ELN)
			ATYPE = []
			ELECTRONS_i = 0
			for atom_i in range(NA):
			    ATYPE.append( ASYMB[atom_i] )
			    ELECTRONS_i += ELN[atom_i] 
				
				
			
			if( options.out_xyz ):
			    if( options.verbose ):
				print "      Writing xyz file of fchk geometry ",options.out_xyz
			    xmol.write_xyz(ASYMB,R,options.out_xyz)

			#
			# Updat atomic information
			#
			#json_data = jsonapy.append_atomic(json_data,ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R)

			if( options.verbose ):
			    print "       Parsing Z-matrix file to creat input files for run  "
			    
			
	                os.chdir(struct_dir)
			    
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
			    DIH_TAG = tag_dih(options,RING_CONNECT, RING_NUMB,  DIH_ID, DIH_VAL, DIH_ATOMS)
			    write_dihlist(dlist_name, RING_NUMB, DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS )
			    
			dlist_exists = file_io.file_exists( dlist_name )
			if ( dlist_exists ):                    
			    if( options.verbose ):
				print "       Writing input files for all loop dihedrals "
				print "         Nodes ",options.nnodes
				print "         Processors  ",options.npros
				
			    qm_kywd_o = options.qm_kywd 
			    options.qm_kywd = " popt=Zmat  nosym "
			    

			    write_input(options, json_data, struct_dir ,job_name , DIH_ID,DIH_TAG,DIH_VAL, DIH_ATOMS, zmatrix,pbs_templ,fix_templ,work_dir)
			    options.qm_kywd = qm_kywd_o
			    
	            else:
                    	print " zmatrix file not optimized ",zmat_fchk
	
                os.chdir(work_dir)
        
    
if __name__=="__main__":
    main() 

