#! /usr/bin/env python
"""
Read and write json files 
"""

# Dr. Travis Kemper
# NREL
# 12/09/2013
# travis.kemper@nrel.gov

def read_jsondata(json_name):
    """
    Read in json file 
    """
    import file_io
    import json 
    
    success = 0
    
    if(  file_io.file_exists( json_name) ):
	
	f = open(json_name, 'r')
	json_data = json.load(f)
	f.close()
	
	success = 1
	
    else:
	# initialize
	json_data = {}
    
    return json_data,success


def read_meta(json_data):
    """
    Return meta data from OPV database json file 
    """
    import json , numpy 
    
    success = 0

    #print json_data['metadata']["atomic
    metadata_found,atomicdata_found = check_atomic(json_data)
    
	
    if (  metadata_found  ):
	print " Reading meta data from json file "
	mol_dir = json_data['metadata']['mol_dir']
	accuracy = json_data['metadata']['accuracy']
	donor_substituents = json_data['metadata']['donor_substituents']
	donors = json_data['metadata']['donors']
	acceptor_substituents = json_data['metadata']['acceptor_substituents']
	basis = json_data['metadata']['basis']
	method = json_data['metadata']['method']
	
	terminal_substituents = json_data['metadata']['terminal_substituents']
	number = json_data['metadata']['number']
	spacers = json_data['metadata']['spacers']
	terminals = json_data['metadata']['terminals']
	tag = json_data['metadata']['tag']
	n_units = json_data['metadata']['n']
	spacer_substituents = json_data['metadata']['spacer_substituents']
	acceptors = json_data['metadata']['acceptors']

	success = 1
	    
    return (mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,success)


def check_atomic(json_data):
    """
    Test if atomic data exisits in  OPV database json file 
    """
    
    metadata_found = 0 
    atomicdata_found = 0
    
    #
    # Check for metadata section 
    #
    for data in json_data:
	if( data == 'metadata' ):
	    metadata_found = 1
	
    if( metadata_found ):
	for meta_data in json_data['metadata']:
	    if ( meta_data == "atomic" ):
		atomicdata_found = 1
		
    return (metadata_found,atomicdata_found)


def read_atomic(json_data):
    """
    Read atomic data from a  OPV database json file 
    """
    import json , numpy 
    
    success = 0

    ELN = []
    ASYMB = []
    CTYPE = []
    CHARGES = []
    R = []
    UNITNUMB = [] 
    UNITTYPE = []
    
    #print json_data['metadata']["atomic
    metadata_found,atomicdata_found = check_atomic(json_data)
    
	
    if (  atomicdata_found  ):
	print " Reading atomic data from json file with ",json_data['metadata']["atomic"]["natoms"]," atoms "
	
	natoms = json_data['metadata']["atomic"]["natoms"]  
	print "       Number of atoms found in json file ",natoms
	for atom_i in range( natoms ):
	    ELN.append( json_data['metadata']["atomic"]["element"][atom_i]  )
	    ASYMB.append( json_data['metadata']["atomic"]["asymb"][atom_i]  )
	    CTYPE.append( json_data['metadata']["atomic"]["ctype"][atom_i]  )
	    CHARGES.append( json_data['metadata']["atomic"]["q"][atom_i]  )
	    UNITNUMB.append( json_data['metadata']["atomic"]["unitnumb"][atom_i]  )
	    UNITTYPE.append( json_data['metadata']["atomic"]["unittype"][atom_i]  )
	    pos_str = json_data['metadata']["atomic"]["pos"][atom_i].split()
	    r_i = numpy.array( [float(pos_str[0]),float(pos_str[1]),float(pos_str[2])] )
	    R.append( r_i )
	    
	    success = 1
	    
    return (ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,success)

def append_atomic(json_data,ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R):
    """
    Append atomic data to a  OPV database json file 
    """
    import file_io
    
    
    metadata_found,atomicdata_found = check_atomic(json_data)
    
    if( not metadata_found ):
	# Add metadata section if not in json_data 
	meta_data = {}
	json_data["metadata"] = meta_data 
	
	
    if( not atomicdata_found ):
	# Add metadata section if not in json_data 
	atomic_data = {}
        json_data['metadata']["atomic"] = atomic_data
	
        atomic_data["natoms"] = len(ELN) 
        atomic_data["element"] = []
        atomic_data["asymb"] = []
        atomic_data["q"] = []
        atomic_data["ctype"] = []    
        atomic_data["unitnumb"] = []
        atomic_data["pos"] = []
        atomic_data["unittype"] = []
	
	
        for atom_i in range( len(ELN) ):
            atomic_data["element"].append( ELN[atom_i] )
            atomic_data["asymb"].append( ASYMB[atom_i] )
            atomic_data["ctype"].append( CTYPE[atom_i] )
            atomic_data["q"].append(CHARGES[atom_i])
            atomic_data["unitnumb"].append( UNITNUMB[atom_i] )
            pos_str = str( "%f %f %f "% (R[atom_i][0],R[atom_i][1],R[atom_i][2]) )
            atomic_data["pos"].append( pos_str )
            atomic_data["unittype"].append( UNITTYPE[atom_i] )
	    
    else:
            
	
        for atom_i in range( len(ELN) ):
            json_data['metadata']["atomic"]["element"][atom_i]  = ELN[atom_i] 
            json_data['metadata']["atomic"]["asymb"][atom_i]  = ASYMB[atom_i] 
            json_data['metadata']["atomic"]["ctype"][atom_i]  =CTYPE[atom_i] 
            json_data['metadata']["atomic"]["q"][atom_i]  = CHARGES[atom_i]
            json_data['metadata']["atomic"]["unitnumb"][atom_i]  = UNITNUMB[atom_i] 
            pos_str = str( "%f %f %f "% (R[atom_i][0],R[atom_i][1],R[atom_i][2] ) )
            json_data['metadata']["atomic"]["pos"][atom_i]  = pos_str 
            json_data['metadata']["atomic"]["unittype"][atom_i]  = UNITTYPE[atom_i] 
    
    return json_data


def read_qm_tor(json_data):
    """
    Read torisonal potential calculation data from  OPV database json file 
    """
    dih_id_list = []
    cent_min_list  = []
    cent_max_list = []
    cent_step_list = []
    a_k_list = []
    a_i_list = [] 
    a_j_list = []
    a_l_list = []
    
    qmtor_found = 0
    
    #
    # Check for metadata section 
    #
    for data in json_data:
	if( data == 'metadata' ):
	    
	    for meta_data in json_data['metadata']:
	        if ( meta_data == "qm_tor_data" ):
		    qmtor_found = 1
		    
		    dih_id_list = json_data['metadata']["qm_tor_data"]["cent_id"]
		    cent_min_list = json_data['metadata']["qm_tor_data"]["cent_min"]
		    cent_max_list = json_data['metadata']["qm_tor_data"]["cent_max"]
		    cent_step_list = json_data['metadata']["qm_tor_data"]["cent_step"]
		    a_k_list = json_data['metadata']["qm_tor_data"]["a_k"]
		    a_i_list = json_data['metadata']["qm_tor_data"]["a_i"]
		    a_j_list = json_data['metadata']["qm_tor_data"]["a_j"]
		    a_l_list = json_data['metadata']["qm_tor_data"]["a_l"]
		    
		
    return ( dih_id_list ,cent_min_list ,cent_max_list ,cent_step_list,a_k_list, a_i_list, a_j_list, a_l_list, qmtor_found )


def read_ff_tor(json_data):
    """
    Read force-field torisonal potential calculation data from  OPV database json file 
    """
    dih_id_list = []
    cent_min_list  = []
    cent_max_list = []
    cent_step_list = []

    a_k_list = []
    a_i_list = []
    a_j_list = []
    a_l_list = []

    ff_type_list = []
    
    success  = 0
    
    #
    # Check for metadata section 
    #
    for data in json_data:
	if( data == 'metadata' ):
	    
	    for meta_data in json_data['metadata']:
	        if ( meta_data == "ff_tor_data" ):
		    success = 1
		    
		    dih_id_list = json_data['metadata']["ff_tor_data"]["cent_id"]
		    cent_min_list = json_data['metadata']["ff_tor_data"]["cent_min"]
		    cent_max_list = json_data['metadata']["ff_tor_data"]["cent_max"]
		    cent_step_list = json_data['metadata']["ff_tor_data"]["cent_step"]
		    a_k_list = json_data['metadata']["ff_tor_data"]["a_k"]
		    a_i_list = json_data['metadata']["ff_tor_data"]["a_i"]
		    a_j_list = json_data['metadata']["ff_tor_data"]["a_j"]
		    a_l_list = json_data['metadata']["ff_tor_data"]["a_l"]
		    ff_type_list  = json_data['metadata']["ff_tor_data"]["ff_type"]
		    
		
    return ( dih_id_list ,cent_min_list ,cent_max_list ,cent_step_list,a_k_list, a_i_list, a_j_list, a_l_list, ff_type_list, success )
