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


def read_connections(json_data):
    """
    Read atomic data from a  OPV database json file 
    """
    import json , numpy  , sys 

    debug = 0
    success = 0

    BONDS = []

    atomicdata_found = 0
    connectionsdata_found = 0
    
    #
    # Check for metadata section 
    #
    for data in json_data:
	if( data == 'metadata' ):
	    metadata_found = 1

    if( metadata_found ):
	for meta_data in json_data['metadata']:
	    if ( meta_data == "connections" ):
		connectionsdata_found = 1

    if( connectionsdata_found ):
        
	nbonds = json_data['metadata']["connections"]["nbonds"]  
	for bond_i in range( nbonds ):
            bond_col =  json_data['metadata']["connections"]["bonds"][bond_i].split()
            BONDS.append( [int(bond_col[0]),int(bond_col[1]) ] )
            if(debug):
                print bond_i
                print bond_col
                print BONDS[bond_i][0],BONDS[bond_i][1]
                print ""
            
    if(debug):
        sys.exit("read_connections debug 1 ")
        
    return BONDS
    
def read_atomic(json_data):
    """
    Read atomic data from a  OPV database json file 
    """
    import json , numpy, sys 
    #from particles import Particle
    #from particles import ParticleContainer

    #from bonds import Bond
    #from bonds import BondContainer

    #from structure import Structure

    verbose = 1 
    
    success = 0

    ELN = []
    ASYMB = []
    CTYPE = []
    CHARGES = []
    R = []
    VEL = []
    ATYPE = []
    AMASS = []
    MOLNUMB = []
    RING_NUMB = []
    RESID = []
    RESN = []
    CHARN = []

    UNITNUMB = [] 
    UNITTYPE = []
    
    #print json_data['metadata']["atomic
    metadata_found,atomicdata_found = check_atomic(json_data)
    
	
    if (  atomicdata_found  ):
	print " Reading atomic data from json file with ",json_data['metadata']["atomic"]["natoms"]," atoms "
	
	natoms = json_data['metadata']["atomic"]["natoms"]  
	print "       Number of atoms found in json file ",natoms

        LV = numpy.zeros( (3,3) )
        
	lv_str = json_data['metadata']["atomic"]["latticevector"].split()
        LV[0][0] = float(lv_str[0])
        LV[0][1] = float(lv_str[1])
        LV[0][2] = float(lv_str[2])
        LV[1][0] = float(lv_str[3])
        LV[1][1] = float(lv_str[4])
        LV[1][2] = float(lv_str[5])
        LV[2][0] = float(lv_str[6])
        LV[2][1] = float(lv_str[7])
        LV[2][2] = float(lv_str[8])
        
        #atoms1 = ParticleContainer()

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
	    vel_str = json_data['metadata']["atomic"]["vel"][atom_i].split()
	    v_i = numpy.array( [float(vel_str[0]),float(vel_str[1]),float(vel_str[2])] )
	    VEL.append( v_i )
	    ATYPE.append( json_data['metadata']["atomic"]["fftype"][atom_i]  )
	    AMASS.append( json_data['metadata']["atomic"]["mass"][atom_i]  )
	    MOLNUMB.append( json_data['metadata']["atomic"]["chain"][atom_i]  )
	    RING_NUMB.append( json_data['metadata']["atomic"]["ring"][atom_i]  )
	    RESID.append( json_data['metadata']["atomic"]["resname"][atom_i]  )
	    RESN.append( json_data['metadata']["atomic"]["residue"][atom_i]  )
	    CHARN.append( json_data['metadata']["atomic"]["chrargegroup"][atom_i]  )
	    
	    success = 1
	    
    return (ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,VEL,ATYPE,AMASS,MOLNUMB,RING_NUMB,RESID,RESN,CHARN,LV,success)

def append_atomic(json_data,ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,VEL,ATYPE,AMASS,MOLNUMB,RING_NUMB,RESID,RESN,CHARN,LV):
    """
    Append atomic data to a  OPV database json file 
    """
    import file_io

    debug = 1 
    
    metadata_found,atomicdata_found = check_atomic(json_data)

    if( debug):
        print  "metadata_found",metadata_found
        print  "atomicdata_found",atomicdata_found
    
    if( not metadata_found ):
	# Add metadata section if not in json_data 
	meta_data = {}
	json_data["metadata"] = meta_data 
	
	
    if( not atomicdata_found ):
	# Add metadata section if not in json_data 
	atomic_data = {}
        json_data['metadata']["atomic"] = atomic_data
	
        atomic_data["natoms"] = len(ELN) 
        lv_string = str( "%f %f %f %f %f %f %f %f %f " % (LV[0,0],LV[0,1],LV[0,2],LV[1,0],LV[1,1],LV[1,2],LV[2,0],LV[2,1],LV[2,2]))
        atomic_data["latticevector"] = lv_string
        atomic_data["element"] = []
        atomic_data["asymb"] = []
        atomic_data["q"] = []
        atomic_data["ctype"] = []    
        atomic_data["unitnumb"] = []
        atomic_data["pos"] = []
        atomic_data["unittype"] = []

        atomic_data["vel"] = []
        atomic_data["fftype"] = []
        atomic_data["mass"] = []
        atomic_data["chain"] = []
        atomic_data["ring"] = []
        atomic_data["resname"] = []
        atomic_data["residue"] = []
        atomic_data["chrargegroup"] = []
	
        for atom_i in range( len(ELN) ):
            atomic_data["element"].append( ELN[atom_i] )
            atomic_data["asymb"].append( ASYMB[atom_i] )
            atomic_data["ctype"].append( CTYPE[atom_i] )
            atomic_data["q"].append(CHARGES[atom_i])
            atomic_data["unitnumb"].append( UNITNUMB[atom_i] )
            pos_str = str( "%f %f %f "% (R[atom_i][0],R[atom_i][1],R[atom_i][2]) )
            atomic_data["pos"].append( pos_str )
            atomic_data["unittype"].append( UNITTYPE[atom_i] )
            
            vel_str = str( "%f %f %f "% (VEL[atom_i][0],VEL[atom_i][1],VEL[atom_i][2]) )
            atomic_data["vel"].append( vel_str )
            atomic_data["fftype"].append( ATYPE[atom_i] )
            atomic_data["mass"].append( AMASS[atom_i] )
            atomic_data["chain"].append( MOLNUMB[atom_i] )
            atomic_data["ring"].append(  RING_NUMB[atom_i] )
            atomic_data["resname"].append( RESID[atom_i] )
            atomic_data["residue"].append( RESN[atom_i] )
            atomic_data["chrargegroup"].append( CHARN[atom_i] )
	 
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
    
            vel_str = str( "%f %f %f "% (VEL[atom_i][0],VEL[atom_i][1],VEL[atom_i][2]) )
            json_data['metadata']["atomic"]["vel"].append( vel_str )
            json_data['metadata']["atomic"]["fftype"].append( ATYPE[atom_i] )
            json_data['metadata']["atomic"]["mass"].append( AMASS[atom_i] )
            json_data['metadata']["atomic"]["chain"].append( MOLNUMB[atom_i] )
            json_data['metadata']["atomic"]["ring"].append(  RING_NUMB[atom_i] )
            json_data['metadata']["atomic"]["resname"].append( RESID[atom_i] )
            json_data['metadata']["atomic"]["residue"].append( RESN[atom_i] )
            json_data['metadata']["atomic"]["chrargegroup"].append( CHARN[atom_i] )
	    
    return json_data


def write_connections(json_data,BONDS):
    """
    Append connections data to a  OPV database json file 
    """
    import file_io

    debug = 1 

    connections_data = {}
    json_data['metadata']["connections"] = connections_data

    connections_data["nbonds"] = len(BONDS) 
    connections_data["bonds"] = []
    

    for bond_i in range( len(BONDS) ):
        b_ij = str( "%d %d "% (BONDS[bond_i][0],BONDS[bond_i][1] ) )
        connections_data["bonds"].append( b_ij )

    return json_data

        
def write_atomic(json_data,ELN,ASYMB,CTYPE,CHARGES,UNITNUMB,UNITTYPE,R,VEL,ATYPE,AMASS,MOLNUMB,RING_NUMB,RESID,RESN,CHARN,LV):
    """
    Append atomic data to a  OPV database json file 
    """
    import file_io

    debug = 1 

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

    atomic_data["vel"] = []
    atomic_data["fftype"] = []
    atomic_data["mass"] = []
    atomic_data["chain"] = []
    atomic_data["ring"] = []
    atomic_data["resname"] = []
    atomic_data["residue"] = []
    atomic_data["chrargegroup"] = []

    for atom_i in range( len(ELN) ):
        atomic_data["element"].append( ELN[atom_i] )
        atomic_data["asymb"].append( ASYMB[atom_i] )
        atomic_data["ctype"].append( CTYPE[atom_i] )
        atomic_data["q"].append(CHARGES[atom_i])
        atomic_data["unitnumb"].append( UNITNUMB[atom_i] )
        pos_str = str( "%f %f %f "% (R[atom_i][0],R[atom_i][1],R[atom_i][2]) )
        atomic_data["pos"].append( pos_str )
        atomic_data["unittype"].append( UNITTYPE[atom_i] )

        vel_str = str( "%f %f %f "% (VEL[atom_i][0],VEL[atom_i][1],VEL[atom_i][2]) )
        atomic_data["vel"].append( vel_str )
        atomic_data["fftype"].append( ATYPE[atom_i] )
        atomic_data["mass"].append( AMASS[atom_i] )
        atomic_data["chain"].append( MOLNUMB[atom_i] )
        atomic_data["ring"].append(  RING_NUMB[atom_i] )
        atomic_data["resname"].append( RESID[atom_i] )
        atomic_data["residue"].append( RESN[atom_i] )
        atomic_data["chrargegroup"].append( CHARN[atom_i] )
	 
    lv_string = str( "%f %f %f %f %f %f %f %f %f " % (LV[0,0],LV[0,1],LV[0,2],LV[1,0],LV[1,1],LV[1,2],LV[2,0],LV[2,1],LV[2,2]))
    atomic_data["latticevector"] = lv_string
	         

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
