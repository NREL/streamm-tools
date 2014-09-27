#! /usr/bin/env python
"""
File manipulations 
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/19/2014
# travis.kemper@nrel.gov

import json
from structureContainer import StructureContainer


def file_exists(filename):
    try:
        with open(filename) as f:
            return True
    except IOError:
        return False


def struc_array_json(struc_array,json_list):
    """
    Loop over list of json files and place each structure container into a list

    Arguments
      json_list (list) list of json files
      
    Returns
      struc_array (list) of structure objects
      
    """
    
    # Read in structure containers from json files 
    if( len(json_list)):
        # Loop over json files 
        json_files = json_list.split(',')
        if( len(json_files) > 0 ):
            # Read index files from args
            json_cnt = 0 
            for json_file in json_files:
		
		print " Reading in ",json_file
                json_cnt += 1 
                struc_i = StructureContainer()
                json_data_i = struc_i.getsys_json(json_file)
		struc_array.append(struc_i)

		
    return struc_array


def struc_array_gromacs(struc_array,gro_list,top_list):

    """
    Loop over list of gro and top files and place each structure container into a list

    Arguments
      gro_list (list) list of gro files
      top_list (list) list of top files
      
    Returns
      struc_array (list) of structure objects
      
    """

    # Read in structure containers from json files 
    if( len(gro_list) > 0  and len(top_list) > 0 ):
        # Loop over json files 
        gro_files = gro_list.split(',')
        top_files = top_list.split(',')
        # Read index files from args
        if( len(gro_files)  == len(top_files) ):
            
            for g_cnt in range( len(gro_files)):
		
                print " Reading in gro file ",gro_files[g_cnt]," with top file ",top_files[g_cnt]
                struc_i = StructureContainer()
                json_data_i = struc_i.get_gromacs(gro_files[g_cnt],top_files[g_cnt])
		struc_array.append(struc_i)
		
    return struc_array


def read_ref(ref_file):
    """
    Read in reference json file and place all the json files from the oligmers tag in a list

    Arguments
      ref_file (json) reference json file with a oligmers tag

    Returns
       json_list (list) 
    
    """
    
    f = open(ref_file, 'r')
    json_ref = json.load(f)
    f.close()

    json_list = []

    for oligo_list in json_ref['oligomers']:
        col = oligo_list #.split()
        json_file = "%s/%s.json"%(col[0],col[1])
        json_list.append(json_file)

    return json_list
        
def parse_jsons(jsons_input):
    """
    Parse a list of json files
    """

    json_list = []
    
    for json_file in jsons_input.split():
        json_list.append(json_file)

    return json_list


def make_json_list(verbose,ref,json):
    """
    Use reference file or list of json files to create a list of json files

    """

    # Get json files from reference file
    if( len(ref) > 0 ):
        if( verbose ):
            log_line = " Getting json list from reference file  %s "%(ref)
            print log_line
        json_list = read_ref(ref)
    
    # Get json files from json option 
    if( len(json) > 0 ):
        if( verbose ):
            log_line = " Getting list of json files from json option  %s "%(json)
            print log_line
        json_list = parse_jsons( json )

    return json_list
    


    def putstruc_json(strucC, paramC, json_data ):  # Move out of class
        """
        Write a structure into json file
                
        Args:
            json_data (json) json data structure 
            json_file (srt) name of json file 

        """        
        # Initialize json data 
        #   
        struc_data = {}        # Data for entire structure  
        particle_data = {}     # Data for particles and positions 
        bonded_data = {}   # Data for connections between particles (two body interactions) 
        ffparameter_data = {}   # Data for connections between particles (two body interactions) 
        
        json_data["structure"] = struc_data
        struc_data["particle"] = particle_data
        struc_data["bonded"] = bonded_data
        struc_data["ffparameter"] = ffparameter_data
    
	# Structure data
        lv_string = str( "%f %f %f %f %f %f %f %f %f " % ( strucC.latvec[0][0], strucC.latvec[0][1], strucC.latvec[0][2], strucC.latvec[1][0], strucC.latvec[1][1], strucC.latvec[1][2], strucC.latvec[2][0], strucC.latvec[2][1], strucC.latvec[2][2]))

	struc_data["latvector"] = lv_string

	# Particle data 
        particle_data["number_id"] = []
        particle_data["type"] = []
        particle_data["position"] = []
        particle_data["mass"] = []
        particle_data["charge"] = []
        particle_data["chain"] = []
        particle_data["ring"] = []
        particle_data["resname"] = []
        particle_data["residue"] = []
        particle_data["linkid"] = []
        particle_data["fftype"] = []
        
        particle_data["symbol"] = []
        particle_data["number"] = []
        particle_data["cov_radii"] = []
        particle_data["vdw_radii"] = []

	# Loop over particles and ad them to json data 
        for  pid, ptclObj in strucC.ptclC:
            particle_data["number_id"].append(pid )
            particle_data["type"].append( ptclObj.type )
            particle_data["position"].append( ptclObj.position )
            particle_data["mass"].append( ptclObj.mass )
            particle_data["charge"].append( ptclObj.charge )
            # Dictionary items 
            particle_data["chain"].append( ptclObj.tagsDict["chain"] )
            particle_data["ring"].append( ptclObj.tagsDict["ring"] )
            particle_data["resname"].append( ptclObj.tagsDict["resname"] )
            particle_data["residue"].append( ptclObj.tagsDict["residue"] )
            particle_data["linkid"].append( ptclObj.tagsDict["linkid"] )        
            particle_data["fftype"].append( ptclObj.tagsDict["fftype"] )        

            particle_data["symbol"].append( ptclObj.tagsDict["symbol"] )        
            particle_data["number"].append( ptclObj.tagsDict["number"] )        
            particle_data["cov_radii"].append( ptclObj.tagsDict["cov_radii"] )        
            particle_data["vdw_radii"].append( ptclObj.tagsDict["vdw_radii"] )        

        bonded_data["bonds"] = []
        for b_i,bondObj in  strucC.bondC:
            pt_i = bondObj.pgid1
            pt_j = bondObj.pgid2
            bonded_data["bonds"].append( [pt_i,pt_j])
            
        bonded_data["angles"] = []
        for a_i,angleObj in  strucC.angleC:
            pt_k = angleObj.pgid1
            pt_i = angleObj.pgid2
            pt_j = angleObj.pgid3
            bonded_data["angles"].append( [pt_k,pt_i,pt_j])
              
        bonded_data["dihedrals"] = []
        for d_i,dihObj in  strucC.dihC:
            pt_k = dihObj.pgid1
            pt_i = dihObj.pgid2
            pt_j = dihObj.pgid3
            pt_l = dihObj.pgid4
            bonded_data["dihedrals"].append( [pt_k,pt_i,pt_j,pt_l])

        ffparameter_data["ljtypes"] = []
        for lj_p, ljObj_p  in ljtypC_p:
            ptype1 = ljObj_p.ptype1
            mass = ljObj_p.mass
            epsilon = ljObj_p.epsilon
            sigma = ljObj_p.epsilon
            ffparameter_data["ljtypes"].append( [ptype1,mass,epsilon,sigma])
        ffparameter_data["bondtypes"] = []
        for btyp_p, btypObj_p  in btypC_p:
            btype = btypObj_p.type
            ptype1 = btypObj_p.ptype1
            ptype2 = btypObj_p.ptype2
            if( btype == "harmonic" ):
                r0 = btypObj_p.r0
                kb = btypObj_p.kb
                ffparameter_data["bondtypes"].append( [btype,ptype1,ptype2,r0,kb])
                
            #print btyp_p ,btypObj_p.ptype1 ,btypObj_p.ptype2
        ffparameter_data["angletypes"] = []
        for atyp_p, atypObj_p  in atypC_p:
            atype = atypObj_p.type
            ptype1 = atypObj_p.ptype1
            ptype2 = atypObj_p.ptype2
            ptype3 = atypObj_p.ptype3
            if( atype == "harmonic" ):
                theta0 = atypObj_p.theta0
                kb = atypObj_p.kb
                ffparameter_data["angletypes"].append( [atype,ptype1,ptype2,ptype3,theta0,kb])
                
                #print atyp_p ,atypObj_p.ptype1 ,atypObj_p.ptype2,atypObj_p.ptype3
        ffparameter_data["dihtypes"] = []
        for dtyp_p, dtypObj_p  in dtypC_p:
            dtype = dtypObj_p.type
            ptype1 = dtypObj_p.ptype1
            ptype2 = dtypObj_p.ptype2
            ptype3 = dtypObj_p.ptype3
            if( dtype == "harmonic" ):
                d = dtypObj_p.d
                mult = dtypObj_p.mult
                theat_s = dtypObj_p.theat_s
                kb = dtypObj_p.kb
                ffparameter_data["dihtypes"].append( [atype,ptype1,ptype2,ptype3,ptype4,d,mult,theat_s,kb])
            if( dtype == "opls" ):
                k1 = dtypObj_p.k1
                k2 = dtypObj_p.k2
                k3 = dtypObj_p.k3
                k4 = dtypObj_p.k4
                ffparameter_data["dihtypes"].append( [atype,ptype1,ptype2,ptype3,ptype4,k1,k2,k3,k4])
            if( dtype == "rb" ):
                C0 = dtypObj_p.C0
                C1 = dtypObj_p.C1
                C2 = dtypObj_p.C2
                C3 = dtypObj_p.C3
                C4 = dtypObj_p.C4
                C5 = dtypObj_p.C5
                ffparameter_data["dihtypes"].append( [atype,ptype1,ptype2,ptype3,ptype4,C0,C1,C2,C3,C4,C5])
            #print dtyp_p ,dtypObj_p.ptype1 ,dtypObj_p.ptype2,dtypObj_p.ptype3,dtypObj_p.ptype4
            
        return json_data
