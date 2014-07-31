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
    

