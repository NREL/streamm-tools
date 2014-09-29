#! /usr/bin/env python
"""
File manipulations 
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/19/2014
# travis.kemper@nrel.gov

import json
import numpy as np

from structureContainer import StructureContainer
from periodictable import periodictable
from particles     import Particle, ParticleContainer
from bonds         import Bond,     BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer

from parameters    import ParameterContainer
from parameters    import ljtype,   LJtypesContainer
from parameters    import bondtype, BondtypesContainer
from parameters    import angletype,AngletypesContainer
from parameters    import dihtype,  DihtypesContainer



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
                param_i = ParameterContainer()

                struc_i,param_i,json_data_i = getsys_json(struc_i,param_i,json_file)
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
                struc_i, json_data_i = struc_i.get_gromacs(struc_i, gro_files[g_cnt],top_files[g_cnt])
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

    ljtypC_p = paramC.ljtypC
    btypC_p = paramC.btypC
    atypC_p = paramC.atypC
    dtypC_p = paramC.dtypC

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
        ptype4 = dtypObj_p.ptype4
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


def write_json(strucC,dir_id,output_id ):
    """
    Write structure information into a json file
    """

    json_data = {}
    json_data = self.putstruc_json(json_data)


    json_file = dir_id+"/"+output_id + ".json"
    f = open(json_file, 'w')
    json.dump(json_data,f, indent=2)
    f.close()

def getsys_json(strucC,paramC, json_file):
    """
    Read in structure information from json file 

    Args:
        strucC   (StructureContainer) 
        paramC   (ParameterContainer) 
        json_file (json) file with structure information

    Return
        json_data (json structure data )

    """

    # Load periodic table 
    pt = periodictable()

    f = open(json_file, 'r')
    json_data = json.load(f)
    f.close()

    # Place paticle data in sperate data structure 
    struc_data = json_data["structure"]
    particle_data = json_data["structure"]["particle"]
    bonded_data =  json_data["structure"]["bonded"]
    ffparameter_data = json_data["structure"]["ffparameter"]

    # Create structure container for particles


    for p_i in range( len( particle_data["number_id"])):
        r_i = particle_data["position"][p_i]
        atomic_symb = str( particle_data["symbol"][p_i] )
        m_i = float(particle_data["mass"][p_i])
        q_i = float(particle_data["charge"][p_i])
        # Create particle
        pt_i = Particle( r_i,atomic_symb,q_i,m_i )
        # Find needed tags
        chain_i = int( particle_data["chain"][p_i] )
        ring_i = particle_data["ring"][p_i]
        resname_i = particle_data["resname"][p_i]
        residue_i = particle_data["residue"][p_i]
        linkid_i = particle_data["linkid"][p_i]
        fftype_i = particle_data["fftype"][p_i]
        # _i = particle_data[""][p_i]
        # Add particle to structure
        # Get element information
        el = pt.getelementWithSymbol(atomic_symb)
        tagsD={"chain":chain_i,"ring":ring_i,"resname":resname_i,"residue":residue_i,"linkid":linkid_i,"fftype":fftype_i,"symbol":atomic_symb,"number":el.number,"mass":el.mass,"cov_radii":el.cov_radii,"vdw_radii":el.vdw_radii}
        pt_i.setTagsDict(tagsD)
        strucC.ptclC.put(pt_i)

    # Read in bonds
    bondC_i = BondContainer()
    for b_indx in range( len(bonded_data["bonds"] )):
        a_i = int(bonded_data["bonds"][b_indx][0] )
        a_j = int(bonded_data["bonds"][b_indx][1] )
        b_i = Bond( a_i, a_j )            
        strucC.bondC.put(b_i)

    # Read in angles
    angleC_i = AngleContainer()
    for a_indx in range( len(bonded_data["angles"] )):
        a_k = int(bonded_data["angles"][a_indx][0] )
        a_i = int(bonded_data["angles"][a_indx][1] )
        a_j = int(bonded_data["angles"][a_indx][2] )
        a_i = Angle( a_k,a_i, a_j )            
        strucC.angleC.put(a_i)

    # Read in dihedrals
    dihC_i = DihedralContainer()
    for d_indx in range( len(bonded_data["dihedrals"] )):
        a_k = int(bonded_data["dihedrals"][b_indx][0] )
        a_i = int(bonded_data["dihedrals"][b_indx][1] )
        a_j = int(bonded_data["dihedrals"][b_indx][2] )
        a_l = int(bonded_data["dihedrals"][b_indx][3] )
        d_i = Dihedral( a_k,a_i, a_j,a_l )            
        strucC.dihC.put(d_i)

    # Read in lattice vectors
    strucC.latvec = []
    lv_array = struc_data["latvector"].split()
    strucC.latvec.append(  np.array( [float(lv_array[0]),float(lv_array[1]),float(lv_array[2])] ) )
    strucC.latvec.append(  np.array( [float(lv_array[3]),float(lv_array[4]),float(lv_array[5])] ) )
    strucC.latvec.append(  np.array( [float(lv_array[6]),float(lv_array[7]),float(lv_array[8])] ) )



    ljtypC_p = paramC.ljtypC
    btypC_p = paramC.btypC
    atypC_p = paramC.atypC
    dtypC_p = paramC.dtypC

    #print 'ffparameter_data["ljtypes"] ',ffparameter_data["ljtypes"]
    
    for lj_list in ffparameter_data["ljtypes"] :
        print lj_list
        print lj_list[0]
        
        ljObj_p = ljtype( str(lj_list[0]) )
        ljObj_p.setmass(float(lj_list[1]) )
        ljObj_p.setparam( float(lj_list[2]),float(lj_list[3]) )
        ljtypC_p.put(ljObj_p)
        #print ljObj_p

        #    sys.exit("ljtypes debug ")
        
    for btyp_list  in ffparameter_data["bondtypes"]:
        btype = str(btyp_list[0])
        btypObj_p = bondtype( str(btyp_list[1]) , str(btyp_list[2]), btype )
        if( btype == "harmonic" ):
            btypObj_p.setharmonic(float(btyp_list[3]),float(btyp_list[4]) )
        btypC_p.put(btypObj_p)
        
    for atyp_list  in ffparameter_data["angletypes"]:

        
        atype = str(atyp_list[0])
        atypObj_p = angletype( str(atyp_list[1]) , str(atyp_list[2]) , str(atyp_list[3]) , atype )
        if( atype == "harmonic" ):
            atypObj_p.setharmonic(float(atyp_list[4]),float(atyp_list[5]) )
        atypC_p.put(atypObj_p)
        
    for dtyp_list in ffparameter_data["dihtypes"]:
        dtype = str(dtyp_list[0])
        dtypObj_p = dihtype( str(dtyp_list[1]) , str(dtyp_list[2]) , str(dtyp_list[3]) , str(dtyp_list[4]) , dtype )
        if( dtype== "harmonic" ):
            dtypObj_p.setharmonic(float(dtyp_list[5]),float(dtyp_list[6]),float(dtyp_list[8]),float(dtyp_list[7]) )
        if( dtype== "opls" ):
            dtypObj_p.setopls(float(dtyp_list[5]),float(dtyp_list[6]),float(dtyp_list[7]),float(dtyp_list[8]) )
        if( dtype== "rb" ):
            dtypObj_p.setharmonic(float(dtyp_list[5]),float(dtyp_list[6]),float(dtyp_list[7]),float(dtyp_list[8]),float(dtyp_list[9]),float(dtyp_list[10]) )

            
        dtypC_p.put(dtypObj_p)
        

    return strucC,paramC,json_data
    


def get_gromacs(strucC, gro_file,top_file):  # Move out of class
    """
    Read in structure information from gromacs files

    Args:
        gro_file (str) gro file
        top_file (str) top file


    """
    sys.exit("get_gromacs in file io needs to be updated ")
    
    import gromacs , elements

    GTYPE,R,VEL,LV = gromacs.read_gro(gro_file)
    ATYPE,RESN,RESID,GTYPE,CHARN,CHARGES,AMASS,BONDS,ANGLES,DIH,MOLNUMB,MOLPNT,MOLLIST = gromacs.read_top(top_file)
    ASYMB,ELN  = elements.mass_asymb(AMASS)


    for p_i in range( len( ASYMB)):
        r_i =  [ float(R[p_i][0]), float(R[p_i][1]), float(R[p_i][2]) ] 
        atomic_symb = str( ASYMB[p_i] )
        m_i = float(AMASS[p_i])
        q_i = float(CHARGES[p_i])
        # Create particle
        pt_i = Particle( r_i,atomic_symb,q_i,m_i )
        # Find needed tags
        chain_i = int( MOLNUMB[p_i] )
        ring_i = 0 #particle_data["ring"][p_i]
        resname_i = RESID[p_i]
        residue_i = RESN[p_i]
        linkid_i = "UNKNOWN" #particle_data["linkid"][p_i]
        fftype_i =ATYPE[p_i]
        gtype_i = GTYPE[p_i]
        # _i = particle_data[""][p_i]
        # Add particle to structure 
        tagsD = {"chain":chain_i,"ring":ring_i,"resname":resname_i,"residue":residue_i,"linkid":linkid_i,"fftype":fftype_i,"gtype":gtype_i}
        pt_i.setTagsDict(tagsD)
        strucC.ptclC.put(pt_i)

    # Read in bonds
    for b_indx in range( len(BONDS )):
        a_i = int(BONDS[b_indx][0] )
        a_j = int(BONDS[b_indx][1] )
        b_i = Bond( a_i + 1 , a_j + 1  )            
        strucC.bondC.put(b_i)

    # Read in lattice vectors
    strucC.latvec = []
    strucC.latvec.append(  np.array( [float(LV[0][0]),float(LV[0][1]),float(LV[0][2])] ) )
    strucC.latvec.append(  np.array( [float(LV[1][0]),float(LV[1][1]),float(LV[1][2])] ) )
    strucC.latvec.append(  np.array( [float(LV[2][0]),float(LV[2][1]),float(LV[2][2])] ) )

    return strucC
