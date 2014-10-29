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

# Stream modules 
import lammps , gromacs , xmol

sperator_line = "---------------------------------------------------------------------\n"


def file_exists(filename):
    try:
        with open(filename) as f:
            return True
    except IOError:
        return False


def struc_array_json(struc_array,param_array,json_list):
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
		param_array.append(param_i)

		
    return struc_array,param_array


def struc_array_gromacs(struc_array,gro_list,top_list,paramC):

    """
    Loop over list of gro and top files and place each structure container into a list

    Arguments
      gro_list (list) list of gro files
      top_list (list) list of top files
      
    Returns
      struc_array (list) of structure objects
      
    """
    import gromacs

    verbose = True

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
                param_i = ParameterContainer()
                struc_i = gromacs.read_gro(struc_i,gro_files[g_cnt])
                struc_i,param_i,ljmixrule = gromacs.read_top(struc_i,param_i,top_files[g_cnt])

                if( verbose ):
                    print sperator_line
                    print sperator_line

                    print " Gromacs read in gro file %s with top file %s "%(gro_files[g_cnt],top_files[g_cnt])
                    print  struc_i
                    print param_i

                paramC += param_i

		struc_array.append(struc_i)
		
    return struc_array,paramC

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

    # Place paticle data in sperate data structure
    if( 'oligomers' in json_ref):
        for oligo_list in json_ref['oligomers']:
            col = oligo_list #.split()
            json_file = "%s/%s.json"%(col[0],col[1])
            json_list.append(json_file)
    else:
        print " No oligomers found in  %s "%(ref_file)     


    # Place paticle data in sperate data structure
    if( 'systems' in json_ref):
        for oligo_list in json_ref['systems']:
            col = oligo_list #.split()
            json_file = "%s/%s.json"%(col[0],col[1])
            json_list.append(json_file)
    else:
        print " No systems found in  %s "%(ref_file)     

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
        particle_data["position"].append( ptclObj.position )  # In angstroms 
        particle_data["mass"].append( ptclObj.mass )          # In AMU  
        particle_data["charge"].append( ptclObj.charge )       # In e 
        # Dictionary items 
        particle_data["chain"].append( ptclObj.tagsDict["chain"] )
        particle_data["ring"].append( ptclObj.tagsDict["ring"] )
        particle_data["resname"].append( ptclObj.tagsDict["resname"] )
        particle_data["residue"].append( ptclObj.tagsDict["residue"] )
        particle_data["linkid"].append( ptclObj.tagsDict["linkid"] )        
        particle_data["fftype"].append( ptclObj.tagsDict["fftype"] )        

        particle_data["symbol"].append( ptclObj.tagsDict["symbol"] )        
        particle_data["number"].append( ptclObj.tagsDict["number"] )        
        particle_data["cov_radii"].append( ptclObj.tagsDict["cov_radii"] )        # In angstroms 
        particle_data["vdw_radii"].append( ptclObj.tagsDict["vdw_radii"] )        # In angstroms 

    bonded_data["bonds"] = []
    for b_i,bondObj in  strucC.bondC:
        pt_i = bondObj.pgid1
        pt_j = bondObj.pgid2
        type_i = bondObj.type 
        bonded_data["bonds"].append( [type_i, pt_i,pt_j])

    bonded_data["angles"] = []
    for a_i,angleObj in  strucC.angleC:
        pt_k = angleObj.pgid1
        pt_i = angleObj.pgid2
        pt_j = angleObj.pgid3
        type_i = angleObj.type 
        bonded_data["angles"].append( [type_i,pt_k,pt_i,pt_j])

    bonded_data["dihedrals"] = []
    for d_i,dihObj in  strucC.dihC:
        pt_k = dihObj.pgid1
        pt_i = dihObj.pgid2
        pt_j = dihObj.pgid3
        pt_l = dihObj.pgid4
        type_i = dihObj.type 
        bonded_data["dihedrals"].append( [type_i,pt_k,pt_i,pt_j,pt_l])

    ljtypC_p = paramC.ljtypC
    btypC_p = paramC.btypC
    atypC_p = paramC.atypC
    dtypC_p = paramC.dtypC

    ffparameter_data["ljtypes"] = []
    for lj_p, ljObj_p  in ljtypC_p:
        ptype1 = ljObj_p.ptype1
        mass = ljObj_p.get_mass()
        epsilon = ljObj_p.get_epsilon()
        sigma = ljObj_p.get_sigma()
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
    debug = False 
    ffparameter_data["dihtypes"] = []
    for dtyp_p, dtypObj_p  in dtypC_p:
        dtype = dtypObj_p.type
        ptype1 = dtypObj_p.ptype1
        ptype2 = dtypObj_p.ptype2
        ptype3 = dtypObj_p.ptype3
        ptype4 = dtypObj_p.ptype4

        if(debug):
            print "  putstruc_json dih type ",dtyp_p,dtype
            
        if( dtype == "harmonic" ):
            d = dtypObj_p.d
            mult = dtypObj_p.mult
            theat_s = dtypObj_p.theat_s
            kb = dtypObj_p.kb
            ffparameter_data["dihtypes"].append( [dtype,ptype1,ptype2,ptype3,ptype4,d,mult,theat_s,kb])
        if( dtype == "opls" ):
            k1 = dtypObj_p.k1
            k2 = dtypObj_p.k2
            k3 = dtypObj_p.k3
            k4 = dtypObj_p.k4
            ffparameter_data["dihtypes"].append( [dtype,ptype1,ptype2,ptype3,ptype4,k1,k2,k3,k4])
        if( dtype == "rb" ):
            C0 = dtypObj_p.C0
            C1 = dtypObj_p.C1
            C2 = dtypObj_p.C2
            C3 = dtypObj_p.C3
            C4 = dtypObj_p.C4
            C5 = dtypObj_p.C5
            ffparameter_data["dihtypes"].append( [dtype,ptype1,ptype2,ptype3,ptype4,C0,C1,C2,C3,C4,C5])
        #print dtyp_p ,dtypObj_p.ptype1 ,dtypObj_p.ptype2,dtypObj_p.ptype3,dtypObj_p.ptype4

    if( debug): sys.exit(" dih type write debug ")

    return json_data


def write_json(strucC,paramC,dir_id,output_id ):
    """
    Write structure information into a json file
    """

    json_data = {}
    json_data = putstruc_json(strucC,paramC,json_data)


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
    if( 'structure' in json_data):
        struc_data = json_data["structure"]        
        # Read in lattice vectors
        if( 'latvector' in struc_data):
            strucC.latvec = []
            lv_array = struc_data["latvector"].split()
            strucC.latvec.append(  np.array( [float(lv_array[0]),float(lv_array[1]),float(lv_array[2])] ) )
            strucC.latvec.append(  np.array( [float(lv_array[3]),float(lv_array[4]),float(lv_array[5])] ) )
            strucC.latvec.append(  np.array( [float(lv_array[6]),float(lv_array[7]),float(lv_array[8])] ) )
        else:
            print " No latvector data "
            
        if( 'particle' in struc_data):
            particle_data = struc_data["particle"]

            # Create structure container for particles
            for p_i in range( len( particle_data["number_id"])):


                r_i = particle_data["position"][p_i]
                m_i = float(particle_data["mass"][p_i])
                q_i = float(particle_data["charge"][p_i])

                if( 'symbol' in particle_data):
                    atomic_symb = str( particle_data["symbol"][p_i] )
                else:
                    atomic_symb = str( particle_data["type"][p_i] )
                    
                type_i = str( particle_data["type"][p_i] )
                # Create particle
                pt_i = Particle( r_i,type_i,q_i,m_i )
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
                tagsD={"symbol":atomic_symb,"chain":chain_i,"ring":ring_i,"resname":resname_i,"residue":residue_i,"linkid":linkid_i,"fftype":fftype_i,"symbol":atomic_symb,"number":el.number,"mass":el.mass,"cov_radii":el.cov_radii,"vdw_radii":el.vdw_radii}
                pt_i.setTagsDict(tagsD)
                strucC.ptclC.put(pt_i)

        else:
            print " No particle structure data "        

        if( 'bonded' in struc_data):
            bonded_data =  struc_data["bonded"]

            # Read in bonds
            for b_indx in range( len(bonded_data["bonds"] )):
                type_i = str(bonded_data["bonds"][b_indx][0] )
                a_i = int(bonded_data["bonds"][b_indx][1] )
                a_j = int(bonded_data["bonds"][b_indx][2] )
                r_i = 0.0
                b_i = Bond( a_i, a_j,r_i,type_i )            
                strucC.bondC.put(b_i)

            # Read in angles
            for a_indx in range( len(bonded_data["angles"] )):
                type_i = str(bonded_data["angles"][a_indx][0] )
                a_k = int(bonded_data["angles"][a_indx][1] )
                a_i = int(bonded_data["angles"][a_indx][2] )
                a_j = int(bonded_data["angles"][a_indx][3] )
                theta_i = 0.0
                a_i = Angle( a_k,a_i, a_j,theta_i,type_i )            
                strucC.angleC.put(a_i)

            # Read in dihedrals
            for d_indx in range( len(bonded_data["dihedrals"] )):
                type_i = str(bonded_data["dihedrals"][d_indx][0] )
                a_k = int(bonded_data["dihedrals"][d_indx][1] )
                a_i = int(bonded_data["dihedrals"][d_indx][2] )
                a_j = int(bonded_data["dihedrals"][d_indx][3] )
                a_l = int(bonded_data["dihedrals"][d_indx][4] )
                theta_i = 0.0
                d_i = Dihedral( a_k,a_i, a_j,a_l,theta_i,type_i )  
                strucC.dihC.put(d_i)

        else:
            print " No bonded structure data "        

        if( 'ffparameter' in struc_data):
            ffparameter_data =  struc_data["ffparameter"]

            ljtypC_p = paramC.ljtypC
            btypC_p = paramC.btypC
            atypC_p = paramC.atypC
            dtypC_p = paramC.dtypC

            for lj_list in ffparameter_data["ljtypes"] :
                ljObj_p = ljtype( str(lj_list[0]) )
                ljObj_p.setmass(float(lj_list[1]) )
                ljObj_p.setparam( float(lj_list[2]),float(lj_list[3]) )
                ljtypC_p.put(ljObj_p)

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
                    dtypObj_p.setrb(float(dtyp_list[5]),float(dtyp_list[6]),float(dtyp_list[7]),float(dtyp_list[8]),float(dtyp_list[9]),float(dtyp_list[10]) )


                dtypC_p.put(dtypObj_p)
        
        else:
            print " No ffparameter structure data "        

        
    else:
        print " No structure data "        

        
    return strucC,paramC,json_data
    

"""
No longer used as new functions

gromacs.read_top
gromacs.read_gro

are working 

def get_gromacs(strucC, gro_file,top_file):  # Move out of class
    Read in structure information from gromacs files

    Args:
        gro_file (str) gro file
        top_file (str) top file


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
"""

def create_search(f_id,f_symb,f_chain,f_ring,f_resname,f_residue,f_linkid,f_fftype):
    """
    Create a dictionary to pass to particle search
    """
    
    search_i = {}

    if( len( f_symb ) ):
        search_i["symbol"] = []
        for id_s in f_symb.split():
            search_i["symbol"].append(id_s)
    if( len( f_chain ) ):
        search_i["chain"] = []
        for id_s in f_chain.split():
            search_i["chain"].append(id_s)
    if( len( f_ring ) ):
        search_i["ring"] = []
        for id_s in f_ring.split():
            search_i["f_ring"].append(id_s)
    if( len( f_resname ) ):
        search_i["resname"] = []
        for id_s in f_resname.split():
            search_i["resname"].append(id_s)
    if( len( f_residue ) ):
        search_i["residue"] = []
        for id_s in f_residue.split():
            search_i["residue"].append(id_s)
    if( len( f_linkid  ) ):
        search_i["linkid"] = []
        for id_s in f_linkid.split():
            search_i["linkid"].append(id_s)
    if( len( f_fftype ) ):
        search_i["fftype"] = []
        for id_s in f_fftype.split():
            search_i["fftype"].append(id_s)
            
    return search_i
    
def getstrucC(struc_o,param_o, in_json, in_gro , in_top,in_itp, in_data, in_xmol,xmol_format ):
    """
    Read in Structure data from simulation output files
    """
    verbose = True
    #
    # Read in json file
    #
    if( len(in_json) ):
        if(  verbose ):
            print  "     - Reading in ",in_json
        json_data_i = getsys_json(struc_o,param_o, in_json)
    #
    # Read gro file 
    #
    if( len(in_gro) ):
        if( verbose ): print  "     GROMACS .gro file ",in_gro
        struc_o = gromacs.read_gro(struc_o,in_gro)
    #
    # Read in top file
    #
    if( len(in_top) ):
        if( verbose ): print  "     GROMACS .top file ",in_top
        struc_o,param_o,ljmixrule = gromacs.read_top(struc_o,param_o,in_top)
    #
    # Read in itp file
    #
    if( len(in_itp) ):
        if( verbose ): print  "     - Reading in ",in_itp
        param_o = gromacs.read_itp(param_o, in_itp, ljmixrule)
    # 
    # Read lammps data file 
    #
    if( len(in_data) ):
        if( verbose ): print  "     LAMMPS data file ",in_data            
        struc_o,param_o = lammps.read_lmpdata(struc_o,param_o,in_data)
    # 
    # Read xmol file 
    #
    if( len(in_xmol) ):
        if( verbose ): print  "     xmol data file ",in_xmol
        ptclC_array = xmol.read(in_xmol,xmol_format)
        # Set last set of positions as the particles in the structure 
        struc_o.ptclC = ptclC_array[:-1]

        for pid,pt_i in  struc_o.ptclC:
            print pid,pt_i.type, pt_i.postion
        sys.exit("debug xmol read in 1")
    #
    # HOOMD input file 
    #
    #a = hoomd_xml.hoomd_xml(sys.argv[1])

    return struc_o,param_o

def write_calcinfo(calc_type,loc_dir,job_dir,job_name,status):
    """
    Write info for calculation
    """

    calc_info_list = [calc_type,loc_dir,job_dir,job_name,status]

    return calc_info_list

    
def read_calcinfo(calc_info_list):
    """
    Read info for calculation
    """

    calc_type  = calc_info_list[0].strip()
    loc_dir   = calc_info_list[1].strip()
    job_dir  = calc_info_list[2].strip()
    job_name  = calc_info_list[3].strip()
    status    = calc_info_list[4].strip()

    return calc_type,loc_dir,job_dir,job_name,status
    
