#! /usr/bin/env python
"""

Radial distribution  code

 length - Angstroms
 mass   - AMU
 volume - Angstroms^3
"""

# Dr. Travis Kemper
# NREL
# Initial Date 6/30/2014
# travis.kemper@nrel.gov


# Scott's new classes 
from particles import Particle
from particles import ParticleContainer
from bonds     import BondContainer
from structureContainer import StructureContainer
import mpiNREL
import json
import numpy as np
import datetime


def get_options():
    """
    Set options
    """
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-p","--ptime", dest="ptime", default=False,action="store_true", help="Print performance information  ")
    parser.add_option("-o","--output_id", dest="output_id", default="rdf",type="string",help=" prefix for output files  ")
    
    # Files with particle information 
    parser.add_option("-j","--in_json", dest="in_json", type="string", default="", help="Input json file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    
    # Trajectory file
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")    
    
    # Bins
    parser.add_option("--r_cut", dest="r_cut", type=float, default=10.0, help=" Cut off radius in angstroms ")
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.10, help=" Bin size in angstroms")
    

    # Searchable properties
    # I
    parser.add_option("--id_i", dest="id_i", type="string", default="", help=" select atoms of group i by id number  ")    
    parser.add_option("--symb_i", dest="symb_i", type="string", default="", help=" select atoms of group i by (atomic) symbol   ")    
    parser.add_option("--chains_i", dest="chains_i", type="string", default="", help="select atoms of group i by chain number  ")    
    parser.add_option("--ring_i", dest="ring_i", type="string", default="", help="select atoms of group i by particles in a ring   ")    
    parser.add_option("--resname_i", dest="resname_i", type="string", default="", help="select atoms of group i by residue name  ")    
    parser.add_option("--residue_i", dest="residue_i", type="string", default="", help="select atoms of group i by resudue number  ")    
    parser.add_option("--linkid_i", dest="linkid_i", type="string", default="", help="select atoms of group i by  link id ")    
    parser.add_option("--fftype_i", dest="fftype_i", type="string", default="", help="select atoms of group i by force field type  ")
    # J
    parser.add_option("--id_j", dest="id_j", type="string", default="", help=" select atoms of group i by id number  ")    
    parser.add_option("--symb_j", dest="symb_j", type="string", default="", help=" select atoms of group i by (atomic) symbol   ")    
    parser.add_option("--chains_j", dest="chains_j", type="string", default="", help="select atoms of group i by chain number  ")    
    parser.add_option("--ring_j", dest="ring_j", type="string", default="", help="select atoms of group i by particles in a ring   ")    
    parser.add_option("--resname_j", dest="resname_j", type="string", default="", help="select atoms of group i by residue name  ")    
    parser.add_option("--residue_j", dest="residue_j", type="string", default="", help="select atoms of group i by resudue number  ")    
    parser.add_option("--linkid_j", dest="linkid_j", type="string", default="", help="select atoms of group i by  link id ")    
    parser.add_option("--fftype_j", dest="fftype_j", type="string", default="", help="select atoms of group i by force field type  ")    
    
    # Frames

    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=0, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    
        
    (options, args) = parser.parse_args()
        
    return options, args
   

def getsys_json(json_file):
    """
    Return a new Structure object with partcleID's in input list
        
    Args:
        ptclIDList (list) global particles ID's for which to return structure

    Return:
        New Structure() object. IDs in new object are unique
    """
    import json 

    f = open(json_file, 'r')
    json_data = json.load(f)
    f.close()

    # Place paticle data in sperate data structure 
    particle_data = json_data["structure"]["particle"]

    # Create structure container for particles 
    struc_i = ParticleContainer()

    for p_i in range( len( particle_data["number_id"])):
        r_i = particle_data["position"][p_i]
        atomic_symb = str( particle_data["type"][p_i] )
        m_i = float(particle_data["mass"][p_i])
        q_i = float(particle_data["charge"][p_i])
        # Create particle
        pt_i = Particle( r_i,atomic_symb,q_i,m_i )
        # Find needed tags
        chain_i = particle_data["chain"][p_i]
        ring_i = particle_data["ring"][p_i]
        resname_i = particle_data["resname"][p_i]
        residue_i = particle_data["residue"][p_i]
        linkid_i = particle_data["linkid"][p_i]
        fftype_i = particle_data["fftype"][p_i]
        # _i = particle_data[""][p_i]
        # Add particle to structure 
        tagsD = {"chain":chain_i,"ring":ring_i,"resname":resname_i,"residue":residue_i,"linkid":linkid_i,"fftype":fftype_i}
        pt_i.setTagsDict(tagsD)
        struc_i.put(pt_i)

    twobody_data =  json_data["structure"]["twobody"]
    bond_i = BondContainer()
    for b_indx in rang( len(twobody_data["bonds"] )):
        a_i = twobody_data["bonds"][0]
        a_j = twobody_data["bonds"][0]
        self.bondC.put(b_i)

    # Put paticles and bonds into system 
    system_i = StructureContainer(ptclC=struc_i,bondC=bond_i)
    

    return system_i

def  create_search(f_id,f_symb,f_chain,f_ring,f_resname,f_residue,f_linkid,f_fftype):
    """
    Create a dictionary to pass to particle search
    """
    
    search_i = {}

    if( len( f_symb ) ):
        search_i["type"] = []
        for id_s in f_symb.split():
            search_i["type"].append(id_s)
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
            search_i["f_resname"].append(id_s)
    if( len( f_residue ) ):
        search_i["residue"] = []
        for id_s in f_residue.split():
            search_i["f_residue"].append(id_s)
    if( len( f_linkid  ) ):
        search_i["linkid"] = []
        for id_s in f_linkid.split():
            search_i["f_linkid"].append(id_s)
    if( len( f_fftype ) ):
        search_i["fftype"] = []
        for id_s in f_fftype.split():
            search_i["f_fftype"].append(id_s)
            
    return search_i
    
def main():
    """
    Calculate radial distribution (RDF) based on specified atom types

    Input:
        - json file ".json" containing all the information for structure
        - frames to average RDF over
            - lammps xyz file ".xyz"
        - calculation specifications 
            - atomic groups to calculate RDF between
            - which frames to use 
    Output:
        - data file ".dat" containing the rdf
        - log file ".log" containing calculation information
            the amount of information in the log file can be increased using the verbose option ( -v ) 
        
    """

    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    options, args = get_options()

    # Read in system data from json file
    if( len(options.in_json)):
        system_i = getsys_json(options.in_json)

    # Read in system data from json file
    # 
    # if( len(in_data)):
    #    system_i = getsys_lmpdata(options.in_data)
        
    # Get paticle and bond structures
    struc_o = system_i.ptclC
    bond_i  = system_i.bondC
    
    p.barrier()
    
    #
    # Open output files 
    #
    if( rank == 0 ):

        log_file = options.output_id + ".log"
        log_out = open(log_file,"w") 

        dat_file = options.output_id + ".dat"
        dat_out = open(dat_file,"w") 
        dat_out.write("#   Input ")

    p.barrier()

     
    # Print system properties
    if( rank == 0 ):
        sys_prop = system_i.printprop()
        print sys_prop
        log_out.write(str(sys_prop))
        
    # Create sub lists for RDF groups i and j
    search_i = create_search(options.id_i,options.symb_i,options.chains_i,options.ring_i,options.resname_i,options.residue_i,options.linkid_i,options.fftype_i)    
    list_i = struc_o.getParticlesWithTags(search_i)

    search_j = create_search(options.id_j,options.symb_j,options.chains_j,options.ring_j,options.resname_j,options.residue_j,options.linkid_j,options.fftype_j)    
    list_j = struc_o.getParticlesWithTags(search_j)
    
    if( rank == 0 ):
        if( options.verbose ):
            log_line = "  Particles of group i: %d "%(len(list_i))
            log_out.write(log_line)
            print log_line
            log_line = "---------------------------------------------------------------------"
            log_out.write(log_line)
            print log_line
            log_line = " id# ; type ; ff type ; ring # ; residue name ; residue # ; link id  "
            log_out.write(log_line)
            print log_line
            log_line = "---------------------------------------------------------------------"
            log_out.write(log_line)
            print log_line
            for p_i, ptcl in struc_o(list_i):
                log_line = "   %d %s %s %s %s %s %s "%(p_i,ptcl.type,ptcl.tagsDict["fftype"],ptcl.tagsDict["ring"],ptcl.tagsDict["resname"],ptcl.tagsDict["residue"],ptcl.tagsDict["linkid"])
                log_out.write(log_line)
                print log_line
            log_line = "  Particles of group i: %d "%(len(list_j))
            log_out.write(log_line)
            print log_line
            log_line = "---------------------------------------------------------------------"
            log_out.write(log_line)
            print log_line
            log_line = " id# ; type ; ff type ; ring # ; residue name ; residue # ; link id  "
            log_out.write(log_line)
            print log_line
            log_line = "---------------------------------------------------------------------"
            log_out.write(log_line)
            print log_line
            for p_j, ptcl in struc_o(list_j):
                log_line = "   %d %s %s %s %s %s %s "%(p_j,ptcl.type,ptcl.tagsDict["fftype"],ptcl.tagsDict["ring"],ptcl.tagsDict["resname"],ptcl.tagsDict["residue"],ptcl.tagsDict["linkid"])
                log_out.write(log_line)
                print log_line
            
    # Calculate rdf relate values
    n_bins = int(options.r_cut/options.bin_size)
    sq_r_cut = options.r_cut**2
    rdf_cnt_ij = np.zeros(n_bins+1)    

    # Record initial time

    if( rank == 0  ): 
	t_i = datetime.datetime.now()

    if( len(options.in_lammpsxyz) ):
        
        # Initialize line count
        NP = len( struc_i )
        line_cnt = 0
        frame_cnt = 0
        pt_cnt = 0 
        # Read in file line by line 
        with open(options.in_lammpsxyz) as f:
            for line in f:
                line_cnt += 1
                col = line.split()
                
            if( line_cnt == 1  ):
                # Read first line to record number of particles NP
                NP = int(col[0])
                frame_cnt += 1
                pt_cnt = 0 
            elif( line_cnt > 2 and len(col) >= 4 ):
                # Read lines and update particle postions 
                pt_cnt += 1
                pt_i = struc_i[pt_cnt]
                pt_i.postion =  [ col[1],col[2],col[3] ]
                atomic_symb = col[0]
                #
                #r_i = [ col[1],col[2],col[3] ]
                #part_i = Particle( r_i,atomic_symb,q_i,m_i )
                #struc_i.put(part_i)

            if( line_cnt > 1  and line_cnt > NP + 2):
                pt_cnt = 0 
                if( options.frame_o >= frame_cnt and frame_cnt <= options.frame_f ):
                    if( mod(frame_cnt/options.frame_step) ):
                        rdf_cnt_ij = struc_i.calc_rdf(rdf_cnt_ij)

                frame_cnt += 1
                
            # Calculate rdf for last frame 
            if( options.frame_o >= frame_cnt and frame_cnt <= options.frame_f ):
                if( mod(frame_cnt/options.frame_step) ):
                    rdf_cnt_ij = struc_i.calc_rdf(rdf_cnt_ij)
            

    dat_out.close()
    log_out.close()

if __name__=="__main__":
    main()
   
