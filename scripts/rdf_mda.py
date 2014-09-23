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

from itertools import izip

#MDanalysis
from MDAnalysis import *
from MDAnalysis.core.distances import * ##distance_array
#import MDAnalysis.core.units            # for bulk water density


# Scott's new classes 
from particles import Particle
from particles import ParticleContainer
from bonds     import BondContainer
from structureContainer import StructureContainer
import mpiNREL
import pbcs
import json, math , sys 
import numpy as np
import datetime

import gromacs , lammps


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
    #
    # Input files
    #
    parser.add_option("-j","--in_json", dest="in_json", type="string", default="", help="Input json file, which read in first then over writen by subsequent input files")
    parser.add_option("--in_itp", dest="in_itp",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")

    # Trajectory file
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")    
    parser.add_option("--in_dcd", dest="in_dcd", type="string", default="", help="Input dcd file with atoms listed as atom type numbers")    
    parser.add_option("--in_xtc", dest="in_xtc", type="string", default="", help="Input xtc file with atoms listed as atom type numbers")  

    parser.add_option("--mol_inter",dest="mol_inter", default=False,action="store_true", help="Use only inter molecular rdf's ")
    parser.add_option("--mol_intra",dest="mol_intra", default=False,action="store_true", help="Use only intra molecular rdf's")
    
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
    parser.add_option("--id_j", dest="id_j", type="string", default="", help=" select atoms of group j by id number  ")    
    parser.add_option("--symb_j", dest="symb_j", type="string", default="", help=" select atoms of group j by (atomic) symbol   ")    
    parser.add_option("--chains_j", dest="chains_j", type="string", default="", help="select atoms of group j by chain number  ")    
    parser.add_option("--ring_j", dest="ring_j", type="string", default="", help="select atoms of group j by particles in a ring   ")    
    parser.add_option("--resname_j", dest="resname_j", type="string", default="", help="select atoms of group j by residue name  ")    
    parser.add_option("--residue_j", dest="residue_j", type="string", default="", help="select atoms of group j by resudue number  ")    
    parser.add_option("--linkid_j", dest="linkid_j", type="string", default="", help="select atoms of group j by  link id ")    
    parser.add_option("--fftype_j", dest="fftype_j", type="string", default="", help="select atoms of group j by force field type  ")    
    
    # Frames
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=False,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")

        
    (options, args) = parser.parse_args()
        
    return options, args
   

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
    prtCl_time = False
    debug = False 
    #
    # Formated ouput varables
    #
    sperator_line = "\n---------------------------------------------------------------------"

    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()
 
    options, args = get_options()


    if( rank == 0 ):
        # record initial time 
        t_i = datetime.datetime.now()
        # Open log files 
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w") 

        dat_file = options.output_id + ".dat"
        dat_out = open(dat_file,"w") 
        dat_out.write("#   Input ")

    p.barrier()
    
    #
    #  Initialize blank system 
    # 
    struc_o = StructureContainer()

    if( rank == 0 ):

        log_line = sperator_line
        log_line += "\n  Reading in input files "
        log_line += sperator_line

        log_out.write(log_line)
        print log_line
        
    #
    # Read in json file
    #
    if( len(options.in_json) ):
        if( rank == 0 and options.verbose ):
            print  "     - Reading in ",options.in_json
        json_data_i = struc_o.getsys_json(options.in_json)
    #
    # Read gro file 
    #
    if( len(options.in_gro) ):
        if( options.verbose ): print  "     GROMACS .gro file ",options.in_gro
        struc_o = gromacs.read_gro(struc_o,options.in_gro)
    #
    # Read in top file
    #
    if( len(options.in_top) ):
        if( options.verbose ): print  "     GROMACS .top file ",options.in_top
        struc_o,ljmixrule = gromacs.read_top(struc_o,options.in_top)
    # 
    # Read lammps data file 
    #
    if( len(options.in_data) ):
        if( options.verbose ): print  "     LAMMPS data file ",options.in_data            
        struc_o = lammps.read_lmpdata(struc_o,options.in_data)
        
    p.barrier()
    
    # Get paticle and bond structures
    ptclC_o = struc_o.ptclC
    bondC_o  = struc_o.bondC
        
    # Check options
    if( options.mol_inter and options.mol_intra and rank == 0  ):
	print " Options --mol_inter and --mol_intra are mutually exclusive "
	sys.exit("Error is specified options ")
     
    # Print system properties
    if( rank == 0 ):
        sys_prop = struc_o.printprop()
        print sys_prop
        log_out.write(str(sys_prop))
        
    # Create sub lists for RDF groups i and j
    search_i = create_search(options.id_i,options.symb_i,options.chains_i,options.ring_i,options.resname_i,options.residue_i,options.linkid_i,options.fftype_i)
    if( rank == 0 ):
        if( options.verbose ): print " Searching group i ",search_i
    list_i = ptclC_o.getParticlesWithTags(search_i)
    sum_i = len(list_i)
    
    search_j = create_search(options.id_j,options.symb_j,options.chains_j,options.ring_j,options.resname_j,options.residue_j,options.linkid_j,options.fftype_j)
    if( rank == 0 ):
        if( options.verbose ): print " Searching group j ",search_j
    list_j = ptclC_o.getParticlesWithTags(search_j)
    sum_j = len(list_j)
    
    if( rank == 0 ):
        if( options.verbose ): print " search finished "
    
    if( rank == 0 ):
        if( options.verbose ):
            # log_line += "\n  %d "%()
            log_line  = "\n  Date %s "%(str(t_i))
            log_line += "\n  Number of processors  %d "%(size)
            log_line += "\n  Particles of group i: %d "%(sum_i)
            log_line += sperator_line
            log_line += "\n id# ; type ; ff type ; ring # ; residue name ; residue # ; link id  "
            log_line += sperator_line
            #for p_i, ptcl in ptclC_o(list_i):
            #    log_line += "\n   %d %s %s %s %s %s %s "%(p_i,ptcl.type,ptcl.tagsDict["fftype"],ptcl.tagsDict["ring"],ptcl.tagsDict["resname"],ptcl.tagsDict["residue"],ptcl.tagsDict["linkid"])
            log_line += "\n  Particles of group i: %d "%(sum_j)
            log_line += sperator_line
            log_line += "\n id# ; type ; ff type ; ring # ; residue name ; residue # ; link id  "
            log_line += sperator_line
            #for p_j, ptcl in ptclC_o(list_j):
            #    log_line += "\n   %d %s %s %s %s %s %s "%(p_j,ptcl.type,ptcl.tagsDict["fftype"],ptcl.tagsDict["ring"],ptcl.tagsDict["resname"],ptcl.tagsDict["residue"],ptcl.tagsDict["linkid"])

            log_line += sperator_line
            log_line += "\n Frames "
            log_line += "\n     Initial frame  %d "%(options.frame_o)
            log_line += "\n     Step frame  %d  "%(options.frame_step)
            log_line += "\n     Final frame  %d  "%(options.frame_f)

                
            log_out.write(log_line)
            print log_line
            
	       
    #
    # If multi-core split the number of atoms in the molecule onto each core
    #
    debug = 0 
    if( debug ): print rank, size," splitOnProcs "
    split_list = True
    # Create a list of atomic indices for each processor 
    myChunk_i  = p.splitListOnProcs(list_i)
    debug = 0
    if(debug):                
        print " cpu ",rank ," has atoms ",myChunk_i[0]," - ",myChunk_i[len(myChunk_i)-1],"  \n"
	for atom_i in myChunk_i:
	    print "Processor %d has atom %d "%(rank,atom_i)
	sys.exit(" debug myChunk ")
	
    p.barrier()
         
    # Calculate rdf relate values
    n_bins = int(options.r_cut/options.bin_size)
    dmin = 0.0
    dmax = options.r_cut
    sq_r_cut = options.r_cut**2
    rdf_cnt_ij = np.zeros(n_bins+1)    
    volume_i = []
    volume_sum_i = 0.0 
    rdf_frames = 0
    proc_frame = 0

    # set up rdf
    rdf, edges = numpy.histogram([0], bins=n_bins, range=(dmin, dmax))
    rdf *= 0
    rdf = rdf.astype(numpy.float64)  # avoid possible problems with '/' later on

    if( len(options.in_lammpsxyz) ):
        
        if( rank == 0 ):
            if( options.verbose ): print "  Reading ",options.in_lammpsxyz
            rdft_i = datetime.datetime.now()
            
        # Initialize line count
        NP = len( ptclC_o )
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
                    struc_o.ptclC[pt_cnt].position = np.array(  [ float(col[1]),float(col[2]),float(col[3]) ] )
                    atomic_symb = col[0]
                    
                    #if( pt_cnt == 1 ):
                        #    print struc_o.ptclC[pt_cnt].postion
                        # print "  POS ",pt_cnt, struc_o.ptclC[pt_cnt].postion[0], col[1]
                
                    #
                    #r_i = [ col[1],col[2],col[3] ]
                    #part_i = Particle( r_i,atomic_symb,q_i,m_i )
                    #struc_i.put(part_i)

                if( line_cnt > 1  and line_cnt > frame_cnt*(NP + 2) ):
                    pt_cnt = 0

                    if( options.frame_o <= frame_cnt ):
                        if( frame_cnt <= options.frame_f  or  options.readall_f  ):
                            if( frame_cnt%options.frame_step == 0 ):

                                rdf_frames += 1
                                proc_frame += 1
                                
                                debug = 0
                                if( proc_frame > (size) ):
                                    proc_frame = 1
                                    if( debug ): print rank," proc frame > size ",size +1," reset to ",proc_frame

                                if( split_list ): #proc_frame == (rank+1) ):

                                    if( debug ): print rank," calcuating rdf frame ",rdf_frames,"  proc frame ",proc_frame

                                    struc_o.maxminLatVec()
                                    vol_i = struc_o.getVolume()
                                    volume_i.append( vol_i )
                                    volume_sum_i +=  vol_i

                                    debug = 1
                                    if( debug ):
                                        print "proc",rank+1," analysis of frame %d on processor %d "%(frame_cnt,rank+1)," box from min max %f "%struc_o.latvec[0][0]

                                    #rdf_cnt_ij = struc_o.calc_rdf(rdf_cnt_ij,options.bin_size,list_i,list_j,sq_r_cut)
                                    #
                                    # Loop over list i
                                    #
                                    # for p_i, ptcl_i in struc_o.ptclC(list_i):

                                    print  "       Startin loop "
                                        
                                    for p_i in myChunk_i:
                                        ptcl_i = struc_o.ptclC[p_i]
                                        r_i = np.array( [float(ptcl_i.position[0]),float(ptcl_i.position[1]),float(ptcl_i.position[2] )] )

                                        
                                        if( prtCl_time ):
                                            pt_i = datetime.datetime.now()
                                            
                                        #
                                        # Loop over list j
                                        #
                                        for p_j, ptcl_j in struc_o.ptclC(list_j):
                                            add_ij = True 
                                            if( p_j <= p_i):
                                                add_ij = False 
                                            #
                                            # Check for intra vs. inter molecular
                                            #
                                            if( options.mol_inter ):
                                                if( ptcl_i.tagsDict["chain"]  == ptcl_j.tagsDict["chain"] ): add_ij = False
                                            if( options.mol_intra ):
                                                if( ptcl_i.tagsDict["chain"]  != ptcl_j.tagsDict["chain"] ): add_ij = False
                                            if( add_ij ):
                                                r_j =  np.array( [float(ptcl_j.position[0]),float(ptcl_j.position[1]),float(ptcl_j.position[2])] )
                                                r_ij_sq = pbcs.sq_drij_c(r_i,r_j,struc_o.latvec)
                                                if( r_ij_sq <= sq_r_cut ):
                                                    m_ij = np.sqrt(r_ij_sq)
                                                    bin_index = int( round( m_ij/options.bin_size) )
                                                    rdf_cnt_ij[bin_index] += 2

                                        if( prtCl_time ):
                                            pt_f = datetime.datetime.now()

                                            dt_sec  = pt_f.second - pt_i.second
                                            dt_min  = pt_f.minute - pt_i.minute
                                            if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                            if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                            print  "       Particle  %d  rdf calculated in %d min  %f sec "%(p_i,dt_min,dt_sec)                        
                        
                                # Print output
                                if( rank == 0 and options.verbose ):
                                    rdft_f = datetime.datetime.now()
                                    dt_sec  = rdft_f.second - rdft_i.second
                                    dt_min  = rdft_f.minute - rdft_i.minute
                                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                    if ( dt_min < 0 ): dt_min = 60.0 - dt_min
                                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                    log_line = "\n       Frame  %d  rdf calculated in %d min  %f sec "%(frame_cnt,dt_min,dt_sec)
                                    log_line  += "\n       with volume %f Angstroms^3 estimated from max/min"%(vol_i)
                                    log_out.write(log_line)
                                    print log_line
                                    # Rest intial time for next rdf calcution
                                    rdft_i = datetime.datetime.now()

                    frame_cnt += 1
                
            # Calculate rdf for last frame

                    

            if( options.frame_o <= frame_cnt ):

                if( frame_cnt <= options.frame_f  or  options.readall_f  ):
                    if( frame_cnt%options.frame_step == 0 ):
                        # Index counts 
                        rdf_frames += 1
                        proc_frame += 1

                        debug = 1
                        if( debug ): print " calcuating rdf frame ",rdf_frames
                        if( debug ): print "  proc frame ",proc_frame
                        if( proc_frame > (size  ) ):
                            proc_frame = 1
                            if( debug ): print " proc frame > size ",size +1," reset to ",proc_frame

                        if( split_list ): #proc_frame == (rank+1) ):

                            print " analysis of frame %d on processor %d "%(frame_cnt,rank+1)

                            struc_o.maxminLatVec()                            
                            vol_i = struc_o.getVolume()
                            volume_i.append( vol_i )
                            volume_sum_i +=  vol_i

                            print " box from min max %f "%struc_o.latvec[0][0]

                            #rdf_cnt_ij = struc_o.calc_rdf(rdf_cnt_ij,options.bin_size,list_i,list_j,sq_r_cut)
                            #
                            # Loop over list i
                            #
                            # for p_i, ptcl_i in struc_o.ptclC(list_i):
                            for p_i in myChunk_i:
                                ptcl_i = struc_o.ptclC[p_i]
                                #r_i = np.array( [float(ptcl_i.position[0]),float(ptcl_i.position[1]),float(ptcl_i.position[2] )] )
                                r_i = ptcl_i.position # np.array( [float(ptcl_i.position[0]),float(ptcl_i.position[1]),float(ptcl_i.position[2] )] )

                                
                                if( prtCl_time ):
                                    pt_i = datetime.datetime.now()                                            
                                #
                                # Loop over list j
                                #
                                for p_j, ptcl_j in struc_o.ptclC(list_j):
                                    add_ij = True
                                    if( p_j <= p_i):
                                        add_ij = False 
                                    #
                                    # Check for intra vs. inter molecular
                                    #
                                    if( options.mol_inter ):
                                        if( ptcl_i.tagsDict["chain"]  == ptcl_j.tagsDict["chain"] ): add_ij = False
                                    if( options.mol_intra ):
                                        if( ptcl_i.tagsDict["chain"]  != ptcl_j.tagsDict["chain"] ): add_ij = False
                                    if( add_ij ):
                                        #r_j =  np.array( [float(ptcl_j.position[0]),float(ptcl_j.position[1]),float(ptcl_j.position[2])] )
                                        r_j =  ptcl_j.position #np.array( [float(ptcl_j.position[0]),float(ptcl_j.position[1]),float(ptcl_j.position[2])] )
                                        
                                        r_ij_sq = pbcs.sq_drij_c(r_i,r_j,struc_o.latvec)
                                        if( r_ij_sq <= sq_r_cut ):
                                            m_ij = np.sqrt(r_ij_sq)
                                            bin_index = int( round( m_ij/options.bin_size) )
                                            rdf_cnt_ij[bin_index] += 2


                                if( prtCl_time ):
                                    pt_f = datetime.datetime.now()
                                    
                                    dt_sec  = pt_f.second - pt_i.second
                                    dt_min  = pt_f.minute - pt_i.minute
                                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                    print  "       Particle  %d  rdf calculated in %d min  %f sec "%(p_i,dt_min,dt_sec)                        
                                                                    
                        
                            # Print output
                            if( rank == 0 and options.verbose ):
                                rdft_f = datetime.datetime.now()
                                dt_sec  = rdft_f.second - rdft_i.second
                                dt_min  = rdft_f.minute - rdft_i.minute
                                if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                if ( dt_min < 0 ): dt_min = 60.0 - dt_min
                                if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                log_line = "\n       Frame  %d  rdf calculated in %d min  %f sec "%(frame_cnt,dt_min,dt_sec)
                                log_line  += "\n       with volume %f Angstroms^3 estimated from max/min"%(vol_i)
                                log_out.write(log_line)
                                print log_line
                                # Rest intial time for next rdf calcution
                                rdft_i = datetime.datetime.now()


    if( debug):

        # Read in using mdanalysis
        universe =  Universe(options.in_gro)
        un =  universe.selectAtoms("resname * ")
        coord = un.coordinates()

        print "    coord[0] ",   coord[0]


        # Relabel segments to correspond to lists i and j
        for pid_i in list_i:
            print pid_i,ptclC_o[pid_i].tagsDict["gtype"],ptclC_o[pid_i].position," from u ",universe.atoms[pid_i-1].name #, universe.coord[pid_i-1]
            universe.atoms[pid_i-1].resname = "grpi"
        uni_i = universe.selectAtoms(" resname grpi ")
        coor_i = uni_i.coordinates()

        print "    coor_i[0] ",   coor_i[0]


        sys.exit("gro testing ")

    if( len(options.in_xtc ) and  len(options.in_gro) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s and %s "%(options.in_xtc,options.in_gro)
            rdft_i = datetime.datetime.now()
        # Read in using mdanalysis
        universe =  Universe(options.in_gro, options.in_xtc)
        un =  universe.selectAtoms("resname * ")
        coord = un.coordinates()
        #print " len(coord)",len(coord)

        
        # Relabel segments to correspond to lists i and j
        for pid_i in list_i:
            print pid_i,ptclC_o[pid_i].tagsDict["gtype"],ptclC_o[pid_i].position," from u ",universe.atoms[pid_i-1].name , universe.coord[pid_i-1]
            universe.atoms[pid_i-1].resname = "grpi"
            universe.atoms[pid_i-1].resnum = ptclC_o[pid_i].tagsDict["chain"]  # Set resnum to chain number for inter/intra rdf's
        uni_i = universe.selectAtoms(" resname grpi ")

        for pid_j in list_j:
            print pid_j,ptclC_o[pid_j].tagsDict["gtype"]," from u ",universe.atoms[pid_j-1].name
            universe.atoms[pid_j-1].resname = "grpj"
            universe.atoms[pid_j-1].resnum = ptclC_o[pid_j].tagsDict["chain"]  # Set resnum to chain number for inter/intra rdf's
        uni_j = universe.selectAtoms(" resname grpj ")
        
        n_i = uni_i.numberOfAtoms()
        n_j = uni_j.numberOfAtoms()

        coor_i = uni_i.coordinates()
        coor_j = uni_j.coordinates()

        if( rank == 0 and options.verbose ):
            print " mdanalysis groups "
            print "    uni_i  %d "%n_i
            p_i = 0
            p_j = 0
            for a_i in  uni_i.atoms:
                p_i += 1
                print p_i,a_i.number,a_i.name,a_i.type,a_i.resname,a_i.resid,a_i.mass,coord[a_i.number]
            print "    uni_j  %d "%n_j
            for a_j in  uni_j.atoms:
                print a_j.number,a_j.name,a_j.type,a_j.resname,a_j.resid,a_j.mass,coord[a_j.number]


        # Allocate distance matrix 
        dist = numpy.zeros((n_i,n_j), dtype=numpy.float64)


        # Create include matrix of pairs to include
        #include_ij = [[None]*n_j]*n_i
        # had to use numpy list since regular list was asigning values include_ij[:][p_j] instead of include_ij[p_i][p_j]
        # not sure why 
        include_ij =  numpy.zeros((n_i,n_j), dtype=numpy.int)

        p_i = -1
        
        # Relabel segments to correspond to lists i and j
        for a_i in  uni_i.atoms:
            p_i += 1
            p_j = -1
            for a_j in  uni_j.atoms:
                p_j += 1

                #print p_i,p_j
                #print include_ij

                add_ij = 1
                if( a_j.number <= a_i.number ):
                    add_ij = 0
                #
                # Check for intra vs. inter molecular
                #
                if( options.mol_inter ):
                    if( a_i.resnum  == a_j.resnum ):
                        add_ij = 0
                if( options.mol_intra ):
                    if(  a_i.resnum  != a_j.resnum ):
                        add_ij = 0

                include_ij[p_i][p_j] = add_ij

        debug = False 
        if( debug):
            print " len(include_ij) ", len(include_ij)
            print " len(include_ij[0]) ", len(include_ij[0])
            print "include_ij[0][1]",include_ij[0][1]
            print include_ij[n_i-1][n_j-1]

            sys.exit(" include array size 0 ")
            for p_i in range(n_i):
                for p_j in range(n_j):
                    
                    print " for pair  ",include_ij[p_i][p_j],p_i,p_j,dist[p_i,p_j] 
                                    
            sys.exit(" include array size ")

        volume_sum_i = 0
        rdf_frames = 0 

        # setting the skip for the traj
        # does not seem to work 
        #universe.trajectory.skip = options.frame_step
        
        for ts in universe.trajectory:
            if( options.frame_o <= ts.frame ):
                if( ts.frame <= options.frame_f  or  options.readall_f  ):
                    if( ts.frame%options.frame_step == 0 ):

                        print "Frame %4d with volume %f " % (ts.frame, ts.volume)
                        rdf_frames += 1 
                        volume_sum_i += ts.volume      # correct unitcell volume
                        box = ts.dimensions
                        coor_i = uni_i.coordinates()
                        coor_j = uni_j.coordinates()            
                        distance_array(coor_i,coor_j, box, result=dist)  # use pre-allocated array, box not fully correct!!

                        for p_i in range(n_i):
                            for p_j in range(n_j):

                                #if( debug):
                                #    print " checking ",include_ij[p_i][p_j],p_i,p_j,dist[p_i,p_j] 
                                    
                                if( include_ij[p_i][p_j] == 1 ):
                                    if( dist[p_i,p_j] <= options.r_cut ):
                                        bin_index = int( round( dist[p_i,p_j] / options.bin_size) )
                                        rdf_cnt_ij[bin_index] += 2
                                else:
                                    # Set unincluded pairs to have a seperation beyond the cut off
                                     dist[p_i,p_j] = options.r_cut + 1.0 

                        new_rdf, edges = numpy.histogram(dist, bins=n_bins, range=(dmin, dmax))
                        rdf += new_rdf

        #numframes = universe.trajectory.numframes / universe.trajectory.skip

        volume_sum_i /= rdf_frames    # average volume

        # Normalize RDF
        radii = 0.5*(edges[1:] + edges[:-1])
        vol = (4./3.)*numpy.pi*(numpy.power(edges[1:],3)-numpy.power(edges[:-1], 3))
        # normalization to the average density n/volume_sum_i in the simulation
        density = n_i / volume_sum_i
        # This is inaccurate when solutes take up substantial amount
        # of space. In this case you might want to use
        ## import MDAnalysis.core.units
        ## density = MDAnalysis.core.units.convert(1.0, 'water', 'Angstrom^{-3}')
        norm = density * n_i*n_j * rdf_frames
        rdf /= norm * vol

        outfile = './rdf_mda.dat'
        with open(outfile,'w') as output:
            for radius,gofr in izip(radii, rdf):
                output.write("%(radius)8.3f \t %(gofr)8.3f\n" % vars())
            p.barrier() # Barrier for MPI_COMM_WORLD

    
            
    if( options.verbose and rank == 0 ):
	print "      Wating for all processors to finish  "                                
    p.barrier()			
    if( options.verbose and rank == 0 ):
	print "      Finding averages "
    #
    # Find averages
    #
    debug = 0 
    rdf_cnt = np.zeros(n_bins+1)   
    for bin_index in range( n_bins+1):
	# Sum rdf_cnt of each bin on each processor 
	rdf_cnt[bin_index] = p.allReduceSum(rdf_cnt_ij[bin_index])

        # print " Proc %d rdf_cnt_ij %f sum %f "% (rank,rdf_cnt_ij[bin_index],rdf_cnt[bin_index] )


    
    p.barrier() # Barrier for MPI_COMM_WORLD

    
    #
    # Calculate rdf results 
    # 
    total_cnts = np.sum( rdf_cnt )
    volume_sum = volume_sum_i #p.allReduceSum(volume_sum_i)
    box_vol_ave = volume_sum/float(rdf_frames)
    vol_cut = 4.0*math.pi/3.0*options.r_cut**3
    n_shperes = float(sum_i)*float(rdf_frames)
    sphere_den_j = float(total_cnts)/vol_cut/n_shperes #/2.0  # N_B A^-3
    box_den_i = float(sum_i )/float(box_vol_ave)
    box_den_j = float(sum_j )/float(box_vol_ave)
    
    if( rank == 0 ):
	if( options.verbose ):
            
		
	    print "   Frames ",rdf_frames
	    print "   Total counts ",total_cnts
	    print "   Average box volume ",box_vol_ave
	    print "   Volume of cut-off sphere ",vol_cut
	    print "   Average box density i ",box_den_i," atoms/angstrom^3 "
	    print "   Average box density j ",box_den_j," atoms/angstrom^3 "
	    print "   Average cut-off sphere density ",sphere_den_j," atoms/angstrom^3 "
		
	# Write output 
	#
	dat_out.write("# RDF frames %d %d " %  (options.frame_o,options.frame_f))
	dat_out.write("\n#    Bin-size %f  " % (options.bin_size))
	dat_out.write("\n#    Cut-off %f  " % (options.r_cut))
	dat_out.write("\n#    Frames %d  " % (rdf_frames))
	dat_out.write("\n#    Total_cnts %d  " % (total_cnts))
	dat_out.write("\n#    N_i %d " % (sum_i ))
	dat_out.write("\n#    N_j %d " % (sum_j ))
	dat_out.write("\n#    Average Box Volume %f " % ( box_vol_ave) )
	dat_out.write("\n#    Box density i %f N A^-3 " % (box_den_i ))
	dat_out.write("\n#    Box density j %f N A^-3 " % (box_den_j ))
	dat_out.write("\n#    Sphere volume  %f A^3 " % (vol_cut ))
	dat_out.write("\n#    Average Sphere density  %f N A^3 " % (sphere_den_j ))
	dat_out.write("\n#    ")
	dat_out.write("\n# bin index ; r     ; count_ave/frame ; dr vol ;  dr vol(aprox) ; g_sphere ; g_boxs  ")
	#                bin_index , r_val , dr_cnt_norm      , dr_vol,  dr_vol_apx,     sphere_g, box_g
	
	for bin_index in range( 1,n_bins):
		
	    r_val = options.bin_size*float(bin_index)
	    r_in = r_val - options.bin_size*0.5
	    r_out = r_val + options.bin_size*0.5
    
	    dr_cnt = float( rdf_cnt[bin_index] )
	    dr_cnt_norm =     dr_cnt    /float(rdf_frames)
	    dr_vol = 4.0*math.pi/3.0*( r_out**3 - r_in**3 )
	    dr_vol_apx = 4.0*math.pi*(  r_val**2 )
	    
	    dr_rho = dr_cnt_norm/dr_vol
	    sphere_g = dr_rho/sphere_den_j/float( sum_i )
	    box_g = dr_rho/box_den_j/float( sum_i )
	    
	    dat_out.write("\n  %d %f %f %f %f %f %f " % (bin_index,r_val,dr_cnt_norm,dr_vol,dr_vol_apx,sphere_g,box_g) )
	    

        
        t_f = datetime.datetime.now()
        dt_sec  = t_f.second - t_i.second
        dt_min  = t_f.minute - t_i.minute
        if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
        if ( dt_min < 0 ): dt_min = 60.0 - dt_min
        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
        log_line="\n  Finished time  " + str(t_f)
        log_line+="\n  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "
        if( options.verbose ):
            print log_line
        log_out.write(log_line)

	dat_out.close()
	log_out.close()
	    

if __name__=="__main__":
    main()
   
