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
import json, math 
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
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    
        
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
    
    #
    # Formated ouput varables
    #
    sperator_line = "\n---------------------------------------------------------------------"
    
    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()



    if( rank == 0 ):
        # record initial time 
        t_i = datetime.datetime.now()
            
    options, args = get_options()


    # Read in system data from json file
    if( len(options.in_json)):
        struc_o = StructureContainer()
        struc_o.getsys_json(options.in_json)

        # Get paticle and bond structures
        ptclC_o = struc_o.ptclC
        bondC_o  = struc_o.bondC
    
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
            for p_i, ptcl in ptclC_o(list_i):
                log_line += "\n   %d %s %s %s %s %s %s "%(p_i,ptcl.type,ptcl.tagsDict["fftype"],ptcl.tagsDict["ring"],ptcl.tagsDict["resname"],ptcl.tagsDict["residue"],ptcl.tagsDict["linkid"])
            log_line += "\n  Particles of group i: %d "%(sum_j)
            log_line += sperator_line
            log_line += "\n id# ; type ; ff type ; ring # ; residue name ; residue # ; link id  "
            log_line += sperator_line
            for p_j, ptcl in ptclC_o(list_j):
                log_line += "\n   %d %s %s %s %s %s %s "%(p_j,ptcl.type,ptcl.tagsDict["fftype"],ptcl.tagsDict["ring"],ptcl.tagsDict["resname"],ptcl.tagsDict["residue"],ptcl.tagsDict["linkid"])

            log_line += sperator_line
            log_line += "\n Frames "
            log_line += "\n     Initial frame  %d "%(options.frame_o)
            log_line += "\n     Step frame  %d  "%(options.frame_step)
            log_line += "\n     Final frame  %d  "%(options.frame_f)

                
            log_out.write(log_line)
            print log_line
            
    # Calculate rdf relate values
    n_bins = int(options.r_cut/options.bin_size)
    sq_r_cut = options.r_cut**2
    rdf_cnt_ij = np.zeros(n_bins+1)    
    volume_i = []
    rdf_frames = 0 
    # Record initial time

    if( rank == 0  ): 
	t_i = datetime.datetime.now()

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
                    pt_i = ptclC_o[pt_cnt]
                    pt_i.postion =  [ col[1],col[2],col[3] ]
                    atomic_symb = col[0]
                    #
                    #r_i = [ col[1],col[2],col[3] ]
                    #part_i = Particle( r_i,atomic_symb,q_i,m_i )
                    #struc_i.put(part_i)

                if( line_cnt > 1  and line_cnt > frame_cnt*(NP + 2) ):
                    pt_cnt = 0

                    if( options.frame_o <= frame_cnt and frame_cnt <= options.frame_f ):

                        if( frame_cnt%options.frame_step == 0 ):
                            rdf_frames += 1
                            rdf_cnt_ij = struc_o.calc_rdf(rdf_cnt_ij,options.bin_size,list_i,list_j,sq_r_cut)
                            volume_i.append( struc_o.getVolume() )
                            
                            # Print output
                            if( rank == 0 and options.verbose ):
                                rdft_f = datetime.datetime.now()
                                dt_sec  = rdft_f.second - rdft_i.second
                                dt_min  = rdft_f.minute - rdft_i.minute
                                if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                log_line = "\n       Frame  %d  rdf calculated in %d min  %f sec "%(frame_cnt,dt_min,dt_sec)
                                log_out.write(log_line)
                                print log_line
                                # Rest intial time for next rdf calcution
                                rdft_i = datetime.datetime.now()

                    frame_cnt += 1
                
            # Calculate rdf for last frame 
            if( options.frame_o <= frame_cnt and frame_cnt <= options.frame_f ):
                if( frame_cnt%options.frame_step  == 0 ):
                    rdf_frames += 1
                    rdf_cnt_ij = struc_o.calc_rdf(rdf_cnt_ij,options.bin_size,list_i,list_j,sq_r_cut)
                    volume_i.append( struc_o.getLatVec() )
            
                    if( rank == 0 and options.verbose ):
                        rdft_f = datetime.datetime.now()
                        dt_sec  = rdft_f.second - rdft_i.second
                        dt_min  = rdft_f.minute - rdft_i.minute
                        if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                        log_line = "       Frame  %d  rdf calculated in %d min  %f sec \n"%(frame_cnt,dt_min,dt_sec)
                        log_out.write(log_line)
                        print log_line
                        # Rest intial time for next rdf calcution
                        rdft_i = datetime.datetime.now()



    
    total_cnts = np.sum( rdf_cnt_ij)

    box_vol_ave = np.average( volume_i )
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


    
	    dr_cnt = float( rdf_cnt_ij[bin_index] )
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
        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
        log_line="\n  Finished time  " + str(t_f)
        log_out.write(log_line)
        log_line="\n  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "
        log_out.write(log_line)

	dat_out.close()
	log_out.close()
	    

if __name__=="__main__":
    main()
   