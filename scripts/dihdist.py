#! /usr/bin/env python
"""
Find distribution of dihedral angles

k-i-j-l

"""

# Dr. Travis Kemper
# Initial Date July 2014
# travis.kemper@nrel.gov

from structureContainer import StructureContainer
import mpiNREL
import datetime, sys
import numpy as np

const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    # json files to act on
    parser.add_option("-j","--in_json", dest="in_json", type="string", default="", help="Input json file, which read in first then over writen by subsequent input files")
    parser.add_option("-o","--output_id", dest="output_id", default="trans",type="string",help=" prefix for output files  ")

    parser.add_option("--bin_size", dest="bin_size", type=float, default=2.00, help=" Bin size in degrees ")
    
    # Searchable properties
    # k
    parser.add_option("--id_k", dest="id_k", type="string", default="", help=" select atoms of group k by id number  ")    
    parser.add_option("--symb_k", dest="symb_k", type="string", default="", help=" select atoms of group k by (atomic) symbol   ")    
    parser.add_option("--chains_k", dest="chains_k", type="string", default="", help="select atoms of group k by chain number  ")    
    parser.add_option("--ring_k", dest="ring_k", type="string", default="", help="select atoms of group k by particles in a ring   ")    
    parser.add_option("--resname_k", dest="resname_k", type="string", default="", help="select atoms of group k by residue name  ")    
    parser.add_option("--residue_k", dest="residue_k", type="string", default="", help="select atoms of group k by resudue number  ")    
    parser.add_option("--linkid_k", dest="linkid_k", type="string", default="", help="select atoms of group k by  link id ")    
    parser.add_option("--fftype_k", dest="fftype_k", type="string", default="", help="select atoms of group k by force field type  ")
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
    # l
    parser.add_option("--id_l", dest="id_l", type="string", default="", help=" select atoms of group l by id number  ")    
    parser.add_option("--symb_l", dest="symb_l", type="string", default="", help=" select atoms of group l by (atomic) symbol   ")    
    parser.add_option("--chains_l", dest="chains_l", type="string", default="", help="select atoms of group l by chain number  ")    
    parser.add_option("--ring_l", dest="ring_l", type="string", default="", help="select atoms of group l by particles in a ring   ")    
    parser.add_option("--resname_l", dest="resname_l", type="string", default="", help="select atoms of group l by residue name  ")    
    parser.add_option("--residue_l", dest="residue_l", type="string", default="", help="select atoms of group l by resudue number  ")    
    parser.add_option("--linkid_l", dest="linkid_l", type="string", default="", help="select atoms of group l by  link id ")    
    parser.add_option("--fftype_l", dest="fftype_l", type="string", default="", help="select atoms of group l by force field type  ")

    # Trajectory file
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")
    #
    # Filters
    #
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=False,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")

    #
    # Modify input system file 
    #
    
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

def dih_hist(angle_list,struc_o,bin_size,hist_cnt):
    """
    Loop over dihedrals and calculate histogram 
    """
    
    debug = False
    
    #
    dih_cnt = 0

    for indx_kij in angle_list:
        dih_cnt += 1 
        a_k = indx_kij[0]
        a_i = indx_kij[1]
        a_j = indx_kij[2]
        a_l = indx_kij[3]

        angle_i = struc_o.getDihedral(a_k,a_i,a_j,a_l)
        abs_angle_i = np.absolute(angle_i)

        #if( abs_angle_i > 180.0 - bin_size/2.0 ):
        #    if( debug ):
        #        print " rescaling angle ",abs_angle_i," to ", abs_angle_i  -180.0 
        #    abs_angle_i += -180.0
        if( abs_angle_i <= 180.0 ):
            # Acount for round off errors 
            abs_angle_i += -0.0001
            
        bin_index = int(  abs_angle_i/bin_size ) 
        hist_cnt[bin_index] += 1

        
        if( debug ):
            hist_val = bin_size*float(bin_index) + bin_size/2.0
            if( bin_index == 0 ): hist_val = 0.0 
            if( hist_val == 179.0 ): hist_val = 180.0 
            print " index %d or index %d or index %d or index %d  %f bin_index %d hist_cnt %d hist_val %f "%(a_k-1,a_i-1,a_j-1,a_l-1,angle_i,bin_index, hist_cnt[bin_index],hist_val)

    return hist_cnt

def dihdist():
    """
    Read in files and create new files

    Arguments

    Return
    
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
        
    #        
    # Read options 
    #
    options, args = get_options()
    #
    #  Initialize blank system 
    # 
    struc_o = StructureContainer()
    #
    # Read in json file
    #
    if( len(options.in_json) ):
        if( rank == 0 and options.verbose ):
            print  "     - Reading in ",options.in_json
        json_data_i = struc_o.getsys_json(options.in_json)


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

    #   
    # Create sub lists for RDF groups k-i-j-l
    #
    if( rank == 0 and options.verbose  ):
        print "  searching for list k i j l "
        
    search_k = create_search(options.id_k,options.symb_k,options.chains_k,options.ring_k,options.resname_k,options.residue_k,options.linkid_k,options.fftype_k)
    if( rank == 0 ):
        if( options.verbose ): print " Searching group k ",search_k
    list_k = ptclC_o.getParticlesWithTags(search_k)
    sum_k = len(list_k)
    
    
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
    
    
    search_l = create_search(options.id_l,options.symb_l,options.chains_l,options.ring_l,options.resname_l,options.residue_l,options.linkid_l,options.fftype_l)
    if( rank == 0 ):
        if( options.verbose ): print " Searching group l ",search_l
    list_l = ptclC_o.getParticlesWithTags(search_l)
    sum_l = len(list_l)

    if( rank == 0 and options.verbose  ):
        print "  searching finished "
        print "  Creating list of atomic indies of dihedrals "
    
    
    if( rank == 0 ):
        if( options.verbose ): print " search finished "
    

    

    # Create list of dihderal angles 
    angle_list = struc_o.get_dihatoms(list_k,list_i,list_j,list_l)
    
    
    if( rank == 0 and options.verbose  ):
        print "  list creation finished "
    
    if( rank == 0 ):
        if( options.verbose ):
            # log_line += "\n  %d "%()
            log_line  = "\n  Date %s "%(str(t_i))
            log_line += "\n  Number of processors  %d "%(size)
            
            log_line += sperator_line
            log_line += "\n  Particles of group l: %d "%(sum_l)
            log_line += "\n  Particles of group i: %d "%(sum_i)
            log_line += "\n  Particles of group j: %d "%(sum_j)
            log_line += "\n  Particles of group k: %d "%(sum_k)

            log_line += sperator_line
            log_line += "\n Frames "
            log_line += "\n     Initial frame  %d "%(options.frame_o)
            log_line += "\n     Step frame  %d  "%(options.frame_step)
            log_line += "\n     Final frame  %d  "%(options.frame_f)

            log_out.write(log_line)
            print log_line
            
    p.barrier()

    #
    # Intialize lsits and counts
    #
    angle_max = 180.0
    n_bins = int(angle_max/options.bin_size)
    hist_cnt = np.zeros(n_bins+1)
    #
    # Initialize line count
    #
    frame_cnt = 0
    calc_frames = 0 
    volume_i = []
    volume_sum_i = 0.0 

    if( len(options.in_lammpsxyz) ):
        
        if( rank == 0 ):
            if( options.verbose ): print "  Reading ",options.in_lammpsxyz
            calct_i = datetime.datetime.now()
            
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
                    struc_o.ptclC[pt_cnt].position =  [ float(col[1]),float(col[2]),float(col[3]) ]
                    atomic_symb = col[0]
                    
                if( line_cnt > 1  and line_cnt > frame_cnt*(NP + 2) ):
                    pt_cnt = 0

                    if( options.frame_o <= frame_cnt ):
                        if( frame_cnt <= options.frame_f  or  options.readall_f  ):
                            if( frame_cnt%options.frame_step == 0 ):

                                calc_frames += 1

                                struc_o.maxminLatVec()
                                vol_i = struc_o.getVolume()
                                volume_i.append( vol_i )
                                volume_sum_i +=  vol_i

                                debug = 1
                                if( debug ):
                                    print "proc",rank+1," analysis of frame %d on processor %d "%(frame_cnt,rank+1)," box from min max %f "%struc_o.latvec [0][0]
                                #
                                # Loop over lists
                                #
                                hist_cnt = dih_hist(angle_list,struc_o,options.bin_size,hist_cnt)

                    frame_cnt += 1
                
            # Calculate 

            if( options.frame_o <= frame_cnt ):

                if( frame_cnt <= options.frame_f  or  options.readall_f  ):
                    if( frame_cnt%options.frame_step == 0 ):
                        # Index counts 
                        calc_frames += 1

                        struc_o.maxminLatVec()
                        vol_i = struc_o.getVolume()
                        volume_i.append( vol_i )
                        volume_sum_i +=  vol_i

                        debug = 1
                        if( debug ):
                            print "proc",rank+1," analysis of frame %d on processor %d "%(frame_cnt,rank+1)," box from min max %f "%struc_o.latvec[0][0]

                        #
                        # Loop over lists
                        #
                        hist_cnt = dih_hist(angle_list,struc_o,options.bin_size,hist_cnt)

    if( rank == 0 ):
        hist_sum = np.sum(hist_cnt)
        box_vol_ave = np.average( volume_i )
    

        # Write output 
        #
        dat_line = "#    Frames %d " %  (calc_frames)
        dat_line += "\n#    Bin-size %f  " % (options.bin_size)
        dat_line += "\n#    Total_cnts %d  " % (hist_sum)
        dat_line += "\n#    Average Box Volume %f " % ( box_vol_ave) 
        dat_line += "\n# "
        dat_line += "\n# bin index ; cnt    ; cnt/Total_cnts  "
        dat_out.write(dat_line)
        
	if( options.verbose ):
	    print dat_line

        if( debug ):
            sum_check = 0 
        for bin_index in range( 0,n_bins+1):
            hist_val = options.bin_size*float(bin_index) + options.bin_size/2.0
            if( bin_index == 0 ): hist_val = 0.0
            if( hist_val == 179.0 ): hist_val = 180.0 
                        
            #print " hist_cnt[bin_index] ",bin_index, hist_cnt[bin_index]
            
            val_cnt =  hist_cnt[bin_index] 
            cnt_fnorm =   float(val_cnt)/float(hist_sum) #float(calc_frames)

            
            dat_out.write("\n  %d %f %f %f  " % (bin_index,hist_val,val_cnt,cnt_fnorm) )

            if( debug ): sum_check += val_cnt

        if( debug ): print "sum_check ",sum_check
        
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
        
        log_out.close()
        dat_out.close()
    
	
	    
if __name__=="__main__":
    dihdist()
   
	    
