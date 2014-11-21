#! /usr/bin/env python
"""
Find average properties of a simulation 

"""

# Dr. Travis Kemper
# Initial Date July 2014
# travis.kemper@nrel.gov

import datetime, sys
import numpy as np

from structureContainer import StructureContainer
from parameters import ParameterContainer

import mpiNREL,file_io
import particles 


def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-o","--output_id", dest="output_id", default="trans",type="string",help=" prefix for output files  ")
    #
    # Input files 
    #
    parser.add_option("-j","--in_json", dest="in_json", type="string", default="", help="Input json file, which read in first then over writen by subsequent input files")
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--in_itp", dest="in_itp",  type="string", default="",help="Input gromacs force field parameter file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_xmol", dest="in_xmol", type="string", default="", help="Input xmol file")
    parser.add_option("--xmol_format", dest="xmol_format", type="string", default="atomic_symb", help="Format of xmol file: atomic_symb - first colum atomic symbols; lammps - first colum particle type atomic symbols will be set to carbon; ")
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")
    #
    # Filters
    #
    parser.add_option("--id", dest="id", type="string", default="", help=" select atoms of group by number  ")    
    parser.add_option("--symbol", dest="symbol", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--type", dest="type", type="string", default="", help=" select type  ")    
    parser.add_option("--chains", dest="chains", type="string", default="", help="select atoms of group by chain number  ")    
    parser.add_option("--ring", dest="ring", type="string", default="", help="select atoms of group by particlesn a ring   ")    
    parser.add_option("--resname", dest="resname", type="string", default="", help="select atoms of group by residue name  ")    
    parser.add_option("--residue", dest="residue", type="string", default="", help="select atoms of group by resudue number  ")    
    parser.add_option("--linkid", dest="linkid", type="string", default="", help="select atoms of group by  linkd ")    
    parser.add_option("--fftype", dest="fftype", type="string", default="", help="select atoms of group by force field type  ")
    parser.add_option("--gtype", dest="gtype", type="string", default="", help="select atoms of group by force field type  ")
    #
    # Frame selection 
    #
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=False,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")

    (options, args) = parser.parse_args()
        
    return options, args
   
def aveprop():
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
    if( options.verbose ):
        print " Reading in files"
    
    struc_o = StructureContainer()
    param_o = ParameterContainer()
    
    struc_o,param_o = file_io.getstrucC(struc_o,param_o, options.in_json, options.in_gro , options.in_top, options.in_itp, options.in_data,options.in_xmol,options.xmol_format )
    #struc_o = file_io.getstrucC(struc_o, options.in_json, options.in_gro , options.in_top, options.in_data,options.in_xmol,options.xmol_format )
    #
    # Get paticle and bond structures
    #
    ptclC_o = struc_o.ptclC
    bondC_o  = struc_o.bondC
    
    # 
    # Read in trajectory 
    #
    if( len(options.in_xtc ) and  len(options.in_gro) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s and %s "%(options.in_xtc,options.in_gro)
        universe =  Universe(options.in_gro, options.in_xtc)
    elif(len(options.in_gro) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s "%(options.in_gro)
        universe =  Universe(options.in_gro)
    elif( len(options.in_dcd ) and  len(options.in_data) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s and %s "%(options.in_xtc,options.in_gro)
        universe =  Universe(options.in_dcd, options.in_data)
    elif(len(options.in_dcd) ):
        if( rank == 0 ):
            if( options.verbose ): print "  Reading  %s "%(options.in_dcd)
        universe =  Universe(options.in_dcd)

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
        # sys_prop = struc_o.printprop()
        struc_o.verbose=False
        sys_prop = str(struc_o)
        print sys_prop
        log_out.write(str(sys_prop))
        struc_o.verbose=True
        
    #   
    # Filter particles
    #
    if( options.verbose ):
        print " Filter particles"

    search_o = particles.create_search(options.id,options.type,options.symbol,options.chains,options.ring,options.resname,options.residue,options.linkid,options.fftype,options.gtype)
    if( rank == 0 ):
        if( options.verbose ): print " Filter input by ",search_o,len(search_o)
    list_f = ptclC_o.getParticlesWithTags(search_o)
    sum_f = len(list_f)

    

    for ts in universe.trajectory:
        if( options.frame_o <= ts.frame ):
            if( ts.frame <= options.frame_f  or  options.readall_f  ):
                if( ts.frame%options.frame_step == 0 ):

                    print "Frame %4d with volume %f " % (ts.frame, ts.volume)
                    rdf_frames += 1 
                    volume_sum_i += ts.volume      # correct unitcell volume
                    box = ts.dimensions

    
    if( rank == 0 and options.verbose  ):
        print "  searching finished "
        print "  Creating list of atomic indies of dihedrals "
    

    p.barrier()
    
    if( rank == 0 ):

        # log_line += "\n  %d "%()
        log_line  = "\n  Date %s "%(str(t_i))
        log_line += "\n  Number of processors  %d "%(size)

        log_line += sperator_line
        log_line += "\n  Particles of group : %d "%(sum_f)

        log_line += sperator_line
        log_line += "\n Frames "
        log_line += "\n     Initial frame  %d "%(options.frame_o)
        log_line += "\n     Step frame  %d  "%(options.frame_step)
        log_line += "\n     Final frame  %d  "%(options.frame_f)

        log_out.write(log_line)

            
        if( options.verbose ):
            print log_line
        log_out.close()
        dat_out.close()
    
	
	    
if __name__=="__main__":
    aveprop()
   
	    
