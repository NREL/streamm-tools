#! /usr/bin/env python
"""
Translate between various file types
"""

# Dr. Travis Kemper
# Initial Date April 2014
# travis.kemper@nrel.gov

import sys 
import datetime
from string import replace

from structureContainer import StructureContainer
from parameters import ParameterContainer

import particles 
import mpiBase, lammps , gromacs , file_io , xmol , topology 

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
    #
    # Input files 
    #
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_itp", dest="in_itp",  type="string", default="",help="Input gromacs force field parameter file")
    parser.add_option("--in_fchk", dest="in_fchk", type="string", default="", help="Input gromacs .fchk file ")
    # Need to be added with mdanalysis dependence 
    #parser.add_option("--in_dcd", dest="in_dcd", type="string", default="", help="Input trajectory file in compressed dcd format ")
    #parser.add_option("--in_xtc", dest="in_xtc", type="string", default="", help="Input trajectory file in compressed xtc format ")
    parser.add_option("--in_xmol", dest="in_xmol", type="string", default="", help="Input xmol file")
    parser.add_option("--xmol_format", dest="xmol_format", type="string", default="atomic_symb", help="Format of xmol file: atomic_symb - first colum atomic symbols; lammps - first colum particle type atomic symbols will be set to carbon; ")
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")
    #
    # Output files 
    #
    parser.add_option("--out_data", dest="out_data",type="string",default="",help=" Output Lammps data file ")
    parser.add_option("--out_gro", dest="out_gro", type="string", default="", help=" Output gromacs structure file ")
    parser.add_option("--out_top", dest="out_top", type="string", default="", help=" Output gromacs topology file ")
    parser.add_option("--out_itp", dest="out_itp", type="string", default="", help=" Output gromacs parameter file ")
    parser.add_option("--out_json", dest="out_json", type="string", default="", help=" Output json file ")
    #
    # xmol output 
    #
    parser.add_option("--out_xmol", dest="out_xmol", type="string", default="", help=" Output xmol file of trajectory ")
    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help=" Output xyz file of the final structure ")
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=False,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")
    #
    #
    #
    #
    # Filters
    #
    parser.add_option("--id", dest="id", type="string", default="", help=" select atoms of group by number  ")    
    parser.add_option("--symbol", dest="symbol", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--type", dest="type", type="string", default="", help=" select type  ")    
    parser.add_option("--chains", dest="chains", type="string", default="", help="select atoms of group by chain/molecule number  ")    
    parser.add_option("--ring", dest="ring", type="string", default="", help="select atoms of group by particlesn a ring   ")    
    parser.add_option("--resname", dest="resname", type="string", default="", help="select atoms of group by residue name  ")    
    parser.add_option("--residue", dest="residue", type="string", default="", help="select atoms of group by resudue number  ")    
    parser.add_option("--linkid", dest="linkid", type="string", default="", help="select atoms of group by  linkd ")    
    parser.add_option("--fftype", dest="fftype", type="string", default="", help="select atoms of group by force field type  ")
    parser.add_option("--gtype", dest="gtype", type="string", default="", help="select atoms of group by force field type  ")    
    parser.add_option("--filter_lmptype", dest="filter_lmptype", type="string", default="", help=" filter atoms by lammps type ")
 
    #
    # Modify input system file 
    #
    
    (options, args) = parser.parse_args()
        
    return options, args
   

def main():
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
    p = mpiBase.getMPIObject()

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
    param_o = ParameterContainer() 
    
    struc_o,param_o = file_io.getstrucC(struc_o,param_o, options.in_json, options.in_gro , options.in_top, options.in_itp, options.in_data, options.in_fchk ,options.in_xmol,options.xmol_format )

    p.barrier()
    #
    # Get paticle and bond structures
    #
    ptclC_o = struc_o.ptclC
    bondC_o  = struc_o.bondC

    #   
    # Filter particles
    #
    search_o = particles.create_search(options.id,options.type,options.symbol,options.chains,options.ring,options.resname,options.residue,options.linkid,options.fftype,options.gtype)
    if( rank == 0 ):
        if( options.verbose ): print " Filter input by ",search_o

    if( len(search_o) > 0 ):
        list_f = ptclC_o.getParticlesWithTags(search_o)
        sum_f = len(list_f)
        struc_i = struc_o.getSubStructure(list_f)
    else:
        struc_i = struc_o
        sum_f =  len(struc_i.ptclC)
    ptclC_i = struc_i.ptclC
    # SWS: example of new call struc_i = struc_o.getSubStructure(list_f, particlesOnly=True)
    #
    # Open output files 
    #
    if( rank == 0 ):

        log_file = options.output_id + ".log"
        log_out = open(log_file,"w")
        
        # xmol file 
        if( len(options.out_xmol) ):
            if( options.verbose ): print  "     - Writing  ",options.out_xmol
            append_xmol = True
            # Delete xmol file if exsists
            str_file = open( options.out_xmol, 'w' )
            str_file.close()
    #
    # Print settings 
    #
    if( rank == 0 ):
         # Print initial structure properties
        sys_prop = str( struc_o )
        log_out.write(str(sys_prop))

        log_line  = "\n  Date %s "%(str(t_i))
        log_line += "\n  Number of processors  %d "%(size)
        log_line += "\n  Particles in filtered structure : %d "%(sum_f)
        log_line += sperator_line
        log_line += "\n Frames "
        log_line += "\n     Initial frame  %d "%(options.frame_o)
        log_line += "\n     Step frame  %d  "%(options.frame_step)
        log_line += "\n     Final frame  %d  "%(options.frame_f)

        log_out.write(log_line)
            
        if( options.verbose ):
            print sys_prop
            print log_line
            
    #
    #  Read in trajectories and process output
    #
    if( len(options.in_lammpsxyz) ):
        
        if( rank == 0 ):
            if( options.verbose ): print "  Reading ",options.in_lammpsxyz
            lmpxyz_t_i = datetime.datetime.now()
            
        # Initialize line count
        NP = len( ptclC_o )
        line_cnt = 0
        frame_cnt = 1
        pt_cnt = 0
        read_frames = 0 
        # Read in file line by line 
        with open(options.in_lammpsxyz) as f:
            for line in f:
                line_cnt += 1
                col = line.split()

                if( line_cnt == 1  ):
                    # Read first line to record number of particles NP
                    NP = int(col[0])
                    # frame_cnt += 1
                    pt_cnt = 0
                    if( rank == 0 and options.verbose ):
                        # Record start time of frame read in 
                        fread_t_i = datetime.datetime.now()                    
                elif( line_cnt > 2 and len(col) >= 4 ):
                    # Read lines and update particle postions 
                    pt_cnt += 1
                    pt_o = ptclC_o[pt_cnt]
                    pt_o = struc_o.ptclC[pt_cnt]
                    pt_o.position =  [ col[1],col[2],col[3] ]
                    #struc_o.ptclC[pt_cnt].position =  [ col[1],col[2],col[3] ]
                    atomic_symb = col[0]
                    #
                    #r_i = [ col[1],col[2],col[3] ]
                    #part_i = Particle( r_i,atomic_symb,q_i,m_i )
                    #struc_i.put(part_i)
                    #print " reading lmpxyz ",pt_cnt,pt_o.postion, struc_o.ptclC[pt_cnt].position

                if( line_cnt > 1  and line_cnt > frame_cnt*(NP + 2) ):
                    # End of frame found 
                    pt_cnt = 0
                    frame_cnt += 1
                    
                    if( options.frame_o <= frame_cnt ):
                        if( frame_cnt <= options.frame_f  or  options.readall_f  < 0 ):
                            if( frame_cnt%options.frame_step == 0 ):


                                print " frame ",frame_cnt," read frame ",read_frames

                                
                                # # Delete structure with frame information 
                                del struc_i
                                struc_i = struc_o.getSubStructure(list_f)

                                print " atom 1  struc_i",struc_i.ptclC[1]
                                # Update postion of struc_i
                                #pt_cnt = 0
                                #for p_o, ptcl_o in ptclC_o(list_f):
                                #    pt_cnt += 1
                                #    pt_i = struc_i.ptclC[pt_cnt]
                                #    pt_i.postion =   ptcl_o.postion 
                                
                                if( len(options.out_xmol) ):
                                    # Print xmol 
                                    comment = " Frame %d "%(frame_cnt)
                                    struc_i.write_xmol(options.out_xmol,comment,append_xmol)

                                # Print output
                                if( rank == 0 and options.verbose ):
                                    fread_t_f = datetime.datetime.now()
                                    dt_sec  = fread_t_f.second - fread_t_i.second
                                    dt_min  = fread_t_f.minute - fread_t_i.minute
                                    if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                    if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                    log_line = "\n       Frame  %d  read in %d min  %f sec "%(frame_cnt,dt_min,dt_sec)
                                    log_out.write(log_line)
                                    print log_line
                                    # Rest intial time for next rdf calcution
                                    fread_t_i = datetime.datetime.now()

                                # Index read frames 
                                read_frames += 1
                                    

                
            # Process last frame             
            if( options.frame_o <= frame_cnt ):
                if( frame_cnt <= options.frame_f  or  options.readall_f  < 0 ):
                    
                    if( frame_cnt%options.frame_step == 0 ):
                        read_frames += 1
                        
                        # Delete structure with frame information 
                        del struc_i
                        struc_i = struc_o.getSubStructure(list_f)

                                
                        if( len(options.out_xmol) ):
                            comment = " Frame %d "%(frame_cnt)
                            struc_i.write_xmol(options.out_xmol,comment,append_xmol)

                        # Print output
                        if( rank == 0 and options.verbose ):

                            fread_t_f = datetime.datetime.now()
                            dt_sec  = fread_t_f.second - fread_t_i.second
                            dt_min  = fread_t_f.minute - fread_t_i.minute
                            if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                            if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                            log_line = "\n       Frame  %d  read in %d min  %f sec "%(frame_cnt,dt_min,dt_sec)
                            log_out.write(log_line)
                            print log_line
                            
    #
    #  Find new parameter set for new structure
    #
    if( len(options.out_top) or  len(options.out_itp)  or  len(options.out_data) ):
        ff_charges = False 
        struc_i,cov_nblist, cov_nbindx = topology.create_top(struc_i,ff_charges)
        #.build_covnablist(struc_i)
        norm_dihparam = False
        param_i,struc_i  = topology.set_param(struc_i,param_o,norm_dihparam,cov_nblist, cov_nbindx)
    else:
        param_i = param_o
        
    debug = False
    if( debug ):

        for pid_o, ptclObj_o  in struc_i.ptclC:  
            print  ptclObj_o.tagsDict["fftype"],ptclObj_o.tagsDict["lmptype"]

        sys.exit(" lmptype debug 1 ")

    debug = False
    if( debug ):

        for pid_o, ptclObj_o  in param_i.ljtypC:  
            print "  lj object add type mass ",ptclObj_o.get_ptype1(),ptclObj_o.get_mass()
        sys.exit(" lmptype mass 1 ")


    if( rank == 0 ):
        #
        # Output files 
        #
        prop_struc_i = str(struc_i)
        prop_param_i = str(param_i)
        print prop_struc_i
        print prop_param_i
        #  Write gro file         
        if( len(options.out_gro) ):
            if( options.verbose ): print  "     - Writing  ",options.out_gro
            #for pid, ptclObj  in struc_i.ptclC:
            #     ptclObj.tagsDict["gtype"] =  ptclObj.tagsDict["fftype"]
            gromacs.print_gro(struc_i,options.out_gro)
            #struc_i.write_gro(dir_id,options.out_gro ) # Move out of class
            #gromacs.print_gro(options.out_gro,GTYPE_sys,RESID_sys,RESN_sys,R_sys,LV)

        #  Write top file         
        if( len(options.out_top) ):
            if( len(options.out_itp) > 0 ):
                itp_file = options.out_itp
            else:
                itp_file = "ff-new.itp"
                
            if( options.verbose ): print  "     - Writing %s with itp file %s  "%(options.out_top,itp_file)
            gromacs.print_top(struc_i,options.out_top,itp_file)

        #  Write itp file         
        if( len(options.out_itp) ):
            if( options.verbose ): print  "     - Writing  ",options.out_itp
            gromacs.print_itp(param_i,options.out_itp)

        #  Write xmol file 
        if( len(options.out_xyz) ):
            if( options.verbose ): print  "     - Writing  ",options.out_xyz
            comment = " final structure "
            append = False 
            xmol.write(struc_i.ptclC,options.out_xyz,comment,append)
            xmol.write(struc_i.ptclC,options.out_xyz,comment,append)            

        #  Write Lammps data file
        if( len(options.out_data) ):
            if( options.verbose ): print  "     - Writing  ",options.out_data        
            lammps.write_data(struc_i,param_i,options.out_data)

        # Write json file
        if( len(options.out_json) ):
            if( options.verbose ):
                print  "     - Writing  ",options.out_json
            dir_id =  "./"
            output_id = replace(options.out_json,".json","")
            struc_i.write_json(dir_id,output_id )
        
	if( options.verbose ):
            print " "
        t_f = datetime.datetime.now()
        dt_sec  = t_f.second - t_i.second
        dt_min  = t_f.minute - t_i.minute
        if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
        log_line="\n  Finished time  " + str(t_f)
        log_out.write(log_line)
        log_line="\n  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "
        log_out.write(log_line)
        print log_line
	log_out.close()


if __name__=="__main__":
    main()
   
	    
