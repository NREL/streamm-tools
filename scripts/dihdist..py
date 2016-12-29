#! /usr/bin/env python
"""
Analysis of time series properties 
"""



__author__ = "Travis W. Kemper"
__version__ = "0.2"
__email__ = "travkemp@gmail.com"
__status__ = "Alpha"


import   os, os.path, sys , copy ,shutil, logging
import json, math , sys
import numpy as np
from datetime import datetime
from optparse import OptionParser

import project, buildingblock, simulation, structure
import mpiBase

import stmm_pandas
import stmm_mdanalysis


def timeseries(tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
        
    if( rank == 0 ):
        logging.info('Running on %d procs %s '%(size,datetime.now()))
        logging.info('Reading structure files ')
        
    logging.debug(" proc %d of %d "%(rank,size))

    
    strucC = buildingblock.Container(tag)
    strucC.read_cply(options.cply)

    # Convert to pandas and select sub lists 
    stmm_pandas.convert_strucC_df(strucC)
    s1 = strucC.pd_df['fftype'] == 'C!'
    # Sub selection i 
    sub_i =  strucC.pd_df[s1]
    list_i = sub_i.index
    N_i = len(list_i)
    list_j = sub_i.index
    N_j = len(list_j)

    # CS CB
    s2 = strucC.pd_df['fftype'] == 'CS'
    sub_k =  strucC.pd_df[s2]
    list_k = sub_k.index
    s3 = strucC.pd_df['fftype'] == 'CB'
    sub_l =  strucC.pd_df[s3]
    list_l = sub_l.index
    
    # Create list of dihderal angles 
    strucC.get_dihatoms(list_k,list_i,list_j,list_l)
    strucC.write_dihatoms()
        
    # Select all particles 
    strucC.particles_select = strucC.particles.keys()
    # Make groups by chain for pbcs 
    group_index = strucC.group_prop("chain")
    strucC.group_pbcs()
    strucC.write_xyz()
    prop_dim = strucC.lat.n_dim
    #
    # Write gromacs file .gro for MDanalysis 
    #
    gro_file = "%s.gro"%(tag)
    if( rank == 0 ):
        # Write gromacs gro file for MDANALYSIS 
        gromacs_i  = simulation.GROMACS(tag)
        gromacs_i.add_strucC(strucC)
        gromacs_i.write_gro(gro_file=gro_file)
    p.barrier()


    '''
    # Get posiition and lattice constants from mdanalysis 
    pos_frames,box_frames = stmm_mdanalysis.get_frames(gro_file, options.dcd,options.frame_o,options.frame_step,options.frame_f,options.readall_f,rank)


    n_frames = len(pos_frames)
    for frame_i in range(n_frames):
        logger.info(" Analyizing frame %d/%d "%(frame_i,n_frames))
        strucC.positions = pos_frames[frame_i]
        box_i = box_frames[frame_i]
        strucC.lat.LCtoLV(box_i)
    '''

    strucC.calc_cosdihs()
        
        
    # Normalize
    strucC.ave_cosdihs = np.zeros(strucC.n_dih)
    strucC.std_cosdihs = np.zeros(strucC.n_dih)
    for d_i in range(strucC.n_dih):
        strucC.ave_cosdihs[d_i] = np.average(strucC.cosdihs[d_i])
        strucC.std_cosdihs[d_i] = np.std(strucC.cosdihs[d_i])
        
                
    # Print results
    dat_line = "# n_dih %d  \n"%(strucC.n_dih)
    dat_line += "# d_i , <cosine_kijl> , sigma(cosine_kijl) \n"
    for d_i in range(strucC.n_dih):
        cosine_kijl  = strucC.ave_cosdihs[d_i]
        std_kijl  = strucC.std_cosdihs[d_i]
        ang_rad = np.arccos(cosine_kijl )
        ang_deg = np.rad2deg( ang_rad )
        
        dat_line += " %d  %f  %f %f \n"%(d_i,cosine_kijl,std_kijl,ang_deg)
    dat_file = open(str("dih_%s.dat"%(tag)),"w")
    dat_file.write(dat_line)
    dat_file.close()

if __name__=="__main__":
    
    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--cply", dest="cply", type="string", default="min1.cply", help="Input cply file")
    parser.add_option("--dcd", dest="dcd", type="string", default="nvt1_dump.dcd", help="Input dcd file")
    # Frames
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=10000, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    parser.add_option("--readall_f", dest="readall_f", default=True,action="store_true", help=" Read to end of trajectory file (negates frame_f value)")
    # Bins
    parser.add_option("--r_cut", dest="r_cut", type=float, default=20.0, help=" Cut off radius in angstroms ")
    parser.add_option("--bin_size", dest="bin_size", type=float, default=0.10, help=" Bin size in angstroms")
    # 
    (options, args) = parser.parse_args()
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()


    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    #logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # 
    # logging.basicConfig(level=logging.DEBUG,
    #                format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    #                datefmt='%m-%d %H:%M',
    #                filemode='w')

    if( len(args) < 1 ):
        calc_tag = 'rgy'
    else:
        calc_tag =  args[0]
        
    calc_tag = 't2'

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)


    hdlr = logging.FileHandler('%s.log'%(calc_tag),mode='w')
    hdlr.setLevel(logging.DEBUG)
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)

    start_time = datetime.now()
    logger.info('Started %s '%(start_time))
    

    #rdf(calc_tag,options)
    timeseries(calc_tag,options,p)
    
    finish_time = datetime.now()
    delt_t = finish_time - start_time
    
    logger.info('Finished %s in %f seconds '%(finish_time,delt_t.seconds))