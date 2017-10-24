#! /usr/bin/env python
"""

Radial distribution  code

 length - Angstroms
 mass   - AMU
 volume - Angstroms^3

    # g_ij(r) = n_j(r_ij)/ (rho_j 4 pi r^2 dr )
    g_ij(r) = n_j(r_ij)/ (rho_j 4/3 pi( r_out^3 - r_in^3)

    rho_j = N_j / V_ave

    g(r) = n(r) /( rho dV )

    n(r) = 1/Ni sum_i^{N_i} sum_j^{N_j} \gamma( r - r_{ij})

    g(r) = 1/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) / rho_j(i) ]

    Regular density
    g(r) =  1/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) / (  sum_j^{N_j} \gamma( allowed pair ij )/<V> )  ]
    g(r) =  <V>/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) /   sum_j^{N_j} \gamma(  pair ij ) ]
    
    True density 
    g(r) =  1/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) / (  sum_j^{N_j} \gamma( allowed pair ij )/<V> )  ]
    g(r) =  <V>/( dV Ni )  sum_i^{N_i}    [   sum_j^{N_j} \gamma( r - r_{ij}) /   sum_j^{N_j} \gamma( allowed pair ij ) ]


    Nj_i = sum_j^{N_j} \gamma(  pair ij )
    rdf_cnt_p = sum_f^{N_frames}   sum_j^{N_j} \gamma( r - r_{ij})
    sum_j^{N_j} \gamma( r - r_{ij})  = rdf_cnt_p/N_frames
 
"""


__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"


import   os, os.path, sys , copy ,shutil, logging, math, json, csv 
import numpy as np
from datetime import datetime
from optparse import OptionParser

from streamm import *

from MDAnalysis import *
from MDAnalysis.core.distances import * ##distance_array

import pandas as pd 

def grdfs(tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    #
    if( rank == 0 ):
        logging.info('Running on %d procs %s '%(size,datetime.now()))
        logging.info('Reading structure files ')
    #
    logging.debug(" proc %d of %d "%(rank,size))
    if( rank == 0 ):
        logger.info("Reading in structure for %s from %s "%(tag,options.cply))
    #
    # Read cply file 
    #
    strucC = buildingblock.Container()
    strucC.read_cply(options.cply)
    # 
    # Calculate bulk properties 
    # 
    strucC.bonded_nblist.build_nblist(strucC.particles,strucC.bonds )
    strucC.calc_mass()
    strucC.calc_volume()
    strucC.calc_density()
    prop_dim = strucC.lat.n_dim
    den_gcm3 = units.convert_AMUA3_gcm3(strucC.density)
    if( rank == 0 ):
        logger.info("Structure %s mass %f volume %f density %f "%(tag,strucC.mass,strucC.volume,den_gcm3))
        logger.info("Breaking molecules into separate simulations")
    if( rank == 0 ):
        logger.info("Selecting groups %s "%(options.group_id))
    # 
    '''
    group_i_id = 'mol'
    strucC.group_prop(group_i_id,group_i_id)
    groupset_mol = strucC.groupsets[group_i_id]
    groupset_mol.calc_cent_mass()
    groupset_mol.calc_radius()
    groupset_mol.group_pbcs()
    '''
    #
    # Write .gro file for MDanalysis 
    # 
    gro_file = "%s.gro"%(tag)
    if( rank == 0 ):
        # Write gromacs gro file for MDANALYSIS 
        gromacs_i  = gromacs.GROMACS(tag)
        gromacs_i.add_strucC(strucC)
        gromacs_i.write_gro(gro_file=gro_file)
    p.barrier()
    #
    # Read lists
    #
    if( len(options.p_list) > 0 ):
        lfile = open(options.p_list,'rb')
        t_i = lfile.readlines()
        lfile.close()
        p_list = [int(pkey) for pkey in  t_i]
    else:
        p_list = strucC.particles.keys()
    #
    # Select considered particles 
    #
    if( rank == 0 ):
        logger.info("Grouping by %s "%(options.group_id))        
    strucC.group_prop(options.group_id,options.group_id,particles_select=p_list)
    groupset_i = strucC.groupsets[options.group_id]
    groupset_i.calc_cent_mass()
    groupset_i.calc_radius()        
    #
    # Read lists
    #
    if( len(options.glist_i) > 0 ):
        lfile = open(options.glist_i,'rb')
        t_i = lfile.readlines()
        lfile.close()
        glist_i = [int(pkey) for pkey in  t_i]
    else:
        glist_i = groupset_i.groups.keys()
    # 
    if( len(options.glist_j) > 0 ):
        lfile = open(options.glist_j,'rb')
        t_j = lfile.readlines()
        lfile.close()
        glist_j = [int(pkey) for pkey in  t_j]
    else:
        glist_j = groupset_i.groups.keys()
        
    #sub_j =  strucC.pd_df[sub_heavy & D1 ]
    # Get glist and pairs
    pairvalue_ij = groupset_i.find_pairs(glist_i,glist_j,mol_inter=options.mol_inter,mol_intra=options.mol_intra)
    #
    # Find group neighbor list 
    #
    groupset_i.group_nblist.radii_nblist(strucC.lat,groupset_i.properties['cent_mass'],groupset_i.properties['radius'],radii_buffer=options.pairbuffer)
    gkeys_p = p.splitListOnProcs( groupset_i.groups.keys())
    #
    # Write pairs 
    #
    pairs_file = 'pairs_%s.csv'%(options.group_id)
    if( rank == 0 ):
        fout = open(pairs_file,'wb')
        pair_writer = csv.writer(fout,delimiter=',')
        header = ['g_i','g_j','dcm_ij']
        #if( rank == 0 ):
        pair_writer.writerow(header)
        fout.close()
        logger.info('file: output pairs_%s %s '%(options.group_id,pairs_file))


    logger.debug(" Writing %d group pairs on proc %d "%(len(gkeys_p),rank))
    for g_i in gkeys_p:
        group_i = groupset_i.groups[g_i]
        mol_i =group_i.properties['mol']
        # 
        fout = open(pairs_file,'a')
        pair_writer = csv.writer(fout,delimiter=',')
        #
        nb_cnt = groupset_i.group_nblist.calc_nnab(g_i)
        logger.debug(" group %d has %d nieghbors with radius of %f "%(g_i,nb_cnt,group_i.properties['radius']))
        for g_j in groupset_i.group_nblist.getnbs(g_i):
            logger.debug("checking neighbor group %d "%(g_j))
            if( g_j > g_i ):
                group_j  = groupset_i.groups[g_j]
                mol_j =group_j.properties['mol']
                if( mol_i != mol_j):
                    dcm_ij = groupset_i.group_nblist.dist_matrix[g_i,g_j] 
                    row_i = [g_i,g_j,dcm_ij]
                    pair_writer.writerow(row_i)
  
        fout.close()    #
        
    pairs_df = pd.read_csv('pairs_residue.csv')
        
    #Calculate rdfs  
    #
    # gr2_rdf = rdf(gr2_tag,gr2_i,gr2_i_p,gr2_j,gr2_pairs,gro_file,options,p)
    #
    logger.info(" Calculating close_contacts %s "%(tag))
    #
    # for tag_i,calc_j in proj_j.calculations.iteritems():
    # os.chdir(calc_j.dir['scratch'])
    gdr_ij = dict()
    gdr_ij['g_i'] = []
    gdr_ij['g_j'] = []
    gdr_ij['dcm_ij'] = []
    gdr_ij['dr_pi_pj'] = []
    #
    #
    for (g_i,g_j,dcm_ij) in pairs_df[['g_i','g_j','dcm_ij']].values:
        dr_pi_pj = groupset_i.dr_particles(g_i,g_j,options.r_cut,p_list)
        gdr_ij['g_i'].append(g_i)
        gdr_ij['g_j'].append(g_j)
        gdr_ij['dcm_ij'].append(dcm_ij)
        gdr_ij['dr_pi_pj'].append(dr_pi_pj)
    
    df_gdr_ij = pd.DataFrame(gdr_ij)
    file_name = "gdr_ij.csv"
    df_gdr_ij.to_csv(file_name, sep=',', encoding='utf-8')    
        
    
    #bin_r,bin_r_nn,volumes = distbin_pairs(glist_i,glist_j,pairs_ij,gro_file, options.dcd,options.frame_o,options.frame_step,options.frame_f,options.readall_f,options.bin_size,options.r_cut,rank)
    #bin_r,bin_r_nn,volumes = singleframe_gpairs(strucC,options.group_id,glist_i,glist_j,pairvalue_ij,options.bin_size,options.r_cut,rank)
    
    p.barrier()
    
            
if __name__=="__main__":
    
    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--cply", dest="cply", type="string", default="in.cply", help="Input cply file")
    parser.add_option("--p_list", dest="p_list", type="string", default="p_list", help="Input list file")
    parser.add_option("--dcd", dest="dcd", type="string", default="prod1_dump.dcd", help="Input dcd file")
    parser.add_option("--mol_inter",dest="mol_inter", default=False,action="store_true", help="Use only inter molecular rdf's ")
    parser.add_option("--mol_intra",dest="mol_intra", default=False,action="store_true", help="Use only intra molecular rdf's")
    parser.add_option("--groups_inter",dest="groups_inter", default=True,action="store_true", help="Use only inter group rdf's")
    parser.add_option("--group_id", dest="group_id", type="string", default="mol", help="Group id ")
    parser.add_option("--pairbuffer", dest="pairbuffer", type="float", default=2.5, help="Pair buffer ")
    # parser.add_option("--truedensity",dest="truedensity", default=False,action="store_true", help="Use the true density of group j, this makes inter/intra molecular rdfs not components of the total rdf but true independent rdfs")
    # List of groups 
    parser.add_option("--glist_i", dest="glist_i", type="string", default="", help="Input list file")
    parser.add_option("--glist_j", dest="glist_j", type="string", default="", help="Input list file")
    
    
    # Frames
    parser.add_option("--frame_o", dest="frame_o", type=int, default=1, help=" Initial frame to read")
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
        calc_tag = 'rdf'
    else:
        calc_tag =  args[0]
        
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    #
    hdlr = logging.FileHandler('%s.log'%(calc_tag),mode='w')
    hdlr.setLevel(logging.INFO)
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    #
    start_time = datetime.now()
    if( p.getRank() == 0 ):
        logger.info('Started %s '%(start_time))
    #
    grdfs(calc_tag,options,p)
    
    finish_time = datetime.now()
    delt_t = finish_time - start_time
    
    if( p.getRank() == 0 ):
        logger.info('Finished %s in %f seconds '%(finish_time,delt_t.seconds))
