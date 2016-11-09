#! /usr/bin/env python
"""
Analysis of time series properties 
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

def file_test(file_i):
    '''
    Read file into string 
    '''
    try:
        with open(file_i) as F:            
            F.close()
            return False  

    except IOError:
        print " File not found %s "%(file_i)

    return True

def string_load(template_file):
    '''
    Read file into string 
    '''
    file_string = ''
    try:
        with open(template_file) as F:            
            template_lines = F.readlines()
            for line in template_lines:
                file_string += line
            F.close()

    except IOError:
        logger.info(" File not found %s "%(template_file))
        print " File not found %s "%(template_file)
        

    return file_string

def read_group_cplys(group_id):
    
    group_tags = []
    group_file = "group_%s.csv"%(group_id)
    logger.debug("Reading group cplys   file %s "%(group_file))
    with open(group_file, 'rb') as f:
                reader = csv.reader(f)
                rownum = 0
                for row in reader:
                    if rownum == 0:
                        header = row
                    else:
                        group_tags.append(row)
                    rownum += 1
                    
    group_list = []
    for ginfo_i in group_tags:
        gtag_i = ginfo_i[1]
        # print " Reading in ",gtag_i
        bb_i = buildingblock.Container(gtag_i)
        bb_i.read_cply()
        group_list.append(copy.deepcopy(bb_i))
    
    return group_tags,group_list


def read_pairs(group_id):
    
    pair_tags = []
    pairs_file = 'pairs_%s.csv'%(options.group_id)
    logger.debug("Reading pairs from file %s "%(pairs_file))
    with open(pairs_file, 'rb') as f:
                reader = csv.reader(f)
                rownum = 0
                for row in reader:
                    if rownum == 0:
                        header = row
                    else:
                        pair_tags.append([int(row[0]),int(row[1])])
                    rownum += 1
                
    return pair_tags

def read_sims(sims_file):
    
    sim_tags = []
    logger.debug("Reading pairs from file %s "%(sims_file))
    with open(sims_file, 'rb') as f:
                reader = csv.reader(f)
                rownum = 0
                for row in reader:
                    if rownum == 0:
                        header = row
                    else:
                        sim_tags.append(str(row[0]))
                        # print ">read_sims sim_tags",str(row[0])
                    rownum += 1
                
    return sim_tags
        
def pull_groups(tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    #
    # Read in data 
    #
    if( rank == 0 ):
        logging.info('Running pull_groups on %d procs %s '%(size,datetime.now()))
        logging.info('Reading structure files ')
    # 
    logging.debug(" proc %d of %d "%(rank,size))
    if( rank == 0 ):
        logger.info("Reading in structure for %s from %s "%(tag,options.cply))
    struc_o = buildingblock.Container('%s'%(tag))
    struc_o.read_cply(options.cply)
    # 
    if( rank == 0 ):
        logger.info("Building neighbor list from bonds in data file %d  "%(struc_o.n_bonds))
    #
    # Calculate bulk properties 
    #
    struc_o.bonded_nblist.build_nblist(struc_o.particles,struc_o.bonds )
    struc_o.calc_mass()
    struc_o.calc_volume()
    struc_o.calc_density()
    prop_dim = struc_o.lat.n_dim
    den_gcm3 = units.convert_AMUA3_gcm3(struc_o.density)
    if( rank == 0 ):
        logger.info("Structure %s mass %f volume %f density %f "%(tag,struc_o.mass,struc_o.volume,den_gcm3))
        logger.info("Breaking molecules into separate simulations")
    if( rank == 0 ):
        logger.info("Selecting groups %s "%(options.group_id))
    #
    # Select considered particles 
    #
    if( rank == 0 ):
        logger.info("Writing %s group cplys  "%(options.group_id))
    struc_o.propcompile_particles()
    # Read lists
    lfile = open(options.list_i,'rb')
    t_i = lfile.readlines()
    lfile.close()
    list_i = [int(pkey) for pkey in  t_i]
    #
    # Find groups 
    #
    struc_o.group_prop(options.group_id,options.group_id,particles_select=list_i)
    groupset_i = struc_o.groupsets[options.group_id]
    #
    # Calculate group properties 
    #
    groupset_i.calc_cent_mass()
    groupset_i.calc_radius()
    groupset_i.group_pbcs()
    
    # Set radii to fixed value
    et_cut = 10.0 
    groupset_i.properties['radius'] = []
    for gkey,group_i in groupset_i.groups.iteritems():
        groupset_i.properties['radius'].append(et_cut)
    #
    # Write group properties 
    #
    if( rank == 0 ):
        #groupset_i.write_cm_xyz()
        #groupset_i.write_xyzs()
        groupset_i.dump_json()
    #
    group_file = "group_%s.csv"%(options.group_id)
    if( rank == 0 ):
        fout = open(group_file,'wb')
        #group_writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        group_writer = csv.writer(fout,delimiter=',')
        header = ['g_i','tag','mol','x_cent_mas','y_cent_mas','z_cent_mas']
        #if( rank == 0 ):
        group_writer.writerow(header)
        fout.close()
        logger.info('file: output group_%s %s '%(options.group_id,group_file))

    
    gkeys_p = p.splitListOnProcs( groupset_i.groups.keys())
    
    #for g_i,group_i in groupset_i.groups.iteritems():
    for g_i in gkeys_p:
        group_i = groupset_i.groups[g_i]
        fout = open(group_file,'a')
        group_writer = csv.writer(fout,delimiter=',')
        
        tag_i = "%s_%i"%(options.group_id,g_i)
        if( rank == 0 ):
            logger.info("Writing cply of group %s "%(tag_i))
        if( options.hterm ):
            bb_i = group_i.hterm_group()
        else:
            bb_i  = struc_o.getSubStructure(group_i.pkeys,tag=tag_i)
        bb_i.tag = tag_i
        # bb_i.guess_oplsa()
        bb_i.write_cply()
        bb_i.write_xyz()
        mol_i = bb_i.particles[0].properties['mol']
        row_i = [g_i,tag_i,mol_i]
        row_i += group_i.properties['cent_mass']
        group_writer.writerow(row_i)
        fout.close()
    
    #
    # Find group neighbor list 
    #
    groupset_i.group_nblist.radii_nblist(struc_o.lat,groupset_i.properties['cent_mass'],groupset_i.properties['radius'],radii_buffer=1.0)
    #
    # Write pairs 
    #
    pairs_file = 'pairs_%s.csv'%(options.group_id)
    if( rank == 0 ):
        fout = open(pairs_file,'wb')
        pair_writer = csv.writer(fout,delimiter=',')
        header = ['g_i','g_j','mol_i','mol_j']
        #if( rank == 0 ):
        pair_writer.writerow(header)
        fout.close()
        logger.info('file: output pairs_%s %s '%(options.group_id,pairs_file))


    for g_i in gkeys_p:
        group_i = groupset_i.groups[g_i]
        # 
        fout = open(pairs_file,'a')
        pair_writer = csv.writer(fout,delimiter=',')
        
        nb_cnt = groupset_i.group_nblist.calc_nnab(g_i)
        logger.debug(" group %d has %d nieghbors with radius of %f "%(g_i,nb_cnt,group_i.properties['radius']))
        for g_j in groupset_i.group_nblist.getnbs(g_i):
            if( g_j > g_i ):
                group_j  = groupset_i.groups[g_j]
                if( group_i.mol != group_j.mol ):
                    row_i = [g_i,g_j,group_i.properties['mol'],group_j.properties['mol'] ]
                    pair_writer.writerow(row_i)
  
        fout.close()
        

def setup_calc(proj_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    #
    if( rank == 0 ):
        logging.info('Running setup_calc on %d procs %s '%(size,datetime.now()))        
    #
    if( rank == 0 ):
        logger.info("Setting up %s "%(proj_tag))
    #
    #
    if( rank == 0 ):
        logger.info("Reading %s "%(options.t_nw))
    script_nw_str = string_load(options.t_nw)
    #
    if( rank == 0 ):
        logger.info("Reading %s "%(options.t_run))
    script_run_str = string_load(options.t_run)
    #
    if( rank == 0 ):
        logger.info("Setting up resource %s "%("local"))
    #
    # Set up local as resource 
    local = resource.Resource("local")
    # Set default simulation specs 
    local.specs['exe_command'] = './' 
    local.specs['ppn'] = 24
    local.meta['launch_dir'] = local.meta['home_dir'] 
    # Create local directories for project 
    local.make_dir()
    local.dump_json()

    if( rank == 0 ):
        logger.info("Reading in groups ")
    group_tags,group_list = read_group_cplys(options.group_id)
    if( rank == 0 ):
        logger.info("Reading in pairs ")
    pair_tags = read_pairs(options.group_id)
    
    if( rank == 0 ):
        logger.info("Setting up  simulations for %s  "%(proj_tag))
        
    sims_file = "et_sims.csv"
    if( not file_test(sims_file) ):
        sim_tags = read_sims(sims_file)
    else:
        sim_tags = []
    
    if( rank == 0 ):
        logger.info("1st %s read with %d entries "%(sims_file,len(sim_tags)))
        logger.info("Writing %s header "%(sims_file))
        fout = open(sims_file,'wb')
        sims_writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        header = ['tag']
        #if( rank == 0 ):
        sims_writer.writerow(header)
        for sim_tag in sim_tags:
            sims_writer.writerow([sim_tag])
        fout.close()

        
    N_pairs = len(pair_tags)
    pair_indx_list = range(N_pairs)
    pair_indx_list_p = p.splitListOnProcs(pair_indx_list)
    pair_cnt = 0
    #simtags_p = []
    #sims_p = []
    N_pairs_p = len(pair_indx_list_p)
    for pr_i in pair_indx_list_p:
        pw_start_time = datetime.now()
        
        pairs = pair_tags[pr_i]
        g_i = pairs[0]
        g_j = pairs[1]
        tag_i = "et_%s_%d_%d"%(options.group_id,g_i,g_j)

        if( tag_i not in sim_tags):

        
            g_i_info = group_tags[g_i]
            g_j_info = group_tags[g_j]
            mol_i = int(g_i_info[2])
            mol_j = int(g_j_info[2])
            #if( mol_i != mol_j ):
            pair_cnt += 1         

            sim_i = simulation.NWChem(tag_i)
            sim_i.set_resource(local)
            sim_i.template_nw = copy.deepcopy(script_nw_str)
            sim_i.template_script = copy.deepcopy(script_run_str)
            sim_i.write_et(group_list[g_i],group_list[g_j])
            sim_i.write_script()
            sim_i.dump_json()
            sim_i.push()
            # proj_i.add_sim(sim_i)
            #simtags_p.append(tag_i)
            #sims_p.append(sim_i)

            fout = open(sims_file,'a')
            sims_writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            row_i = [tag_i]
            sims_writer.writerow(row_i)
            fout.close()

            pw_finish_time = datetime.now()
            delt_t = pw_finish_time - pw_start_time
            logger.info(" Wrote %s mol_i %d mol_j %d on proc %d time %f s %d/%d (%d)"%(tag_i,mol_i,mol_j,rank,delt_t.seconds,pair_cnt,N_pairs_p,N_pairs))
        else:
            logger.info("Tag %s already written "%(tag_i))

    p.barrier()
    # Add simulations to project
    # sims_all =  p.gatherList(sims_p)
    sim_tags = read_sims(sims_file)
    if( rank == 0 ):
        logger.info("2nd %s read with %d entries "%(sims_file,len(sim_tags)))


    # 
    if( rank == 0 ):
        logger.info("Writing project %s"%(proj_tag))
        proj_i = project.Project(proj_tag)
        proj_i.add_resource(local)
        proj_i.set_resource(local)
        N_sims = len(sim_tags)
        for sim_i_tag in sim_tags:
            #proj_i.add_sim(sim_i)
            sim_i = simulation.load_json(sim_i_tag)
            ind_j = 0 
            for sim_j in proj_i.simulations:
                if( sim_i_tag == sim_j.tag ):
                    proj_i.simulations.pop(ind_j)
                ind_j += 1 
            proj_i.simulations.append(sim_i)
            logger.info("Simulation %s [%d]/%d  with status %s added to project %s "%(sim_i.tag,sim_cnt,N_sims,sim_i.meta['status'],proj_i.tag))
        proj_i.dump_json()
    p.barrier()
    
def check_p(proj_i,p):
    sims_p = []
    # proj_i.check()
    N_sim = len(proj_i.simulations)
    sim_indx_list = range(N_sim)
    sim_indx_list_p = p.splitListOnProcs(sim_indx_list)    
    sim_cnt = 0 
    for sim_indx in sim_indx_list_p:  # proj_i.simulations: #[0:1]:
        sim_i = proj_i.simulations[sim_indx]
        sim_i.check()
        # print " simulation %s has status %s "%(sim_i.tag,sim_i.meta['status'])
        os.chdir(proj_i.home_dir)
        if( sim_i.meta['status'] != 'finished' ):
            sim_i.meta['status'] = 'written'
            sims_p.append(sim_i)
        sim_i.dump_json()
        logger.info("Simulation %s [%d]  with status %s "%(sim_i.tag,sim_cnt,sim_i.meta['status']))
        sim_cnt += 1
    return sims_p

def split_proj(proj_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
        
    if( rank == 0 ):
        logging.info('Running on %d procs %s '%(size,datetime.now()))        
    
    if( rank == 0 ):
        logger.info("Running split_proj function for  %s "%(proj_tag))
        
    sims_file = "et_sims.csv"
    sim_tags = read_sims(sims_file)
    # sim_tags = sim_tags[0:1000]
    N_sims = len(sim_tags)
    if( rank == 0 ):
        logger.info("3rd %s read with %d entries "%(sims_file,len(sim_tags)))

    if( rank == 0 ):
        logger.info("Setting up resource ")
    # 
    # Set up HPC resource 
    #
    peregrine = resource.Resource("peregrine")
    peregrine.meta['type'] = "local"
    peregrine.meta['username'] = "tkemper"    
    peregrine.meta['address'] = "peregrine.nrel.gov"    
    peregrine.meta['home_dir'] = '/home/%s'%(peregrine.meta['username'])
    peregrine.meta['storage_dir'] = '/mss/users/%s'%(peregrine.meta['username'])
    peregrine.meta['scratch_dir'] = '/scratch/%s/%s'%(peregrine.meta['username'],proj_tag)
    peregrine.meta['launch_dir'] = peregrine.meta['scratch_dir'] 
    # Set default simulation specs 
    peregrine.specs['allocation'] = 'orgopv'
    peregrine.specs['walltime'] = 48
    peregrine.specs['nodes'] = int(1)
    peregrine.specs['ppn'] = int(24)
    peregrine.specs['queue'] = 'batch'
    peregrine.specs['feature'] = '24core'
    peregrine.specs['exe_command'] = 'qsub '
    peregrine.dump_json()
        
    '''
    # 
    # This is way too slow needs to be put on nodes 
    #
    #sys.exit("2093r902u3r982uriu2h498y")
    # check_p(proj_i,p)
    if( rank == 0 ):
        logger.info("Splitting tags onto %d proc "%(size))

    sims_p = []
    sim_tags_p = p.splitListOnProcs(sim_tags)
    for sim_i_tag in sim_tags_p:
            #proj_i.add_sim(sim_i)
            sim_i = simulation.load_json(sim_i_tag)
            # sim_i.push()
            sim_i.check()
            # print " simulation %s has status %s "%(sim_i.tag,sim_i.meta['status'])
            os.chdir(peregrine.meta['scratch_dir'])
            if( sim_i.meta['status'] != 'finished' ):
                sim_i.meta['status'] = 'written'
                sims_p.append(sim_i)
    sims_all =  p.gatherList(sims_p)
    if( rank == 0 ):
        logger.info("%d sims to run "%(N_sims))
    
    # proj_i.run()
    # proj_i.check()
    # sims_p = check_p(proj_i,p)
    # Gather sim list from proc 
    # sims_all =  p.gatherList(sims_p)
    '''
    
    n_nodes = int(N_sims/500)
    if( n_nodes < 1 ):
        n_nodes = 1 
    nodes = range(n_nodes)
    
    if( rank == 0 ):
        logger.info("Splitting sims onto %d nodes  "%(n_nodes))
    nodes_p = p.splitListOnProcs(nodes)
    
    
    n = N_sims/n_nodes
    r = N_sims%n_nodes
    b,e = 0, n + min(1, r) # first split
    newseq = []                       # New partitioned list
    for i in nodes:
        newseq.append(sim_tags[b:e])
        r = max(0, r-1)             # use up remainders
        b,e = e, e + n + min(1, r)  # min(1,r) is always 0 or 1
    plist = newseq                  # Make 'parts' number of chunks
    # Error check or return results
    if len(plist) != n_nodes:
        print " " 
        print "Partitioning failed, check data length and #-procs"
        sys.exit(0)
        
    p.barrier()
    if( rank == 0 ):
        logger.info("Writing projects for each node")
    
    for n in nodes_p:
        simtags_n = plist[n]
        N_sims_p = len(simtags_n)
        proj_tag_n = "%s_node_%d"%(proj_tag,n)
        proj_n = project.Project(proj_tag_n)
        
        proj_n.add_resource(peregrine)
        proj_n.set_resource(peregrine)
        proj_n.meta['scratch_dir'] = peregrine.meta['scratch_dir'] 
        proj_n.specs = peregrine.specs

        logger.info("Writing proj %s with %d sims on proc %d "%(proj_n.tag,N_sims_p,rank))
        sim_cnt = 0 
        for sim_i_tag in simtags_n:
            #proj_i.add_sim(sim_i)
            sim_i = simulation.load_json(sim_i_tag)
            #sim_i.push()
            # sim_i.check()
            # print " simulation %s has status %s "%(sim_i.tag,sim_i.meta['status'])
            # os.chdir(proj_n.meta['scratch_dir'])
            #if( sim_i.meta['status'] != 'finished' ):
            #    sim_i.meta['status'] = 'written'
            proj_n.add_sim(sim_i)
            logger.info("Sim %s [%d]/%d (%d) with status %s => %s "%(sim_i.tag,sim_cnt,N_sims_p,N_sims,sim_i.meta['status'],proj_n.tag))
            #else:
            #    logger.info("Sim %s status %s"%(sim_i.tag,sim_i.meta['status']))
            # sim_i.dump_json()
            sim_cnt += 1
        logger.info("Sim %s on proc %d for node %d "%(proj_tag_n,rank,n))
        proj_n.files['input']['pyscript'] = 'run_proj.py'
        proj_n.template_script = string_load('run_etproj_v1.pbs')
        proj_n.write_script()
        proj_n.dump_json()
    
    p.barrier()
    
def run_calc(proj_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()

    
        
    if( rank == 0 ):
        logging.info('Running on %d procs %s '%(size,datetime.now()))        
    
    if( rank == 0 ):
        logger.info("Running  %s "%(proj_tag))
        
    proj_i = project.load_json(calc_tag)
    local = proj_i.resources["local"]

    #sys.exit("2093r902u3r982uriu2h498y")

    proj_i.run()
    # proj_i.check()
    proj_i = check_p(proj_i)
    #proj_i.store()
    # proj_i.analysis()
    os.chdir(proj_i.home_dir)
    proj_i.dump_json()


def read_energies(proj_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()

    if( rank == 0 ):
        logging.info('Running on %d procs %s '%(size,datetime.now()))        

    # Run position update from individual minimizations 
    proj_i = project.load_json(calc_tag)
    proj_i.check()

    group_rows = []
    sim_cnt = 0 
    for sim_i in proj_i.simulations:
        #os.chdir(sim_i.meta['scratch_dir'])
        sim_i.analysis()
        # Get the last run
        run_i = sim_i.run_list[-1]
        # Get the last energy 
        en_i = run_i.toteng_list[-1]
        group_rows.append([sim_cnt,sim_i.tag,en_i])
        
        sim_cnt += 1

    print "Changing to home directory %s "%(proj_i.home_dir)
    os.chdir(proj_i.home_dir)


    group_file = "group_%s.csv"%(options.group_id)
    fout = open(group_file,'wb')
    group_writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
    header = ['g_i','tag','energy']
    #if( rank == 0 ):
    group_writer.writerow(header)
    for row in group_rows:
        group_writer.writerow(row)
    fout.close()
        
    if( rank == 0 ):
        logger.info('output group_%s group_%s.csv'%(options.group_id,options.group_id))
    proj_i.dump_json()

def read_struc(proj_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()

    if( rank == 0 ):
        logging.info('Running on %d procs %s '%(size,datetime.now()))        

    if( rank == 0 ):
        logger.info("Creating empty buildingblock.Container  %s "%('mol_min'))
    
    mol_lmp = buildingblock.Container('mol_min')

    # Run position update from individual minimizations 
    proj_i = project.load_json(calc_tag)
    proj_i.check()

    sim_cnt = 0 
    for sim_i in proj_i.simulations:
        os.chdir(sim_i.meta['scratch_dir'])
        sim_i.read_data(sim_i.files['output']["1"])
        mol_lmp += sim_i.strucC
        print "From sim %s reading data file %s with %d particles into container with %d particles "%(sim_i.tag,sim_i.files['output']["1"],sim_i.strucC.n_particles,mol_lmp.n_particles)
        #if( rank == 0 ):
        #    logger.info("output sim_%d_json %s "%(group_tag,group_file))
        sim_cnt += 1

    if( rank == 0 ):
        logger.info("Reading in structure from %s "%(options.cply))

    if( rank == 0 ):
        logger.info("Reading  %s "%(proj_tag))
        

    print "Changing to home directory %s "%(proj_i.home_dir)
    os.chdir(proj_i.home_dir)
    new_bulk = buildingblock.Container(proj_tag)
    new_bulk.read_cply(cply_file = options.cply)
    new_bulk.write_xyz(xmol_file='initial.xyz')
    new_bulk.positions = mol_lmp.positions
    if( rank == 0 ):
        logger.info("Groups for pbcs %s "%(options.group_id))
    new_bulk.particles_select = new_bulk.particles.keys()
    mols_index = new_bulk.group_prop("mol")
    new_bulk.group_pbcs()

    if( rank == 0 ):
        logger.info("output cply %s.cply "%(new_bulk.tag))
        
    new_bulk.write_cply()
    new_bulk.write_xyz(xmol_file='finial.xyz')
    # 
    proj_i.add_struc(new_bulk)
    proj_i.dump_json()


def write_newdata(calc_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    
    if( rank == 0 ):
        logger.info("Reading %s "%(options.t_in))
        script_in_str = string_load(options.t_in)
        
    if( rank == 0 ):
        logger.info("Reading %s "%(options.t_run))
        script_run_str = string_load(options.t_run)
    
    proj_i = project.load_json(calc_tag)

    if( rank == 0 ):
        logger.info("Writing new data file ")
        
    local = proj_i.resources["local"]

    lmp_new = simulation.LAMMPS(calc_tag)
    lmp_new.set_resource(local)
    lmp_new.read_data(options.data)
    # lmp_new.read_data('prod1.data')
    struc_o = buildingblock.Container(calc_tag)
    struc_o.read_cply()
    lmp_new.strucC.positions = struc_o.positions
    # Update mols numbers 
    for pkey_i, particle_i in struc_o.particles.iteritems():
        lmp_new.strucC.particles[pkey_i].properties["mol"] = particle_i.properties["mol"]
    #struc_o = lmp_new.update_cply('%s.cply'%(calc_tag))
    lmp_new.write_data()
    if( rank == 0 ):
        logger.info("output data %s.data "%(calc_tag))

    lmp_new.template_in = copy.deepcopy(script_in_str)
    lmp_new.template_script = copy.deepcopy(script_run_str)
    lmp_new.write_data()

    '''
    lmp_new.write_in()
    lmp_new.write_script()
    lmp_new.push()
    
    proj_i.add_sim(lmp_new)
    proj_i.dump_json()
    
    lmp_new.run()
    lmp_new.check()
    #proj_i.store()
    lmp_new.analysis()
    
    proj_i.dump_json()
    '''
def et(calc_tag,options,p):

    #rdf(calc_tag,options)
    group_file = 'group_%s.csv'%(options.group_id)
    pairs_file = 'pairs_%s.csv'%(options.group_id)
    if( file_test(group_file) and file_test(pairs_file)  ):
        pull_groups(calc_tag,options,p)
    elif( rank == 0 ):
        logger.info('file: output group_%s %s '%(options.group_id,group_file))
        logger.info('file: output pairs_%s %s '%(options.group_id,pairs_file))
       
    #sims_file = "et_sims.csv"
    #if( file_test(sims_file) ):
    #    setup_calc(calc_tag,options,p)
    #else:
    #    if( rank == 0 ):
    #        logger.info('output sim_file %s'%(sims_file))
    #
    #split_proj(calc_tag,options,p)
    #read_energies(calc_tag,options,p)
    #read_struc(calc_tag,options,p) 
    #write_newdata(calc_tag,options,p)
    #run_calc(calc_tag,options,p)    

if __name__=="__main__":
    
    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--cply", dest="cply", type="string", default="", help="Input cply file")
    parser.add_option("--list_i", dest="list_i", type="string", default="", help="Input list file")
    parser.add_option("--t_nw", dest="t_nw", type="string", default="", help="Template nw file")
    parser.add_option("--t_run", dest="t_run", type="string", default="", help="Template run script file")
    parser.add_option("--group_id", dest="group_id", type="string", default="mol", help="Group id ")
    parser.add_option("--hterm", dest="hterm", default=False,action="store_true", help=" Hydrogen terminate groups ")
    
    (options, args) = parser.parse_args()
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if( len(args) < 1 ):
        calc_tag = 'et'
    else:
        calc_tag =  args[0]

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)


    hdlr = logging.FileHandler('%s.log'%(calc_tag),mode='w')
    hdlr.setLevel(logging.INFO)
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)

    start_time = datetime.now()
    if( rank == 0 ):
        logger.info('Started %s '%(start_time))

    et(calc_tag,options,p)

    finish_time = datetime.now()
    delt_t = finish_time - start_time
    if( rank == 0 ):
        logger.info('Finished %s in %f seconds '%(finish_time,delt_t.seconds))
