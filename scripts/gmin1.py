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

def read_calcs(sims_file):
    
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
                        # print ">read_calcs sim_tags",str(row[0])
                    rownum += 1
                
    return sim_tags
        
def set_res(proj_tag):

    # Set up local as resource 
    local = resource.Resource("local")
    # Set default simulation specs 
    local.properties['exe_command'] = './' 
    local.properties['ppn'] = 24
    local.properties['nproc'] = 24
    local.dir['launch'] = local.dir['home'] 
    # Create local directories for project 
    local.make_dir()
    local.dump_json()

    # 
    # Set up HPC resource 
    #
    peregrine = resource.Resource("peregrine")
    peregrine.meta['type'] = "local"
    peregrine.ssh['username'] = "tkemper"    
    peregrine.ssh['address'] = "peregrine.nrel.gov"    
    peregrine.dir['storage'] = '/mss/users/%s'%(peregrine.ssh['username'])
    peregrine.dir['scratch'] = '/scratch/%s/%s'%(peregrine.ssh['username'],proj_tag)
    peregrine.dir['home'] = peregrine.dir['scratch']  #'/home/%s'%(peregrine.ssh['username'])
    peregrine.dir['launch'] = peregrine.dir['scratch'] 
    # Set default simulation specs 
    peregrine.properties['allocation'] = 'orgopv'
    peregrine.properties['walltime'] = 48
    peregrine.properties['nodes'] = int(1)
    peregrine.properties['ppn'] = int(24)
    peregrine.properties['nproc'] = 24
    peregrine.properties['queue'] = 'batch'
    peregrine.properties['feature'] = '24core'
    peregrine.properties['exe_command'] = 'qsub '
    peregrine.properties['e-mail'] = 'travis.kemper@nrel.gov'
    peregrine.dump_json()
 
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
    if( len(options.list_i) > 0 ):
        lfile = open(options.list_i,'rb')
        t_i = lfile.readlines()
        lfile.close()
        list_i = [int(pkey) for pkey in  t_i]
    else:
        # Leave list blank and use all the particles 
        list_i = []
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
        for x_i in group_i.properties['cent_mass']:
            row_i.append(x_i)
        group_writer.writerow(row_i)
        fout.close()

def setup_calc(proj_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    #
    if( rank == 0 ):
        logger.info("Reading in groups ")
    group_tags,group_list = read_group_cplys(options.group_id)
        
    #
    if( rank == 0 ):
        logging.info('Running setup_calc on %d procs %s '%(size,datetime.now()))        
        logger.info("Setting up %s "%(proj_tag))
        logger.info("Reading %s "%(options.t_in))
        logger.info("Reading %s "%(options.t_run))
        logger.info("Setting up resource %s "%("local"))
        
    local = resource.Resource('local')
    local.load_json()
    #
    lmp_o = lammps.LAMMPS('lmp_ref')
    lmp_o.files['templates']['t_in'] = options.t_in
    lmp_o.files['templates']['t_run'] = options.t_run

    lmp_o.load_str('templates','t_in')
    lmp_o.load_str('templates','t_run')
    lmp_o.read_param(options.param)

    if( rank == 0 ):
        logger.info("Setting up  simulations for %s  "%(proj_tag))
    sims_file = "%s.csv"%(calc_tag)
    if( not file_test(sims_file) ):
        sim_tags = read_calcs(sims_file)
    else:
        sim_tags = []

    if( rank == 0 ):
        logger.info("1st %s read with %d entries "%(sims_file,len(sim_tags)))
        logger.info("Writing %s header "%(sims_file))
        fout = open(sims_file,'wb')
        sims_writer = csv.writer(fout,delimiter=',')
        header = ['tag']
        #if( rank == 0 ):
        sims_writer.writerow(header)
        for sim_tag in sim_tags:
            sims_writer.writerow([sim_tag])
        fout.close()
        
    N_groups = len(group_tags)
    group_indx_list = range(N_groups)
    group_indx_list_p = p.splitListOnProcs(group_indx_list)
    group_cnt = 0
    for g_i in group_indx_list_p:
        gw_start_time = datetime.now()
        
        tag_i = "gmin_%s_%d"%(options.group_id,g_i)
        if( tag_i not in sim_tags):

            lmp_i = lammps.LAMMPS(tag_i)        
            lmp_i.set_resource(local)
            lmp_i.properties['scratch'] = lmp_i.dir['scratch']
            
            lmp_i.str['t_in'] = copy.deepcopy(lmp_o.str['t_in'])
            lmp_i.str['t_run'] = copy.deepcopy(lmp_o.str['t_run'])
            bb_i  = group_list[g_i]
            
            lmp_i.add_strucC(bb_i)
            
            #lmp_i.make_dir()
            if ( not os.path.isdir(lmp_i.dir['scratch']) ):
                print "Making %s "%(lmp_i.dir['scratch'])
                os.mkdir(lmp_i.dir['scratch'])
                         
            os.chdir(lmp_i.dir['scratch'])

            lmp_i.paramC = copy.deepcopy(lmp_o.paramC)
            lmp_i.strucC.bonded_bonds()
            lmp_i.strucC.bonded_angles()
            lmp_i.strucC.bonded_dih()
            lmp_i.set_ffparam()
            lmp_i.write_data()
            
            lmp_i.replacewrite_prop('t_in','input','in',"%s.in"%(lmp_i.tag))
            lmp_i.properties['input_in'] = lmp_i.files['input']['in']
            lmp_i.replacewrite_prop('t_run','scripts','run',"%s.sh"%(lmp_i.tag))
            
            lmp_i.add_file('output','log',"%s.log"%(lmp_i.tag))

            os.chdir(lmp_i.dir['home'])
            lmp_i.dump_json()
            
            fout = open(sims_file,'a')
            sims_writer = csv.writer(fout,delimiter=',')
            row_i = [tag_i]
            sims_writer.writerow(row_i)
            fout.close()

            gw_finish_time = datetime.now()
            delt_t = gw_finish_time - gw_finish_time
            logger.info(" Wrote %s g_i %d on proc %d time %f s "%(tag_i,g_i,rank,delt_t.seconds))
        else:
            logger.info("Tag %s already written "%(tag_i))
        

def split_proj(calc_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    res_i = resource.Resource(options.proj_res)
    res_i.load_json()
    if( rank == 0 ):
        logging.info('Running on %d procs %s '%(size,datetime.now()))        
        logger.info("Running split_proj function for  %s "%(calc_tag))
        
    calcs_file = "%s.csv"%(calc_tag)
    calc_tags = read_calcs(calcs_file)
    # calc_tags = calc_tags[0:1000]
    N_calcs = len(calc_tags)
    if( rank == 0 ):
        logger.info("3rd %s read with %d entries "%(calcs_file,len(calc_tags)))

    if( rank == 0 ):
        logger.info("Setting up resource ")
    #
    n_jobs = int(N_calcs/options.calcs_node)
    if( n_jobs < 1 ):
        n_jobs = 1 
    jobs = range(n_jobs)
    
    if( rank == 0 ):
        logger.info("Splitting calcs onto %d jobs  "%(n_jobs))
    jobs_p = p.splitListOnProcs(jobs)
    
    n = N_calcs/n_jobs
    r = N_calcs%n_jobs
    b,e = 0, n + min(1, r) # first split
    newseq = []                       # New partitioned list
    for i in jobs:
        newseq.append(calc_tags[b:e])
        r = max(0, r-1)             # use up remainders
        b,e = e, e + n + min(1, r)  # min(1,r) is always 0 or 1
    plist = newseq                  # Make 'parts' number of chunks
    # Error check or return results
    if len(plist) != n_jobs:
        print " " 
        print "Partitioning failed, check data length and #-procs"
        sys.exit(0)
        
    p.barrier()
    if( rank == 0 ):
        logger.info("Writing projects for each node")
        
    proj_o = project.Project(calc_tag)
    proj_o.set_resource(res_i)
        
    
    for n in jobs_p:
        calctags_n = plist[n]
        N_calcs_p = len(calctags_n)
        calc_tag_n = "%s_node_%d"%(calc_tag,n)
        proj_n = project.Project(calc_tag_n)
        proj_n.set_resource(res_i)
        
        logger.info("Writing proj %s with %d calcs on proc %d "%(proj_n.tag,N_calcs_p,rank))
        calc_cnt = 0 
        for calc_i_tag in calctags_n:
            calc_i = lammps.LAMMPS(calc_i_tag)
            calc_i.load_json()
            print " calculation %s has status %s "%(calc_i.tag,calc_i.meta['status'])
            # os.chdir(proj_n.meta['scratch_dir'])
            #if( calc_i.meta['status'] != 'finished' ):
            #    calc_i.meta['status'] = 'written'
            #proj_i.add_calc(calc_i)
            proj_n.calculations[calc_i_tag] = calc_i
            logger.info("Calc %s [%d]/%d (%d) with status %s => %s "%(calc_i.tag,calc_cnt,N_calcs_p,N_calcs,calc_i.meta['status'],proj_n.tag))
            #else:
            #    logger.info("Calc %s status %s"%(calc_i.tag,calc_i.meta['status']))
            # calc_i.dump_json()
            calc_cnt += 1
        logger.info("Calc %s on proc %d for node %d "%(calc_tag_n,rank,n))
        proj_n.files['input']['pyscript'] = options.proj_in
        proj_n.files['templates']['run'] = options.proj_run
        proj_n.load_str('templates','run')
        proj_n.properties['streamm_command'] = 'python %s %s > %s.out '%(options.proj_in,proj_n.tag,proj_n.tag)
        proj_n.replacewrite_prop('run','input','run',"%s.pbs"%(proj_n.tag))
        proj_n.dump_json()
        # Add node project to global project 
        proj_o.calculations[proj_n.tag] = proj_n
    
    p.barrier()
    proj_o.dump_json()
    
    proj_file = "proj_%s.json"%(calc_tag)
    if( rank == 0 ):
        logger.info('file: output proj_file %s'%(proj_file))
    

def proj_run(calc_tag,options,p):
    proj_i = project.Project(calc_tag)
    proj_i.load_json()
    proj_i.check()
    proj_i.run()
    os.chdir(proj_i.dir['home'])
    
def proj_analysis(calc_tag,options,p):
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()
    #
    proj_i = project.Project(calc_tag)
    proj_i.load_json()
    #proj_i.check()

    gmin = buildingblock.Container('gmin')
    
    # proj_i.analysis()
    for proj_key,proj_j in proj_i.calculations.iteritems():
        print "Running analysis on project %s "%(proj_key)
        for calc_key,calc_i in proj_j.calculations.iteritems():
            os.chdir(calc_i.dir['scratch'])
            calc_i.analysis()
            #calc_i.read_data(calc_i.files['output']["data_1"])
            gmin += calc_i.strucC
            
            os.chdir(calc_i.dir['home'])
            calc_i.dump_json()

    os.chdir(proj_i.dir['home'])
    os.getcwd()
    logger.info("Reading in original structure %s "%(options.cply))
    struc_o = buildingblock.Container('%s'%(calc_tag))
    struc_o.read_cply(options.cply)
    struc_o.write_xyz(xyz_file='initial.xyz')
    proj_i.files['output']['initial_xyz'] = 'initial.xyz'
    # Update positions from group calcs
    if( struc_o.n_particles == gmin.n_particles ):
        struc_o.positions = gmin.positions
    else:
        logger.warning("Initial structure has  %d particles and updated group calc structure has %d "%( struc_o.n_particles , gmin.n_particles))
        sys.exit() 
    #
    # Find groups 
    #
    logger.info(" Applying %s group pbcs "%(options.group_id))
    struc_o.group_prop(options.group_id,options.group_id)
    groupset_i = struc_o.groupsets[options.group_id]
    groupset_i.group_pbcs()
    struc_o.write_xyz()
    proj_i.files['output']['xyz'] = '%s.xyz'%(calc_tag)
    if( rank == 0 ):
        logger.info("file: output xyz %s "%(proj_i.files['output']['xyz']))
    struc_o.write_cply()
    proj_i.files['output']['cply'] = '%s.cply'%(calc_tag)
    if( rank == 0 ):
        logger.info("file: output cply %s "%(proj_i.files['output']['cply']))
    #
    proj_i.dump_json()
                                
def gmin(calc_tag,options,p):

    set_res(calc_tag)
    group_file = 'group_%s.csv'%(options.group_id)
    
    if( file_test(group_file)  ):
        pull_groups(calc_tag,options,p)
    elif( rank == 0 ):
        logger.info('file: output group_%s %s '%(options.group_id,group_file))

    calcs_file = "%s.csv"%(calc_tag)
    if( file_test(calcs_file) ):
        setup_calc(calc_tag,options,p)
    elif( rank == 0 ):
        logger.info('file: output calc_file %s'%(calcs_file))
              
    proj_file = "proj_%s.json"%(calc_tag)
    if( file_test(proj_file) ):
        split_proj(calc_tag,options,p)
    elif( rank == 0 ):
        logger.info('file: output proj_file %s'%(proj_file))

    proj_run(calc_tag,options,p)    
    proj_analysis(calc_tag,options,p)    
    
if __name__=="__main__":
    
    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--cply", dest="cply", type="string", default="", help="Input cply file")
    parser.add_option("--list_i", dest="list_i", type="string", default="", help="Input list file")
    parser.add_option("--itp", dest="itp", type="string", default="", help="Input itp file")
    parser.add_option("--param", dest="param", type="string", default="", help="Input param file")
    parser.add_option("--t_in", dest="t_in", type="string", default="", help="Template in file")
    parser.add_option("--t_run", dest="t_run", type="string", default="", help="Template run script file")
    parser.add_option("--group_id", dest="group_id", type="string", default="mol", help="Group id ")
    parser.add_option("--hterm", dest="hterm", default=False,action="store_true", help=" Hydrogen terminate groups ")
    parser.add_option("--calcs_node", dest="calcs_node", type="int", default=10, help="Number of calculations per node ")
    parser.add_option("--proj_res", dest="proj_res", type="string", default="local", help="Project resource tag ")
    parser.add_option("--proj_in", dest="proj_in", type="string", default="run_proj.py", help="Project python input tag ")
    parser.add_option("--proj_run", dest="proj_run", type="string", default="streamm.sh", help="Project run script tag ")
    
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
    logger.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if( len(args) < 1 ):
        calc_tag = 'gmin'
    else:
        calc_tag =  args[0]

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
    
    gmin(calc_tag,options,p)
    
    finish_time = datetime.now()
    delt_t = finish_time - start_time
    
    logger.info('Finished %s in %f seconds '%(finish_time,delt_t.seconds))
