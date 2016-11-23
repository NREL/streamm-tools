'''
Set of luigi tasks for the opv project

version 0.3

'''
import   os, os.path, sys , copy ,shutil, logging, math, json, csv 
import numpy as np
from datetime import datetime
from optparse import OptionParser

# import pybel

from streamm import *

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)


formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

def set_res(options):

    # Set up local as resource 
    local = resource.Resource("local")
    # Set default simulation specs 
    local.properties['exe_command'] = './' 
    local.properties['ppn'] = options.ppn
    local.properties['nproc'] =  options.ppn
    local.properties['feature'] = '%dcore'%( options.ppn)
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
    peregrine.dir['home'] = '/home/%s'%(peregrine.ssh['username'])
    peregrine.dir['storage'] = '/mss/users/%s'%(peregrine.ssh['username'])
    peregrine.dir['scratch'] = '/scratch/%s/%s'%(peregrine.ssh['username'],options.tag)
    peregrine.dir['launch'] = peregrine.dir['scratch'] 
    # Set default simulation specs 
    peregrine.properties['allocation'] = 'orgopv'
    peregrine.properties['walltime'] = 48
    peregrine.properties['nodes'] = int(1)
    peregrine.properties['ppn'] = int( options.ppn)
    peregrine.properties['nproc'] =  options.ppn
    peregrine.properties['queue'] = 'batch'
    peregrine.properties['feature'] = '%dcore'%( options.ppn)
    peregrine.properties['exe_command'] = 'qsub '
    peregrine.properties['e-mail'] = 'travis.kemper@nrel.gov'
    peregrine.dump_json()
    
def opt_struc(options):

    # Set up resources
    set_res(options)

    tag_i = options.tag
    bbdir = options.bbdir
    bbtype = options.bbtype
    tdir = options.tdir
    bbtype = options.bbtype
    bbtype = options.bbtype
    basis = '6-31G(d)'
    method = 'b3lyp'
    nstates = 12
    
    local = resource.Resource('local')
    local.load_json()
    peregrine = resource.Resource('peregrine')
    peregrine.load_json()

    # Copy cply to scratch dir
    gaus_i = gaussian.Gaussian(tag_i)
    gaus_i.set_resource(local)
    gaus_i.dir['BuildingBlocks'] = "%s/%s"%(bbdir,bbtype)
    gaus_i.dir['scratch'] = local.dir['home']
    gaus_i.properties['scratch']  = gaus_i.dir['scratch'] 
    gaus_i.dir['templates'] = tdir
    gaus_i.cp_file('input','cply','%s.cply'%gaus_i.tag,'BuildingBlocks','scratch')
    gaus_i.cp_file('templates','run','gaussian_peregrine.sh','templates','scratch')
    gaus_i.cp_file('templates','com','gaussian_tddft.com','templates','scratch')
    gaus_i.load_str('templates','com')
    gaus_i.load_str('templates','run')
    gaus_i.strucC.tag = tag_i
    gaus_i.strucC.read_cply()
    # gaus_i.properties['commands'] = 'HF/3-21G OPT'
    gaus_i.properties['method'] = method
    gaus_i.properties['basis'] = basis
    gaus_i.properties['nstates'] = nstates
    gaus_i.properties['charge'] = 0
    gaus_i.properties['spin_mult'] = 1
    gaus_i.properties['coord'] = gaus_i.strucC.write_coord()
    gaus_i.replacewrite_prop('com','input','com','%s.com'%(gaus_i.tag))
    gaus_i.properties['input_com'] = gaus_i.files['input']['com']
    gaus_i.replacewrite_prop('run','scripts','run','%s.sh'%(gaus_i.tag))
    gaus_i.files['output']['log'] = "%s.log"%(gaus_i.tag)
    gaus_i.files['output']['fchk'] = "%s.fchk"%(gaus_i.tag)
    
    # gaus_i.check()
    gaus_i.run()
    gaus_i.check()
    gaus_i.analysis()
    gaus_i.dump_json()
    

if __name__=="__main__":
    
    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--ppn", dest="ppn", type="int", default=24, help="Processors per node for nwchem calculations ")
    parser.add_option("--bbdir", dest="bbdir", type="string", default="/home/tkemper/BuildingBlocks-cand", help="BuildingBlocks dir")
    parser.add_option("--tdir", dest="tdir", type="string", default="/home/tkemper/streamm-tools/templates", help="Template dir")
    parser.add_option("--tag", dest="tag", type="string", default="BDT", help="Input tag")
    parser.add_option("--bbtype", dest="bbtype", type="string", default="donors", help="Input tag")

    # bbtype = 'donors'
    # bbdir = "/home/tkemper/BuildingBlocks-cand"
    # tdir="/home/tkemper/streamm-tools/templates"
    
    (options, args) = parser.parse_args()


    if( len(args) < 1 ):
        calc_tag = 'et'
    else:
        calc_tag =  args[0]


    hdlr = logging.FileHandler('%s.log'%(calc_tag),mode='w')
    hdlr.setLevel(logging.DEBUG)
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)


    opt_struc(options)
     
     
