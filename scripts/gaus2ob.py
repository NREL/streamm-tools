'''
Set of luigi tasks for the opv project

version 0.3

'''
import   os, os.path, sys , copy ,shutil, logging, math, json, csv 
import numpy as np
from datetime import datetime
from optparse import OptionParser

import pybel,openbabel

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
    
def convert_pos_bohr_ang(struc_o):
    for pkey_o  in struc_o.particles.keys():
        pos_o = struc_o.positions[pkey_o]       
        pos_i = [ units.convert_bohr_ang(v_i) for v_i in pos_o ] 
        struc_o.positions[pkey_o] = pos_i
        
def opt_struc(options):

    # Set up resources
    set_res(options)

    tag_o = options.tag
    bbdir = options.bbdir
    bbtype = options.bbtype

    local = resource.Resource('local')
    local.load_json()
    peregrine = resource.Resource('peregrine')
    peregrine.load_json()

    # Copy cply to scratch dir
    gaus_i = gaussian.Gaussian("%s_tddft_ob"%tag_o)
    gaus_i.set_resource(peregrine)
    gaus_i.dir['BuildingBlocks'] = "%s/%s"%(bbdir,bbtype)
    gaus_i.dir['scratch'] = '/scratch/%s/%s'%(peregrine.ssh['username'],gaus_i.tag)

    gaus_i.files['output']['log'] = "%s.log"%(tag_o)
    gaus_i.files['output']['fchk'] = "%s.fchk"%(tag_o)

    gaus_i.struc_o = buildingblock.Container(tag_o)
    gaus_i.struc_o.read_cply()
    gaus_i.files['input']['cply'] = "%s.cply"%(tag_o)
    gaus_i.struc_o.lat_cubic(1000.0)
    gaus_i.struc_o.write_xyz()
    gaus_i.files['input']['xyz'] = "%s.xyz"%(tag_o)

    gaus_i.strucC = copy.deepcopy(gaus_i.struc_o)
    gaus_i.strucC.tag  = tag_o

    gaus_i.check()
    gaus_i.analysis()

    convert_pos_bohr_ang(gaus_i.strucC)

    print gaus_i.strucC.n_bonds
    gaus_i.strucC.lat_cubic(1000.0)
    gaus_i.strucC.write_xyz()    
    gaus_i.strucC.write_cply()
    gaus_i.files['output']['xyz'] = "%s.xyz"%(gaus_i.strucC.tag)

    gaus_i.cp_file('output','cply',"%s.cply"%(gaus_i.strucC.tag),'scratch','BuildingBlocks')
    # gaus_i.files['output']['cply'] = "%s.cply"%(gaus_i.strucC.tag)

    #gaus_j.strucC.bonded_nblist.guess_nblist(gaus_j.strucC.lat,gaus_j.strucC.particles,gaus_j.strucC.positions,"cov_radii",radii_buffer=1.25)
    #gaus_j.strucC.bonded_bonds()
    xyz_str = gaus_i.strucC.write_xyz_str()
    pybelmol = pybel.readstring('xyz',xyz_str)

    desc = pybelmol.calcdesc(['atoms','bonds','HBA1','HBA2','MW','MR','logP','TPSA','HBD','sbonds','abonds','tbonds','dbonds'])
    for pkey,prop_val in  desc.iteritems():
        gaus_i.properties[pkey] = prop_val

    molstr  = pybelmol.write('mol')
    gaus_i.properties['molstr'] = molstr

    mass = pybelmol.exactmass
    gaus_i.properties['mass'] = mass

    nrings = 0
    rings = []
    for ring in openbabel.OBMolRingIter(pybelmol.OBMol):
        nrings = nrings + 1
        rings.append(ring.Size())
    gaus_i.properties['nrings'] = nrings
    ring_size = ",".join(map(str,rings))
    gaus_i.properties['rings'] = ring_size
    nrbonds = pybelmol.OBMol.NumRotors()
    gaus_i.properties['nrbonds'] = nrbonds


    # Calc fingerprint  
    fpbits = pybelmol.calcfp().bits
    if( len(fpbits) > 0 ):
        fplong = reduce(lambda x,y: x|y,map(lambda x: 1<<(x-1),fpbits))
    else:
        fpbits = 001
        fplong = fpbits
        
    gaus_i.properties['fpbits'] = fpbits
    gaus_i.properties['fplong'] = fplong

    # Do not calc fp4 it will remove hydrogens 
    #fp4bits = pybelmol.calcfp(fptype='FP4').bits
    #gaus_i.properties['fp4bits'] = fp4bits
    #fp4long = reduce(lambda x,y: x|y,map(lambda x: 1<<(x-1),fp4bits))
    #gaus_i.properties['fp4long'] = fp4long

    output = pybel.Outputfile("sdf", "%s.sdf"%(gaus_i.tag))
    output.write(pybelmol)
    output.close()

    gaus_i.dump_json()
    gaus_i.cp_file('output','json',"%s_%s.json"%(gaus_i.prefix,gaus_i.tag),'scratch','BuildingBlocks')
    # gaus_i.meta['status'] = 'finished'
    gaus_i.store()


if __name__=="__main__":

    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--tag", dest="tag", type="string", default="BDT", help="Input tag")
    parser.add_option("--bbdir", dest="bbdir", type="string", default="/home/tkemper/BuildingBlocks-cand_opt", help="BuildingBlocks dir")
    parser.add_option("--bbtype", dest="bbtype", type="string", default="donors", help="Input tag")
    parser.add_option("--ppn", dest="ppn", type="int", default=24, help="Processors per node for nwchem calculations ")
    # bbtype = 'donors'
    # bbdir = "/home/tkemper/BuildingBlocks-cand"
    # tdir="/home/tkemper/streamm-tools/templates"
    (options, args) = parser.parse_args()


    if( len(args) < 1 ):
        calc_tag = 'BDT_tddft_ob'
    else:
        calc_tag =  args[0]


    hdlr = logging.FileHandler('%s.log'%(calc_tag),mode='w')
    hdlr.setLevel(logging.DEBUG)
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)


    opt_struc(options)
     
     
