"""
This script adds canidate building blocks to the opv database 

And depends on the streamm v0.3 module found at
https://github.com/NREL/streamm-tools/tree/devel_v0.3.0.12

"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

import sys,  datetime, os , string, copy, math, string, csv, json, logging
import numpy as np
from datetime import datetime
from optparse import OptionParser
import tempfile, psycopg2
import openbabel,pybel
#from mpi4py import MPI
#import pandas as pd

from streamm import (mpiBase,resource,structure,buildingblock,calculation,gaussian,lammps,units)
from streamm import periodictable
#import functionalize as func 

logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


def convert_struc_pymol(bb_i):
    # Get mol str     
    xyz_str = bb_i.write_xyz_str()
    pybelmol = pybel.readstring('xyz',xyz_str)
    molstr  = pybelmol.write('mol')
    # properties = ob_properties(pybelmol)
    desc = pybelmol.calcdesc(['atoms','bonds'])    
    # Check properties from openbabel
    if( bb_i.n_particles != desc['atoms'] ):
            print "%s has an Error in atoms "%(bb_i.properties['ctag'])
            sys.exit(1)
    for p_index in range(pybelmol.OBMol.NumAtoms()):
        p_i = pybelmol.OBMol.GetAtom(p_index+1)
        particle_i = bb_i.particles[p_index]
        if( particle_i.properties['number'] != p_i.GetAtomicNum() ):
            print " ",particle_i.properties['symbol'],particle_i.properties['number'] , p_i.GetAtomicNum()
            print "%s has an Error in particle_i  number "%(bb_i.properties['ctag'])
        #
        particle_i.properties['AtomicMass'] = p_i.GetAtomicMass()
        particle_i.properties['ExactMass'] = p_i.GetExactMass()
        particle_i.properties['FormalCharge'] = p_i.GetFormalCharge()
        particle_i.properties['Isotope'] = p_i.GetIsotope()
        particle_i.properties['SpinMultiplicity'] = p_i.GetSpinMultiplicity()
        particle_i.properties['Valence'] = p_i.GetValence()
        particle_i.properties['Hyb'] = p_i.GetHyb()
        particle_i.properties['ImplicitValence'] = p_i.GetImplicitValence()
        particle_i.properties['HvyValence'] = p_i.GetHvyValence()
        particle_i.properties['HeteroValence'] = p_i.GetHeteroValence()
        #
        # print p_i.GetAtomicNum(),p_i.GetFormalCharge(),p_i.GetIsotope(),p_i.GetSpinMultiplicity(),p_i.GetValence(),p_i.GetHyb(),p_i.GetImplicitValence(),p_i.GetHvyValence(),p_i.GetHeteroValence()

    set_bonds = True
    if( set_bonds ):
        # This segfaults when trying to delete bonds 
        if( bb_i.n_bonds > 0 ):
            b_range = range(1,pybelmol.OBMol.NumBonds()+1)
            for b_index in b_range:
                b_index = 0 
                deleted = pybelmol.OBMol.DeleteBond(pybelmol.OBMol.GetBond(b_index))
            for bkey_i, bond_i  in bb_i.bonds.iteritems():
                # print bond_i.pkey1,bond_i.pkey2,bond_i.properties['bondorder']
                #
                b_i = bond_i.pkey1 + 1 
                b_j = bond_i.pkey2 + 1
                #
                b_n = bond_i.properties['bondorder'] 
                # 1,2,3 - (single,double,triple) bond
                if( b_i < b_j ):
                    pybelmol.OBMol.AddBond(b_i, b_j, b_n)   # atoms indexed from 1
    # 
    desc = pybelmol.calcdesc(['atoms','bonds'])    
    if( bb_i.n_bonds != desc['bonds'] ):
            print "%s has an Error in bonds "%(bb_i.properties['ctag'])
            print bb_i.n_bonds , desc['bonds']
            for b_index in range(pybelmol.OBMol.NumBonds()):
                b_i = pybelmol.OBMol.GetBond(b_index)
                b1 = b_i.GetBeginAtomIdx() -1
                b2 = b_i.GetEndAtomIdx()-1
                bo = b_i.GetBondOrder()
                print ">NumBonds ",b1,b2,bo
                #
            bb_i.write_xyz()
            #fout = open('bond_errors.csv','a')
            #writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            #writer.writerow( [bb_i.tag] )
            #fout.close()            
            sys.exit(2)
                       
    bb_i.properties['molstr'] = molstr
    #
    #
    return pybelmol


def set_res(calc_tag,options):
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
    peregrine.dir['scratch'] = '/scratch/%s/%s'%(peregrine.ssh['username'],calc_tag)
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
    

def countcplys(bbdir,moleculetype):
    '''
    Count the number of cply files in a directory
    '''
    cply_dir = bbdir + '/'+ moleculetype+ '/'
    files_all =  os.listdir( cply_dir )
    ftags_list = [ str(f[:-5]) for f in files_all if(f[-5:] == ".cply" ) ]
    return len(ftags_list),ftags_list

def write_ftag(filename,ftags_list):
    '''
    '''
    fout = open(filename,'wd')
    writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
    writer.writerow( ["ftag"] )
    for ftag in ftags_list:
         writer.writerow([ftag])
    fout.close()

def readcply(bbdir,moleculetype,ftag):
    '''
    Read building block from cply file
    '''    
    cply_file = bbdir + '/' + moleculetype+ '/'+ftag + '.cply'
    bb_i = buildingblock.Container()
    bb_i.read_cply(cply_file)
    return bb_i 

def load_json(bbdir,moleculetype,bb_i):
    '''
    Read building block from cply file
    '''
    json_file = bbdir + '/' + moleculetype+ '/calc_'+bb_i.tag + '_tddft_ob.json'
    try:
        with open(json_file) as f:            
            json_data = json.load(f)
            f.close()
            bb_i.properties.update(json_data['properties'])
    except IOError:
        logger.warning(" File not found %s in %s "%(json_file,os.getcwd()))
        print " File not found %s in %s "%(json_file,os.getcwd())

def countmols(schema,table,a_ctag):
    #
    # Check for the number of entries of a ctag in table 
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    #connection.autocommit = True
    cursor = connection.cursor()    
    tag_query = "SELECT COUNT(1) FROM   <schema>.<table> WHERE ctag = '<ctag>';"
    tag_query = tag_query.replace('<schema>',schema)
    tag_query = tag_query.replace('<table>',table)  
    # Find entries
    atag_query = copy.deepcopy(tag_query)
    sql = atag_query.replace('<ctag>',a_ctag)
    cursor.execute(sql)
    type_entries = cursor.fetchone()[0]
    logger.debug('%d entries foung for ctag %s '%(type_entries,a_ctag))
    # commit and close connection
    connection.commit()
    connection.close()
    # 
    return type_entries


def execute_sql(sql):
    '''
    Execute a SQL line with psycopg2
    '''
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()
    cursor.execute(sql)
    connection.commit()
    cursor.close()
    connection.close()


def insert_prop(schema,table,ctag,prop_key,prop_val,command):
    '''
    Update property 
    '''
    sql = "<command> <schema>.<table> SET <prop_key> = <prop_val> WHERE ctag = '<ctag>';"
    sql = sql.replace('<command>',command)
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',table)
    sql = sql.replace('<ctag>',ctag)
    sql = sql.replace('<prop_key>',prop_key)
    sql = sql.replace('<prop_val>',str(prop_val))
    execute_sql(sql)
    return
    #

def convert_belmol_strucC(schema,moleculetype,mol_id):
    '''
    Get data from database to create a streamm.buildingblock container object
    '''
    # Open connection
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()    
    # Set up base sql query 
    sql_root = "SELECT <prop_key> FROM  <schema>.<moleculetype> WHERE mol_id = <mol_id>;"
    sql_root = sql_root.replace('<schema>',schema)
    sql_root = sql_root.replace('<moleculetype>',moleculetype)
    sql_root = sql_root.replace('<mol_id>',str(mol_id))
    # print " sql_root ",sql_root
    # Get tag string
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','ctag')
    cursor.execute(sql)
    ctag = cursor.fetchone()[0]
    connection.commit()
    del sql
    # Get mol string
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','mol')
    cursor.execute(sql)
    molstr = cursor.fetchone()[0]
    connection.commit()
    del sql
    try:
        mol = pybel.readstring("sdf",molstr)
    except IOError:
        print "Load failed for %d" % (mol_id)
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','iupac')
    cursor.execute(sql)
    iupac = cursor.fetchone()[0]
    connection.commit()
    del sql
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','name')
    cursor.execute(sql)
    name = cursor.fetchone()[0]
    connection.commit()
    del sql
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','common_tag')
    cursor.execute(sql)
    common_tag = cursor.fetchone()[0]
    connection.commit()
    del sql
    # Create structure 
    strucC = buildingblock.Container(ctag)
    strucC.properties['common_tag'] = common_tag
    strucC.properties['ctag'] = ctag
    strucC.properties['name'] = name
    strucC.properties['iupac'] = iupac
    strucC.properties['deptag'] = ctag
    strucC.properties['moltype'] = moleculetype
    #
    #
    #
    for p_index in range(mol.OBMol.NumAtoms()):
        p_i = mol.OBMol.GetAtom(p_index+1)
        el_properties = periodictable.element_number(p_i.GetAtomicNum())
        BBatom_i = buildingblock.BBatom(el_properties['symbol'])
        # BBatom_i.properties.update(  )
        BBatom_i.tag = BBatom_i.properties['symbol']
        pos_i = [ p_i.GetX(), p_i.GetY(), p_i.GetZ() ]
        #BBatom_i.properties["mass"] = p_i.GetAtomicMass() 
        #BBatom_i.properties["mass"] = p_i.GetExactMass() 
        # BBatom_i.properties["charge"] = p_i.GetFormalCharge() 
        # BBatom_i.properties["label"] = p_i.GetAtomicMass() 
        # BBatom_i.properties["fftype"] = p_i.GetAtomicMass() 
        strucC.add_partpos(BBatom_i,pos_i, deepcopy = True)
    #
    # Get mol lists 
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','func_list')
    cursor.execute(sql)
    func_list = cursor.fetchone()[0]
    connection.commit()
    del sql
    for pkey in func_list:
        strucC.particles[pkey].properties['bbid'] = 'R'
    del func_list    
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','term_list')
    cursor.execute(sql)
    term_list = cursor.fetchone()[0]
    connection.commit()
    del sql
    for pkey in term_list:
        strucC.particles[pkey].properties['bbid'] = 'T'
    del term_list
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','sub_list')
    cursor.execute(sql)
    sub_list = cursor.fetchone()[0]
    connection.commit()
    del sql
    for pkey in sub_list:
        strucC.particles[pkey].properties['bbid'] = 'S'
    del sub_list
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','bond_list')
    cursor.execute(sql)
    bond_list = cursor.fetchone()[0]
    connection.commit()
    del sql
    #
    for bond in bond_list:
        pkey1 = bond[0]
        pkey2 = bond[1]
        bond_i = structure.Bond(pkey1,pkey2)
        strucC.add_bond(bond_i)
        # print ">bond_list ", strucC.n_bonds,pkey1,pkey2
    #
    strucC.bonded_nblist.build_nblist(strucC.particles,strucC.bonds )
    del bond_list
    # Close connection
    cursor.close()
    connection.close()
    # 
    set_bondorder = True
    if( set_bondorder ):
        for b_index in range(mol.OBMol.NumBonds()):
            b_i = mol.OBMol.GetBond(b_index)
            b1 = b_i.GetBeginAtomIdx() -1
            b2 = b_i.GetEndAtomIdx()-1
            bo = b_i.GetBondOrder()
            if( b1 > b2 ):
                t1 = copy.deepcopy(b1)
                b1 = b2
                b2 = t1 
            # print ">NumBonds ",b1,b2,bo
            for bkey,bond_i in strucC.bonds.iteritems():
                if( bond_i.pkey1 == b1 and  bond_i.pkey2 == b2 ):
                    bond_i.properties['bondorder'] = bo
                if( bond_i.pkey1 == b2 and  bond_i.pkey2 == b1 ):
                    bond_i.properties['bondorder'] = bo
        #for bkey,bond_i in strucC.bonds.iteritems():
        #    print ">set_bondorder ",bkey,bond_i.pkey1,bond_i.pkey2,bond_i.properties['bondorder'] 
    #
    return strucC
    
def convert_belmol_strucC2(cursor,schema,moleculetype,mol_id):
    '''
    Get data from database to create a streamm.buildingblock container object
    '''
    # Set up base sql query 
    sql_root = "SELECT <prop_key> FROM  <schema>.<moleculetype> WHERE mol_id = <mol_id>;"
    sql_root = sql_root.replace('<schema>',schema)
    sql_root = sql_root.replace('<moleculetype>',moleculetype)
    sql_root = sql_root.replace('<mol_id>',str(mol_id))
    # print " sql_root ",sql_root
    # Get tag string
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','ctag')
    cursor.execute(sql)
    ctag = cursor.fetchone()[0]
    del sql
    # Get mol string
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','mol')
    cursor.execute(sql)
    molstr = cursor.fetchone()[0]
    del sql
    try:
        mol = pybel.readstring("sdf",molstr)
    except IOError:
        print "Load failed for %d" % (mol_id)
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','iupac')
    cursor.execute(sql)
    iupac = cursor.fetchone()[0]
    del sql
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','name')
    cursor.execute(sql)
    name = cursor.fetchone()[0]
    del sql
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','common_tag')
    cursor.execute(sql)
    common_tag = cursor.fetchone()[0]
    del sql
    # Create structure 
    strucC = buildingblock.Container(ctag)
    strucC.properties['common_tag'] = common_tag
    strucC.properties['ctag'] = ctag
    strucC.properties['name'] = name
    strucC.properties['iupac'] = iupac
    strucC.properties['deptag'] = ctag
    strucC.properties['moltype'] = moleculetype
    #
    #
    #
    for p_index in range(mol.OBMol.NumAtoms()):
        p_i = mol.OBMol.GetAtom(p_index+1)
        el_properties = periodictable.element_number(p_i.GetAtomicNum())
        BBatom_i = buildingblock.BBatom(el_properties['symbol'])
        # BBatom_i.properties.update(  )
        BBatom_i.tag = BBatom_i.properties['symbol']
        pos_i = [ p_i.GetX(), p_i.GetY(), p_i.GetZ() ]
        #BBatom_i.properties["mass"] = p_i.GetAtomicMass() 
        #BBatom_i.properties["mass"] = p_i.GetExactMass() 
        # BBatom_i.properties["charge"] = p_i.GetFormalCharge() 
        # BBatom_i.properties["label"] = p_i.GetAtomicMass() 
        # BBatom_i.properties["fftype"] = p_i.GetAtomicMass() 
        strucC.add_partpos(BBatom_i,pos_i, deepcopy = True)
    #
    # Get mol lists 
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','func_list')
    cursor.execute(sql)
    func_list = cursor.fetchone()[0]
    del sql
    for pkey in func_list:
        strucC.particles[pkey].properties['bbid'] = 'R'
    del func_list    
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','term_list')
    cursor.execute(sql)
    term_list = cursor.fetchone()[0]
    del sql
    for pkey in term_list:
        strucC.particles[pkey].properties['bbid'] = 'T'
    del term_list
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','sub_list')
    cursor.execute(sql)
    sub_list = cursor.fetchone()[0]
    del sql
    for pkey in sub_list:
        strucC.particles[pkey].properties['bbid'] = 'S'
    del sub_list
    #
    sql = copy.deepcopy(sql_root)
    sql = sql.replace('<prop_key>','bond_list')
    cursor.execute(sql)
    bond_list = cursor.fetchone()[0]
    del sql
    #
    for bond in bond_list:
        pkey1 = bond[0]
        pkey2 = bond[1]
        bond_i = structure.Bond(pkey1,pkey2)
        strucC.add_bond(bond_i)
        # print ">bond_list ", strucC.n_bonds,pkey1,pkey2
    #
    strucC.bonded_nblist.build_nblist(strucC.particles,strucC.bonds )
    del bond_list
    # 
    set_bondorder = True
    if( set_bondorder ):
        for b_index in range(mol.OBMol.NumBonds()):
            b_i = mol.OBMol.GetBond(b_index)
            b1 = b_i.GetBeginAtomIdx() -1
            b2 = b_i.GetEndAtomIdx()-1
            bo = b_i.GetBondOrder()
            if( b1 > b2 ):
                t1 = copy.deepcopy(b1)
                b1 = b2
                b2 = t1 
            # print ">NumBonds ",b1,b2,bo
            for bkey,bond_i in strucC.bonds.iteritems():
                if( bond_i.pkey1 == b1 and  bond_i.pkey2 == b2 ):
                    bond_i.properties['bondorder'] = bo
                if( bond_i.pkey1 == b2 and  bond_i.pkey2 == b1 ):
                    bond_i.properties['bondorder'] = bo
        #for bkey,bond_i in strucC.bonds.iteritems():
        #    print ">set_bondorder ",bkey,bond_i.pkey1,bond_i.pkey2,bond_i.properties['bondorder'] 
    #
    return strucC
    

'''
for pid, ptclObj  in strucC.particles.iteritems():
   print pid,ptclObj.properties['symbol'],ptclObj.properties['bbid']

print 'mol modstyle 1 0 VDW 0.200000 12.000000'
print 'mol modstyle 2 0 VDW 0.400000 12.000000'

for pid, ptclObj  in strucC.particles.iteritems():
    sel_j = 'mol modselect 1 0 index %s \n'%(pid)
    sel_j += 'mol modselect 2 0 index  '
    nb_list_i = strucC.bonded_nblist.getnbs(pid)
    for p_j in  nb_list_i[0:-1]:
        sel_j += "%d or index "%(p_j)
    
    sel_j += "%d  "%(strucC.bonded_nblist.getnbs(pid)[-1])
    print sel_j

    
   print pid,ptclObj.properties['symbol'],ptclObj.properties['bbid']

        
for b_index in range(mol.OBMol.NumBonds()):
    b_i = mol.OBMol.GetBond(b_index)
    # print b_i.GetBeginAtomIdx() ,b_i.GetEndAtomIdx(),b_i.GetBondOrder ()
    pkey1 = b_i.GetBeginAtomIdx() -1
    pkey2 = b_i.GetEndAtomIdx() -1 
    bond_i = structure.Bond(pkey1,pkey2)
    bond_i.border  = b_i.GetBondOrder()
    strucC.add_bond(bond_i)
    
# Check bonds 
for k,b in strucC.bonds.iteritems():
    if( [b.pkey1,b.pkey2] not in bond_list and [b.pkey2,b.pkey1] not in bond_list ):
        print "Bad bond found in mol file %d - %d  "%( b.pkey1,b.pkey2)
'''

def ob_properties(pybelmol):
    '''
    Find openbabel properties of a molecule and return them as a dictionary
    '''
    properties = dict()
    #
    desc = pybelmol.calcdesc(['atoms','bonds','HBA1','HBA2','MW','MR','logP','TPSA','HBD','sbonds','abonds','tbonds','dbonds'])
    for pkey,prop_val in  desc.iteritems():
        properties[pkey] = prop_val
    # 
    #properties['natoms'] = properties['atoms']
    #del properties['atoms']
    #properties['nbonds'] = properties['bonds']
    #del properties['bonds']
    
    mass = pybelmol.exactmass
    properties['mass'] = mass
    #
    nrings = 0
    rings = []
    for ring in openbabel.OBMolRingIter(pybelmol.OBMol):
        nrings = nrings + 1
        rings.append(ring.Size())
    properties['nrings'] = nrings
    ring_size = ",".join(map(str,rings))
    properties['rings'] = ring_size
    nrbonds = pybelmol.OBMol.NumRotors()
    properties['nrbonds'] = nrbonds
    #
    fpbits = pybelmol.calcfp().bits
    if( len(fpbits) > 0 ):
        fplong = reduce(lambda x,y: x|y,map(lambda x: 1<<(x-1),fpbits))
    else:
        fpbits = 001
        fplong = fpbits
    properties['fplong'] = fplong
    # 
    fp4bits = pybelmol.calcfp(fptype='FP4').bits
    if( len(fp4bits) > 0 ):
        fp4long = reduce(lambda x,y: x|y,map(lambda x: 1<<(x-1),fp4bits))
    else:
        fp4bits = 001
        fp4long = fp4bits
    properties['fp4long'] = fp4long
    #
    return properties
    
   
def set_prop_names_desc():
    '''
    '''
    prop_name = dict()
    prop_name['TPSA'] = 'tpsa'
    prop_name['dbonds'] = 'dbonds'
    prop_name['HBA2'] = 'hba2'
    prop_name['fpbits'] = ''
    prop_name['logP'] = 'logp'
    prop_name['bonds'] = 'nbonds'
    prop_name['HBD'] = 'hbd'
    prop_name['rings'] = ''
    prop_name['nrbonds'] = 'nrbonds'
    prop_name['atoms'] = 'natoms'
    prop_name['tbonds'] = 'tbonds'
    prop_name['MW'] = 'weight'
    prop_name['mass'] = 'mass'
    prop_name['nrings'] = 'nrings'
    prop_name['fplong'] = 'fingerprint'
    prop_name['fp4long'] = 'fingerprint-fp4'
    prop_name['MR'] = 'mr'
    prop_name['HBA1'] = 'hba1'
    prop_name['abonds'] = 'abonds'
    prop_name['sbonds'] = 'sbonds'
    #
    prop_desc = dict()
    prop_desc['TPSA'] = 'OpenBabel Total Polar Surface Area'
    prop_desc['dbonds'] = 'OpenBabel Double Bond Count'
    prop_desc['HBA2'] = 'OpenBabel Hydrogen Bonding Ability Parameter 2'
    prop_desc['fpbits'] = ''
    prop_desc['logP'] = 'OpenBabel Octinol-Water Partition Coefficent'
    prop_desc['bonds'] = 'OpenBabel Num Bonds'
    prop_desc['HBD'] = 'OpenBabel Hydrogen Bond Donors'
    prop_desc['rings'] = ''
    prop_desc['nrbonds'] = 'OpenBabel Rotatable Bond Count'
    prop_desc['atoms'] = 'OpenBabel Num Atoms'
    prop_desc['tbonds'] = 'OpenBabel Triple Bond Count'
    prop_desc['MW'] = 'OpenBabel Molecular Weight'
    prop_desc['mass'] = 'OpenBabel Exact Mass'
    prop_desc['nrings'] = 'Number of Rings'
    prop_desc['fplong'] = 'OpenBabel FP2/Daylight'
    prop_desc['fp4long'] = 'OpenBabel FP4'
    prop_desc['MR'] = 'OpenBabel Molecular Refectivity'
    prop_desc['HBA1'] = 'OpenBabel Hydrogen Bonding Ability Parameter 1'
    prop_desc['abonds'] = 'OpenBabel Aromatic Bond Count'
    prop_desc['sbonds'] = 'OpenBabel Single Bond Count'
    return prop_name,prop_desc
      
def count_bbid(bblockC_i):
    #
    bblockC_i.term_cnt = 0
    bblockC_i.func_cnt = 0
    bblockC_i.sub_cnt = 0
    for pkey_i, particle_i  in bblockC_i.particles.iteritems():
        if( particle_i.properties["bbid"] == "R" ):
             bblockC_i.func_cnt += 1 
        elif( particle_i.properties["bbid"] == "T" ):
             bblockC_i.term_cnt += 1 
        elif( particle_i.properties["bbid"] == "S" ):
             bblockC_i.sub_cnt += 1

def convert_strucC_openbabel_direct(strucC):
    '''
    Convert structure container to open babel mol
    '''
    
    mol = openbabel.OBMol()
    for pkey_i, particle_i  in strucC.particles.iteritems():

        pos_i = strucC.positions[pkey_i]
        a = mol.NewAtom()
        a.SetAtomicNum(particle_i.properties["number"])   # carbon atom
        a.SetVector(pos_i[0],pos_i[1],pos_i[2]) # coordinates

    # This is wrong b_n needs to be the bond order 
    if( strucC.n_bonds > 0 ):
        for bkey_i, bond_i  in strucC.bonds.iteritems():
            b_n = bkey_i + 1
            # 1,2,3 - (single,double,triple) bond
            b_i = bond_i.pkey1 + 1 
            b_j = bond_i.pkey2 + 1
            if( b_i < b_j ):
                mol.AddBond(b_i, b_j, b_n)   # atoms indexed from 1
            
    return mol

             
def pybel_opt(strucC):
    '''
    Convert structure container to open babel mol
    '''
    #mol = convert_strucC_openbabel(strucC)
    pybelmol = convert_strucCxyz_pybel(strucC)
    pybelmol.localopt()
    
    for pkey_i, particle_i  in strucC.particles.iteritems():
        a = pybelmol.atoms[pkey_i]
        strucC.positions[pkey_i][0] = a.vector.GetX()
        strucC.positions[pkey_i][1] = a.vector.GetY()
        strucC.positions[pkey_i][2] = a.vector.GetZ()
    
    
    return pybelmol

                                  
def functionalize_fullsym(bblockC_i,func_list,blank_str="H",blank_deptag = "R0",verbose=False,obopt=True):
    '''
    Functionalize building block
    using configs="fullsym"
    
    This allows for only  B_(A,A), B_(C,C), B_(A,C) ...
    but not redundant 'symetric' configurations  B_(C,A)
    type configurations

    Arguments:
            bblockC_i (Container) Buildingblock container 1 
            func_list (list of Containers) functional groups to attached to bblockC_i
            blank_str (str) name of non subtituted group 
    Returns:
            config_list (list of Containers) of functionalized bblockC_i
            
    '''    
    
    #
    debug = False   
    if( debug ):
        print "!!!!!!!!!!!!!!!!!!!!!!!  functionalize_fullsym debug !!!!!!!!!!!!!!!!!!!!!!!"
    count_bbid(bblockC_i)
    n_pos =  bblockC_i.func_cnt
    #
    logging.info(' Functionalizing  %s with %d functionalizable positions '%(bblockC_i.tag,n_pos))
    logging.info('                          %d term positions '%(bblockC_i.term_cnt))
    logging.info('                          %d sub positions '%(bblockC_i.sub_cnt))
    f_overlap = open('overlap_%s.text'%(bblockC_i.tag),'w')
    f_overlap.write('# Initial tag, Being add to tag, being added tag, position , final tag not produced  \n')
    # Add 
    config_list = []
    bbC_i = copy.deepcopy(bblockC_i)
    bbC_i.properties['common_tag'] += "_("
    bbC_i.tag += "_"
    bbC_i.properties['ctag'] = bbC_i.tag
    bbC_i.properties['deptag'] += "_"
    bbC_i.properties['name']  += "_"
    config_list.append(bbC_i)
    #
    func_pos_list = []
    func_pos_list.append(-1)
    #
    for pos_i in range(n_pos):
        config_list_i = []
        func_pos_list_i = []
        #
        if( debug ):
            config_cnt = 0
            print "pos_i ",pos_i
            print " func_pos_list ",func_pos_list
        #
        for bb_indx in range(len(config_list)):
            bblockconfig_i = config_list[bb_indx]
            last_func = func_pos_list[bb_indx]
            if( debug ):
                print "  bbC_i.tag",bblockconfig_i.tag
                print "  last_func",last_func
            #
            for funcx_i in range(len(func_list)):
                if( funcx_i >= last_func ):
                    bb_i = copy.deepcopy(bblockconfig_i)
                    func_i = func_list[funcx_i] 
                    bb_j = copy.deepcopy(func_i)
                    #   
                    # bb_k = buildingblock.attach(bb_k,func_j,"R",0,"R",0)
                    #
                    bb_i = bb_i.prepattach("R",0,0,-1,0.0,debug = debug)
                    # bb_j = func_i.prepattach("R",0,0,1,angle_rad,debug = debug)
                    #
                    bbC_i,bbC_j =  buildingblock.shiftprep(bb_i,bb_j,debug = debug )
                    bbC_i =  buildingblock.attachprep(bbC_i,bbC_j,debug = debug )
                    #
                    if( obopt ):
                        # pybelmol = ob.pybel_opt(bbC_i)
                        pybelmol =convert_struc_pymol(bbC_i)
                        pybelmol.localopt()
                        for pkey_i, particle_i  in bbC_i.particles.iteritems():
                            a = pybelmol.atoms[pkey_i]
                            bbC_i.positions[pkey_i][0] = a.vector.GetX()
                            bbC_i.positions[pkey_i][1] = a.vector.GetY()
                            bbC_i.positions[pkey_i][2] = a.vector.GetZ()
                    #
                    if( pos_i > 0 ):
                        bb_j.properties['common_tag'] = ","+bb_j.properties['common_tag']
                    #
                    bbC_i.tag += bb_j.tag
                    bbC_i.properties['ctag'] +=  bb_j.properties['ctag']
                    bbC_i.properties['deptag'] += bb_j.properties['deptag']
                    bbC_i.properties['name'] += bb_j.properties['name']
                    bbC_i.properties['common_tag'] += bb_j.properties['common_tag']
                    bbC_i.properties['iupac'] += bb_j.properties['iupac']
                    #
                    if( debug ):
                        config_cnt += 1
                        print " join count",config_cnt,bbC_i.tag,bbC_i.properties['ctag']
                    # bb_k.tag += ")"
                    config_list_i.append(copy.deepcopy(bbC_i))
                    func_pos_list_i.append(funcx_i)
                    logger.info(" Created config %s %s "%(bbC_i.properties['ctag'],bbC_i.properties['common_tag']))
        #
        config_list = config_list_i
        func_pos_list = func_pos_list_i
    # Add ) cap to end of tag 
    for bb_indx in range(len(config_list)):
        bb_i = config_list[bb_indx]
        # Set substitution positions to blanks
        for pos_blank in range(bb_i.sub_cnt):
            bb_i.properties['common_tag'] += ','+blank_str
            bb_i.properties['ctag'] += blank_deptag
        bb_i.properties['common_tag'] += ")"    
    #
    f_overlap.close()
    #
    if( debug ):
        for bb_indx in range(len(config_list)):
            bb_i = config_list[bb_indx]
            print bb_i.tag
        print "!!!!!!!!!!!!!!!!!!!!!!!  functionalize_fullsym debug !!!!!!!!!!!!!!!!!!!!!!!"
        # sys.exit("ifnuw8tb27bm48")
    #
    return config_list

def convert_pos_bohr_ang(struc_o):
    for pkey_o  in struc_o.particles.keys():
        pos_o = struc_o.positions[pkey_o]       
        pos_i = [ units.convert_bohr_ang(v_i) for v_i in pos_o ] 
        struc_o.positions[pkey_o] = pos_i
        
def proc_gaus():

    prop_name,prop_desc = set_prop_names_desc()

    bbdir = "/home/tkemper/BuildingBlocks-cand"
    bbtype = 'functional_groups'
    moleculetype = bbtype
    tag_o = 'OAcDA'

    peregrine = resource.Resource('peregrine')
    peregrine.load_json()

    # Copy cply to scratch dir
    gaus_i = gaussian.Gaussian("%s_tddft"%tag_o)
    gaus_i.set_resource(peregrine)
    gaus_i.dir['BuildingBlocks'] = "%s/%s"%(bbdir,bbtype)
    gaus_i.dir['scratch'] = '/scratch/%s/%s'%(peregrine.ssh['username'],gaus_i.tag)
    print gaus_i.dir['scratch']
    print gaus_i.dir['storage']

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
    gaus_i.files['output']['cply'] = "%s.cply"%(gaus_i.strucC.tag)
    gaus_i.files['output']['xyz'] = "%s.xyz"%(gaus_i.strucC.tag)

    bb_i = gaus_i.strucC

    schema = 'opv_cand'
    #strprop_dic = {'iupac':'IUPAC','name':'name','common_tag':'common_tag','mol':'molstr'}
    strprop_dic = {'name':'name','common_tag':'common_tag','mol':'molstr'}
    intprop_dic = {'mol_id':'mol_id'}
    listprop_dic = {'term_list':'term_list','func_list':'func_list','sub_list':'sub_list','bond_list':'bond_list'}

    bb_i.proc_bbid()
    bb_i.properties['mol_id']  = mol_id
    bb_i.properties['bond_list'] = []
    for bkey,b_i in bb_i.bonds.iteritems():
        bb_i.properties['bond_list'].append([b_i.pkey1,b_i.pkey2])

    pybelmol = convert_struc_pymol(bb_i)

    print bb_i.properties['molstr'] 

    #if( True ):
    type_entries = countmols(schema,moleculetype,bb_i.properties['ctag'])
    print " type_entries %d "%(type_entries)

    if( type_entries == 0 ):
        # Create entry 
        connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
        cursor = connection.cursor()       
        sql = "INSERT INTO <schema>.<table> (mol_id,ctag,iupac,name,common_tag,mol) VALUES (%s,%s,%s,%s,%s,%s)"
        sql = sql.replace('<schema>',schema)
        sql = sql.replace('<table>',moleculetype)
        print sql,mol_id,bb_i.properties['ctag'],bb_i.properties['IUPAC'],bb_i.properties['name'],bb_i.properties['common_tag'],bb_i.properties['molstr'] 
        cursor.execute(sql,(mol_id,bb_i.properties['ctag'],bb_i.properties['IUPAC'],bb_i.properties['name'],bb_i.properties['common_tag'],bb_i.properties['molstr'] ))
        connection.commit()
        cursor.close()
        connection.close()
    else:
        command = 'UPDATE'
        for prop_key,pkey in strprop_dic.iteritems():
            prop_val = '\'%s\''%bb_i.properties[pkey]
            insert_prop(schema,moleculetype,bb_i.properties['ctag'],prop_key,prop_val,command)

    command = 'UPDATE'
    for prop_key,pkey in listprop_dic.iteritems():
        prop_val = '\'%s\''%bb_i.properties[pkey]
        prop_val = prop_val.replace(']','}')
        prop_val = prop_val.replace('[','{')
        print  prop_key,pkey,prop_val
        insert_prop(schema,moleculetype,bb_i.properties['ctag'],prop_key,prop_val,command)

    table = "%s_calc_prop"%(moleculetype)
    properties = ob_properties(pybelmol)
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()
    cursor.execute("BEGIN;")
    for prop_key,prop_val in properties.iteritems():
        name_i = prop_name[prop_key]
        desc_i = prop_desc[prop_key]
        if( len(name_i) and len(desc_i) ):
            print prop_key,name_i,desc_i,prop_val
            sql = "INSERT INTO <schema>.<table> (mol_id,property_name,property_method,property_value) VALUES (%s,%s,%s,%s)"
            sql = sql.replace('<schema>',schema)
            sql = sql.replace('<table>',table)
            print sql,mol_id,name_i,desc_i,prop_val
            cursor.execute(sql,(mol_id,name_i,desc_i,prop_val))
    cursor.execute("COMMIT")
    connection.commit()
    cursor.close()
    connection.close()

                
def add_struc(calc_tag,options):

    local = resource.Resource('local')
    local.load_json()
    peregrine = resource.Resource('peregrine')
    peregrine.load_json()
    #
    #script_dir = '/Users/tkemper/Projects/opv-project/scripts'
    #os.chdir(script_dir)
    #
    bbdir ="/home/tkemper/BuildingBlocks-cand_opt/" # options.bbdir  #
    #bbdir = '/Users/tkemper/Projects/opv-project/BuildingBlocks-cand_opt'
    schema = 'opv_cand'
    #strprop_dic = {'iupac':'IUPAC','name':'name','common_tag':'common_tag','mol':'molstr'}
    strprop_dic = {'name':'name','common_tag':'common_tag','mol':'molstr'}
    intprop_dic = {'mol_id':'mol_id'}
    listprop_dic = {'term_list':'term_list','func_list':'func_list','sub_list':'sub_list','bond_list':'bond_list'}
    os.chdir(bbdir)
    moleculetype = 'acceptors'
    ftag_id = 0 
    #
    moleculetypes = ['donors','functional_groups','acceptors'] # 'spacers']
    for moleculetype in moleculetypes:
        print "moleculetype %s "%(moleculetype)
        #
        # if( True ):
        n_ftags,ftags_list = countcplys(bbdir,moleculetype)
        output = 'ftags_%s.csv'%(moleculetype)
        write_ftag(output,ftags_list)
        #
        for ftag_id in range(len(ftags_list)):
            #if( True ):
            ftag = ftags_list[ftag_id]
            print '%s %s '%(ftag_id,ftag)
            # ftag_id = 0 
            ftag = ftags_list[ftag_id]

            bb_i = buildingblock.Container(ftag)
            cply_file = bbdir + '/'+ moleculetype + '/'+ftag + '.cply'
            bb_i.read_cply(cply_file)
            bb_i.proc_bbid()
            load_json(bbdir,moleculetype,bb_i)
            bb_i.properties['mol_id']  = ftag_id
            bb_i.properties['bond_list'] = []
            for bkey,b_i in bb_i.bonds.iteritems():
                bb_i.properties['bond_list'].append([b_i.pkey1,b_i.pkey2])
            #if( True ):
            type_entries = countmols(schema,moleculetype,bb_i.properties['ctag'])
            print " type_entries %d "%(type_entries)
            if( type_entries > 0 ):
                command = 'UPDATE'
            else:
                command = 'INSERT'
            #if( True ):
            for prop_key,pkey in strprop_dic.iteritems():
                prop_val = '\'%s\''%bb_i.properties[pkey]
                insert_prop(schema,moleculetype,bb_i.properties['ctag'],prop_key,prop_val,command)
            #if( True ):
            for prop_key,pkey in intprop_dic.iteritems():
                prop_val = bb_i.properties[pkey]
                print  prop_key,pkey,prop_val
                insert_prop(schema,moleculetype,bb_i.properties['ctag'],prop_key,prop_val,command)
            #if( True ):
            for prop_key,pkey in listprop_dic.iteritems():
                prop_val = '\'%s\''%bb_i.properties[pkey]
                prop_val = prop_val.replace(']','}')
                prop_val = prop_val.replace('[','{')
                print  prop_key,pkey,prop_val
                insert_prop(schema,moleculetype,bb_i.properties['ctag'],prop_key,prop_val,command)

             

def cal_obprops(calc_tag,options):
    '''
     for prop_key,prop_val in properties.iteritems():
        print "prop_desc[\'%s\'] = \'\'"%(prop_key)
    for prop_key,prop_val in properties.iteritems():
        print "prop_name[\'%s\'] = \'\'"%(prop_key)
    '''

    local = resource.Resource('local')
    local.load_json()
    peregrine = resource.Resource('peregrine')
    peregrine.load_json()

    #script_dir = '/Users/tkemper/Projects/opv-project/scripts'
    #os.chdir(script_dir)

    schema = 'opv_cand'
    
    prop_name,prop_desc = set_prop_names_desc()
    
    moleculetypes = ['functional_groups','donors','acceptors'] # 'spacers']
    for moleculetype in moleculetypes:
        print "moleculetype %s "%(moleculetype)
        table = "%s_calc_prop"%(moleculetype)

        sql = "ALTER SEQUENCE <moleculetype>_calc_prop_id_seq RESTART WITH 1;"
        sql = sql.replace('<moleculetype>',moleculetype)
        print sql
        sql = "DELETE FROM <moleculetype>_calc_prop ;"
        sql = sql.replace('<moleculetype>',moleculetype)
        print sql
        '''
        connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
        cursor = connection.cursor()
        cursor.execute(sql)
        connection.commit()
        cursor.close()
        connection.close()        
        '''
    for moleculetype in moleculetypes:
        print "moleculetype %s "%(moleculetype)
        table = "%s_calc_prop"%(moleculetype)

        connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
        cursor = connection.cursor()
        sql = "SELECT mol_id FROM  <schema>.<table>;"
        sql = sql.replace('<schema>',schema)
        sql = sql.replace('<table>',moleculetype)
        cursor.execute(sql)
        mol_ids = cursor.fetchall()
        connection.commit()
        cursor.close()
        connection.close()
        print "mol_ids %s "%(str(mol_ids))

        for mol_return in mol_ids:
            (mol_id,) = mol_return
            print "mol_id %d "% mol_id

            connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
            cursor = connection.cursor()
            sql = "SELECT mol FROM  <schema>.<table> WHERE mol_id = <mol_id>;"
            sql = sql.replace('<schema>',schema)
            sql = sql.replace('<table>',moleculetype)
            sql = sql.replace('<mol_id>',str(mol_id))
            cursor.execute(sql)
            molstr = cursor.fetchone()[0]
            connection.commit()
            sql = "SELECT ctag FROM  <schema>.<table> WHERE mol_id = <mol_id>;"
            sql = sql.replace('<schema>',schema)
            sql = sql.replace('<table>',moleculetype)
            sql = sql.replace('<mol_id>',str(mol_id))
            cursor.execute(sql)
            ctag = cursor.fetchone()[0]
            connection.commit()
            cursor.close()
            connection.close()
            print "ctag %s "% ctag

            add_props = False 
            try:
                pybelmol = pybel.readstring("sdf",molstr)
                add_props = True 
            except IOError:
                print "Load failed for %d" % (mol_id)
            if( add_props ):
                logger.info("Adding properties for mol_id %d ctag %s "%(mol_id,ctag))
                #
                properties = ob_properties(pybelmol)
                #
                connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                cursor = connection.cursor()
                cursor.execute("BEGIN;")
                for prop_key,prop_val in properties.iteritems():
                    name_i = prop_name[prop_key]
                    desc_i = prop_desc[prop_key]
                    if( len(name_i) and len(desc_i) ):
                        print prop_key,name_i,desc_i,prop_val
                        sql = "INSERT INTO <schema>.<table> (mol_id,property_name,property_method,property_value) VALUES (%s,%s,%s,%s)"
                        sql = sql.replace('<schema>',schema)
                        sql = sql.replace('<table>',table)
                        print sql,mol_id,name_i,desc_i,prop_val
                        cursor.execute(sql,(mol_id,name_i,desc_i,prop_val))
                #
                cursor.execute("COMMIT")
                connection.commit()
                cursor.close()
                connection.close()
          
def func_all(calc_tag,options,p):
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()

    local = resource.Resource('local')
    local.load_json()
    peregrine = resource.Resource('peregrine')
    peregrine.load_json()
    bbdir ="/home/tkemper/BuildingBlocks-cand_opt" # options.bbdir  #
    bbdir = '/Users/tkemper/Projects/opv-project/BuildingBlocks-cand_opt'
    #script_dir = '/Users/tkemper/Projects/opv-project/scripts'
    #os.chdir(script_dir)
    prop_name,prop_desc = set_prop_names_desc()


    strprop_dic = {'name':'name','common_tag':'common_tag','mol':'molstr'}
    listprop_dic = {'term_list':'term_list','func_list':'func_list','sub_list':'sub_list','bond_list':'bond_list'}
    
    schema = 'opv_cand'
    funcmoltype = 'functional_groups'
    configs = 'fullsym'
    # Prep functional groups    
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()
    sql = "SELECT mol_id FROM  <schema>.<table>;"
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',funcmoltype)
    cursor.execute(sql)
    mol_ids = cursor.fetchall()
    connection.commit()
    cursor.close()
    connection.close()
    print "mol_ids %s "%(str(mol_ids))    
    angle_rad = 90.0*math.pi/180.0

    dir_i = bbdir +"/"+funcmoltype
    if ( not os.path.isdir(dir_i) ):
        print "Making %s "%(dir_i)
        os.mkdir(dir_i)
    os.chdir(dir_i)


    fg_list_prep = []
    #print " !!!!!!!!!!!!!!!! hack  1 !!!!!!!!!!!!!!!!"
    #mol_ids = [(6,)]
    
    for mol_return in mol_ids:
        (mol_id,) = mol_return
        bb_i = convert_belmol_strucC(schema,funcmoltype,mol_id)
        bb_i.write_xyz()
        logger.info("Prepattch mol_id %d ctag %s common_tag %s "%(mol_id,bb_i.properties['ctag'],bb_i.properties['common_tag']))
        bb_i_prep = bb_i.prepattach(bbid_i="R",n_i=0,Xn_i=0,dir=1,yangle=angle_rad)
        bb_i_prep.write_xyz("%s_prep.xyz"%(bb_i_prep.tag))
        fg_list_prep.append(bb_i_prep)
        
    logger.info("All %s read in "%(funcmoltype))
    # sys.exit()
    
    moleculetypes = ['donors','acceptors'] # 'spacers']    
    for moleculetype in moleculetypes:
        print "moleculetype %s "%(moleculetype)
        func_moltype = "%s_func_%s"%(moleculetype,configs)
        func_mol_id = 0

        core_dir = bbdir +"/"+moleculetype
        if ( not os.path.isdir(core_dir) ):
            print "Making %s "%(core_dir)
            os.mkdir(core_dir)
        os.chdir(core_dir)
        
        core_func_dir = bbdir +"/"+func_moltype
        if ( not os.path.isdir(core_func_dir) ):
            print "Making %s "%(core_func_dir)
            os.mkdir(core_func_dir)
        os.chdir(core_func_dir)
        #
        connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
        cursor = connection.cursor()
        sql = "SELECT mol_id FROM  <schema>.<table>;"
        sql = sql.replace('<schema>',schema)
        sql = sql.replace('<table>',moleculetype)
        cursor.execute(sql)
        mol_ids = cursor.fetchall()
        connection.commit()
        cursor.close()
        connection.close()

        #print " !!!!!!!!!!!!!!!! hack  2 !!!!!!!!!!!!!!!!"
        #mol_ids = [(3,)]
        print "mol_ids %s "%(str(mol_ids))
        #
        for mol_return in mol_ids:
            (mol_id,) = mol_return
            print "mol_id %d "% mol_id 
            core_i = convert_belmol_strucC(schema,moleculetype,mol_id)
            os.chdir(core_dir)
            core_i.write_xyz()

            # Functionalize core 
            #functionalized_list = functionalize(core_i,fg_list_prep,configs=configs,verbose=False,obopt=self.obopt )

            os.chdir(core_func_dir)
            
            functionalized_list = functionalize_fullsym(core_i,fg_list_prep)
            logger.info(" Building block %s to %s table in %s  functionalized with:"%(core_i.tag,schema,func_moltype))
            for bb_i in functionalized_list:
                logger.info("   - %s %s "%(bb_i.properties['ctag'],bb_i.properties['common_tag']))
            for bb_i in functionalized_list:
                logger.info("Adding %s %s "%(bb_i.properties['ctag'],bb_i.properties['common_tag']))
                bb_i.write_cply()
                bb_i.write_xyz()    
                #
                bb_i.proc_bbid()
                bb_i.properties['bond_list'] = []
                for bkey,b_i in bb_i.bonds.iteritems():
                    bb_i.properties['bond_list'].append([b_i.pkey1,b_i.pkey2])
                # Get mol str     
                xyz_str = bb_i.write_xyz_str()
                pybelmol = pybel.readstring('xyz',xyz_str)
                #
                molstr  = pybelmol.write('mol')
                bb_i.properties['molstr'] = molstr
                #
                properties = ob_properties(pybelmol)
                #
                # Check properties from openbabel
                #
                if( bb_i.n_particles != properties['atoms'] ):
                        print "%s has an Error in atoms "%(bb_i.properties['ctag'])
                        sys.exit(1)
                if( bb_i.n_bonds != properties['bonds'] ):
                        print "%s has an Error in bonds "%(bb_i.properties['ctag'])
                        print bb_i.n_bonds , properties['bonds']
                        sys.exit(2)
                # 
                type_entries = countmols(schema,func_moltype,bb_i.properties['ctag'])
                if( type_entries == 0 ):
                    # Create entry 
                    logger.info("Creating new entry in %s.%s "%(schema,func_moltype))
                    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                    cursor = connection.cursor()       
                    sql = "INSERT INTO <schema>.<table> (ctag,iupac,name,common_tag,mol) VALUES (%s,%s,%s,%s,%s)"
                    sql = sql.replace('<schema>',schema)
                    sql = sql.replace('<table>',func_moltype)
                    cursor.execute(sql,(bb_i.properties['ctag'],bb_i.properties['iupac'],bb_i.properties['name'],bb_i.properties['common_tag'],bb_i.properties['molstr']))
                    connection.commit()
                    cursor.close()
                    connection.close()
                else:
                    logger.info("Entry in %s.%s found updating properties "%(schema,func_moltype))
                    command = 'UPDATE'

                    sql_root = "<command> <schema>.<table> SET <prop_key> = <prop_val> WHERE ctag = '<ctag>';"
                    sql_root = sql_root.replace('<command>',command)
                    sql_root = sql_root.replace('<schema>',schema)
                    sql_root = sql_root.replace('<table>',func_moltype)
                    sql_root = sql_root.replace('<ctag>',bb_i.properties['ctag'])
                    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                    cursor = connection.cursor()                    
                    for prop_key,pkey in strprop_dic.iteritems():
                        prop_val = '\'%s\''%bb_i.properties[pkey]
                        logger.info("  Updating %s "%(prop_key))
                        sql = copy.deepcopy(sql_root)
                        sql = sql.replace('<prop_key>',prop_key)
                        sql = sql.replace('<prop_val>',str(prop_val))
                        # insert_prop(schema,func_moltype,bb_i.properties['ctag'],prop_key,prop_val,command)
                        cursor.execute(sql)
                    connection.commit()
                    cursor.close()
                    connection.close()

                add_props = False
                #if( add_props ):
                logger.info("Checking for  %s in %s.%s "%(bb_i.properties['ctag'],schema,func_moltype))
                
                connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                cursor = connection.cursor()       
                sql = "SELECT COUNT(1) FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                sql = sql.replace('<schema>',schema)
                sql = sql.replace('<table>',func_moltype)
                sql = sql.replace('<ctag>',bb_i.properties['ctag'])
                cursor.execute(sql)
                mol_id_cnt = cursor.fetchone()[0]
                connection.commit()
                cursor.close()
                connection.close()
                if( mol_id_cnt == 0 ):
                    logger.warning(" ctag %s not found in %s.%s "%(bb_i.properties['ctag'],schema,func_moltype))
                else:
                    # Get mol_id
                    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                    cursor = connection.cursor()       
                    sql = "SELECT mol_id FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                    sql = sql.replace('<schema>',schema)
                    sql = sql.replace('<table>',func_moltype)
                    sql = sql.replace('<ctag>',bb_i.properties['ctag'])
                    cursor.execute(sql)
                    mol_id = cursor.fetchone()[0]
                    connection.commit()
                    cursor.close()
                    connection.close()
                    # Add lists 
                    command = 'UPDATE'
                    sql_root = "<command> <schema>.<table> SET <prop_key> = <prop_val> WHERE ctag = '<ctag>';"
                    sql_root = sql_root.replace('<command>',command)
                    sql_root = sql_root.replace('<schema>',schema)
                    sql_root = sql_root.replace('<table>',func_moltype)
                    sql_root = sql_root.replace('<ctag>',bb_i.properties['ctag'])
                    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                    cursor = connection.cursor()                    
                    for prop_key,pkey in listprop_dic.iteritems():
                        prop_val = '\'%s\''%bb_i.properties[pkey]
                        prop_val = prop_val.replace(']','}')
                        prop_val = prop_val.replace('[','{')
                        logger.info("  Updating %s "%(prop_key))
                        sql = copy.deepcopy(sql_root)
                        sql = sql.replace('<prop_key>',prop_key)
                        sql = sql.replace('<prop_val>',str(prop_val))
                        # insert_prop(schema,func_moltype,bb_i.properties['ctag'],prop_key,prop_val,command)
                        cursor.execute(sql)
                    connection.commit()
                    cursor.close()
                    connection.close()
                    # Add properties 
                    prop_table = "%s_calc_prop"%(func_moltype)
                    logger.info("Adding properties for mol_id %d ctag %s to %s"%(mol_id,bb_i.properties['ctag'],prop_table))
                    #
                    command = 'UPDATE'
                    sql_root = "<command> <schema>.<table> SET property_value = <prop_val> WHERE mol_id = <mol_id> AND property_name = \'<property_name>\';"
                    sql_root = sql_root.replace('<command>',command)
                    sql_root = sql_root.replace('<schema>',schema)
                    sql_root = sql_root.replace('<table>',prop_table)
                    sql_root = sql_root.replace('<mol_id>',str(mol_id))
                    # 
                    entry_intable = dict()
                    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                    cursor = connection.cursor()       
                    for prop_key,prop_val in properties.iteritems():
                        name_i = prop_name[prop_key]
                        desc_i = prop_desc[prop_key]
                        entry_intable[name_i] = 0 
                        if( len(name_i) and len(desc_i) ):
                            # See if property is already in table 
                            sql = "SELECT COUNT(1) FROM  <schema>.<table>  WHERE mol_id = <mol_id> AND property_name  = \'<property_name>\';"
                            sql = sql.replace('<schema>',schema)
                            sql = sql.replace('<table>',prop_table)                            
                            sql = sql.replace('<mol_id>',str(mol_id))
                            sql = sql.replace('<property_name>',name_i)
                            cursor.execute(sql)
                            entries = cursor.fetchone()[0]
                            connection.commit()
                            
                            if( entries == 0 ):
                                logger.info("  Property %s not in table will be inserted into"%(name_i))
                                sql = "INSERT INTO <schema>.<table> (mol_id,property_name,property_method,property_value) VALUES (%s,%s,%s,%s)"
                                sql = sql.replace('<schema>',schema)
                                sql = sql.replace('<table>',prop_table)
                                cursor.execute(sql,(mol_id,name_i,desc_i,prop_val))
                                connection.commit()
                            else:
                                logger.info("  Updating %s"%(name_i))
                                sql = copy.deepcopy(sql_root)
                                sql = sql.replace('<property_name>',str(name_i))
                                sql = sql.replace('<prop_val>',str(prop_val))
                                # insert_prop(schema,func_moltype,bb_i.properties['ctag'],prop_key,prop_val,command)
                                cursor.execute(sql)
                                connection.commit()
                    cursor.close()
                    connection.close()




def attachprepedAB(a_i,b_i,obopt=True):
    '''
    Attached prepped buildingblocks a and b 
    '''
    # 
    # Make copies of containers to modify
    # 
    bbC_i = copy.deepcopy(a_i)
    bbC_j = copy.deepcopy(b_i)
    # 
    # bbC_i.tag +=  b_i.tag
    # bbC_i.deptag = a_i.deptag +  b_i.deptag
    bbC_i.tag += b_i.tag
    bbC_i.properties['ctag'] +=  b_i.properties['ctag']
    bbC_i.properties['deptag'] += b_i.properties['deptag']
    bbC_i.properties['name'] += b_i.properties['name']
    bbC_i.properties['common_tag'] += b_i.properties['common_tag']
    bbC_i.properties['iupac'] += b_i.properties['iupac']
    # 
    # bbC_i.bblist = a_i.deptag +',' + b_i.deptag
    # 
    logger.info('- Structure %s - %s => %s (%s)'%(a_i.properties['ctag'],b_i.properties['ctag'],bbC_i.properties['ctag'],bbC_i.properties['common_tag']))
    #
    Xo_i = bbC_i.danglkey
    Xo_j = bbC_j.danglkey
    #
    # Shift  building block j to correct bond length 
    radii_i = bbC_i.particles[Xo_i].properties["cov_radii"]
    radii_j = bbC_j.particles[Xo_j].properties["cov_radii"]
    bond_vec = np.array([radii_i + radii_j,0.0,0.0])
    bbC_j.shift_pos(-1.0*bond_vec)
    #
    # Add j to i
    bbC_i += bbC_j
    #
    # Get updated atom keys to form bond 
    Xp_i = Xo_i                   
    Xp_j = bbC_i.keyupdate[Xo_j]
    #
    # Create bond between  X_i - X_j
    bond_ij = structure.Bond(Xp_i,Xp_j)
    bond_ij.properties['bondorder'] = 1
    bbC_i.add_bond(bond_ij)
    #
    # Remake neighbor list based on updated bonds 
    bbC_i.bonded_nblist = structure.NBlist() 
    bbC_i.bonded_nblist.build_nblist(bbC_i.particles,bbC_i.bonds )
    # Update number of attachment points
    bbC_i.calc_attachments()
    # Optimize with openbabel
    
    if( obopt ):
        # pybelmol = ob.pybel_opt(bbC_i)
        pybelmol =convert_struc_pymol(bbC_i)
        pybelmol.localopt()
        for pkey_i, particle_i  in bbC_i.particles.iteritems():
            a = pybelmol.atoms[pkey_i]
            bbC_i.positions[pkey_i][0] = a.vector.GetX()
            bbC_i.positions[pkey_i][1] = a.vector.GetY()
            bbC_i.positions[pkey_i][2] = a.vector.GetZ()

    #  Return final structure
    # 
    return bbC_i
def CreateDimers(p,table_a,table_b,calc_tag, backbone_type_name):
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()

    local = resource.Resource('local')
    local.load_json()
    peregrine = resource.Resource('peregrine')
    peregrine.load_json()
    bbdir ="/home/tkemper/BuildingBlocks-cand_opt" # options.bbdir  #
    bbdir = '/Users/tkemper/Projects/opv-project/BuildingBlocks-cand_opt'
    listprop_dic = {'term_list':'term_list','func_list':'func_list','sub_list':'sub_list','bond_list':'bond_list'}
    #
    schema = 'opv_cand'
    table_ab = 'candidate_oligomer'
    rotateA = True
    write_xyz = False
    regen = False 
    #
    # script_dir = '/Users/tkemper/Projects/opv-project/scripts'
    # os.chdir(script_dir)
    # 
    prop_name,prop_desc = set_prop_names_desc()
    #
    #
    #
    # Get mol_ids of A's
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()
    sql = "SELECT mol_id FROM  <schema>.<table>;"
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',table_a)
    cursor.execute(sql)
    a_molreturn = cursor.fetchall()
    connection.commit()
    cursor.close()
    connection.close()
    a_mol_ids = []
    for r in reversed(a_molreturn):
        (mol_id,) = r
        a_mol_ids.append(mol_id)
    p.barrier()
    logger.info("a_mol_ids %d "%(len(a_mol_ids)))
    #
    # Get mol_ids of B's
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()
    sql = "SELECT mol_id FROM  <schema>.<table>;"
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',table_b)
    cursor.execute(sql)
    b_molreturn = cursor.fetchall()
    connection.commit()
    cursor.close()
    connection.close()
    b_mol_ids = []
    for r in b_molreturn:
        (mol_id,) = r
        b_mol_ids.append(mol_id)
    p.barrier()

    logger.info("b_mol_ids %d "%(len(b_mol_ids)))

    
    angle_rad = 90.0*math.pi/180.0


    # Set position to attach building blocks 
    bbid_a = "T"
    n_a = 1
    Xn_a = 0  # Number of term atom in neighbor list of cap atom 
    bbid_b = "T"
    n_b = 0
    Xn_b = 0  # Number of term atom in neighbor list of cap atom
    align_a = 1
    align_b = -1

    if( write_xyz ):
        # INSERT needed oligomers
        dir_i = bbdir +"/"+backbone_type_name
        if ( not os.path.isdir(dir_i) ):
            print "Making %s "%(dir_i)
            os.mkdir(dir_i)
        os.chdir(dir_i)


    p.barrier()
    fg_list_prep = []
    a_mol_ids_p = p.splitListOnProcs(a_mol_ids)        
    p.barrier()
    #
    for a_mol_id in a_mol_ids_p:
        
        a_i = convert_belmol_strucC(schema,table_a,a_mol_id)
        
        logger.info("Processing building block A %s on proc %d "%(a_i.properties['ctag'],rank))
                    
        a_i.tag += "_"
        a_i.properties['common_tag'] += "_"
        a_i.properties['ctag'] += "_"
        a_i.properties['deptag'] += "_"
        a_i.properties['name']  += "_"
        a_i.properties['iupac']  += "-"
        # Find keys of attachment points 
        Rkey_a,Xkey_a = a_i.find_XR(bbid_a,n_a,Xn_a)
        a_i.align_bond(Rkey_a,Xkey_a)
        a_i.shift_pos(-1.0*a_i.positions[Xkey_a] )
        if( rotateA ):
            # Find attached heavy atoms
            # align fist heavy neighbor with y axis 
            for key_k in a_i.bonded_nblist.getnbs(Xkey_a):
                    particle_k = a_i.particles[key_k]
                    if( particle_k.properties['number'] != 1 ):                        
                            a_i.align_yaxis(key_k,angle_rad)
                            break
        #
        # Remove atoms at R in building block i
        a_i.del_particle(Rkey_a)
        # Set particle key with dangling bond 
        a_i.danglkey = a_i.keyupdate[Xkey_a]
        a_i.attach_p = a_i.keyupdate[Xkey_a]
        
        for b_mol_id in b_mol_ids:
            #
            b_i = convert_belmol_strucC(schema,table_b,b_mol_id)
            
            b_i.tag += "_"
            b_i.properties['common_tag'] += "_"
            b_i.properties['ctag'] += "_"
            b_i.properties['deptag'] += "_"
            b_i.properties['name']  += "_"
            b_i.properties['iupac']  += "-"

            if( a_i.properties['ctag'] !=  b_i.properties['ctag']  ):
                
                ctag =  a_i.properties['ctag'] +b_i.properties['ctag'] 
                #
                logger.info("Checking for  %s in %s.%s "%(ctag,schema,table_ab))
                #
                connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                cursor = connection.cursor()       
                sql = "SELECT COUNT(1) FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                sql = sql.replace('<schema>',schema)
                sql = sql.replace('<table>',table_ab)
                sql = sql.replace('<ctag>',ctag)
                cursor.execute(sql)
                cand_id_cnt = cursor.fetchone()[0]
                connection.commit()
                cursor.close()
                connection.close()
                logger.info("cand_id_cnt %d "%(cand_id_cnt))
                if( cand_id_cnt == 0 or regen ):
                    # Building  oligomer  
                    logger.info("Building oligomer  %s "%(ctag))
                    # bb_prepped = copy.deepcopy(bb_i)
                    #
                    # Find keys of attachment points
                    #
                    Rkey_b,Xkey_b = b_i.find_XR(bbid_b,n_b,Xn_b)
                    # Align building blocks along bonds of attachment atoms 
                    b_i.align_bond(Xkey_b,Rkey_b)
                    # Remove atoms at R in building block i
                    b_i.del_particle(Rkey_b)
                    # Set particle key with dangling bond 
                    b_i.danglkey = b_i.keyupdate[Xkey_b]
                    b_i.attach_p = b_i.keyupdate[Xkey_b]
                    # 
                    oligomer_i = attachprepedAB(a_i,b_i)
                    # 
                    # Process oligomer 
                    oligomer_i.proc_bbid()
                    oligomer_i.properties['bond_list'] = []
                    for bkey,b_i in oligomer_i.bonds.iteritems():
                        oligomer_i.properties['bond_list'].append([b_i.pkey1,b_i.pkey2])
                    # Get mol str     
                    xyz_str = oligomer_i.write_xyz_str()
                    pybelmol = pybel.readstring('xyz',xyz_str)
                    #
                    molstr  = pybelmol.write('mol')
                    oligomer_i.properties['molstr'] = molstr
                    #
                    properties = ob_properties(pybelmol)
                    #
                    # Check properties from openbabel
                    #                
                    add_oligomer = True 
                    if( oligomer_i.n_particles != properties['atoms'] ):
                            print "%s has an Error in atoms "%(oligomer_i.properties['ctag'])
                            add_oligomer = False
                            sys.exit(1)
                    if( oligomer_i.n_bonds != properties['bonds'] ):
                            print "%s has an Error in bonds "%(oligomer_i.properties['ctag'])
                            print oligomer_i.n_bonds , properties['bonds']
                            add_oligomer = False
                            sys.exit(2)
                    if( add_oligomer and cand_id_cnt == 0  ):
                        # Create entry 
                        logger.info("Creating new entry in %s.%s  proc %d  "%(schema,table_ab,rank))
                        connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                        cursor = connection.cursor()       
                        sql = "INSERT INTO <schema>.<table> (tag,ctag,iupac,name,common_tag,mol,backbone_type_name) VALUES (%s,%s,%s,%s,%s,%s,%s)"
                        sql = sql.replace('<schema>',schema)
                        sql = sql.replace('<table>',table_ab)
                        #print "!!debug0!! ",sql,oligomer_i.properties['ctag'],oligomer_i.properties['iupac'],oligomer_i.properties['name'],oligomer_i.properties['common_tag'],oligomer_i.properties['molstr']
                        cursor.execute(sql,(oligomer_i.tag,oligomer_i.properties['ctag'],oligomer_i.properties['iupac'],oligomer_i.properties['name'],oligomer_i.properties['common_tag'],oligomer_i.properties['molstr'],backbone_type_name))
                        connection.commit()
                        cursor.close()
                        connection.close()
                    else:
                        logger.warning(" Oligomer %s not INSERTED ")
                        #
                    if( add_oligomer ):
                        #if( add_props ):
                        logger.info("Checking for  %s in %s.%s "%(oligomer_i.properties['ctag'],schema,table_ab))
                        connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                        cursor = connection.cursor()       
                        sql = "SELECT COUNT(1) FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                        sql = sql.replace('<schema>',schema)
                        sql = sql.replace('<table>',table_ab)
                        sql = sql.replace('<ctag>',oligomer_i.properties['ctag'])
                        cursor.execute(sql)
                        candidate_oligomer_id_cnt = cursor.fetchone()[0]
                        connection.commit()
                        cursor.close()
                        connection.close()
                        if( candidate_oligomer_id_cnt == 0 ):
                            logger.warning(" ctag %s not found in %s.%s "%(oligomer_i.properties['ctag'],schema,table_ab))
                        else:
                            # Get candidate_oligomer_id
                            connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                            cursor = connection.cursor()       
                            sql = "SELECT candidate_oligomer_id FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                            sql = sql.replace('<schema>',schema)
                            sql = sql.replace('<table>',table_ab)
                            sql = sql.replace('<ctag>',oligomer_i.properties['ctag'])
                            cursor.execute(sql)
                            candidate_oligomer_id = cursor.fetchone()[0]
                            connection.commit()
                            cursor.close()
                            connection.close()
                            logger.info("Found %s with candidate_oligomer_id %d "%(oligomer_i.properties['ctag'],candidate_oligomer_id))
                            # Add lists 
                            command = 'UPDATE'
                            sql_root = "<command> <schema>.<table> SET <prop_key> = <prop_val> WHERE candidate_oligomer_id = '<candidate_oligomer_id>';"
                            sql_root = sql_root.replace('<command>',command)
                            sql_root = sql_root.replace('<schema>',schema)
                            sql_root = sql_root.replace('<table>',table_ab)
                            sql_root = sql_root.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
                            connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                            cursor = connection.cursor()                    
                            if( regen ):
                                # Update mol string as other properties should not have changed
                                prop_key = 'mol'
                                pkey = 'molstr'
                                prop_val = '\'%s\''%oligomer_i.properties[pkey]
                                logger.info("  UPDATE %s "%(prop_key))
                                sql = copy.deepcopy(sql_root)
                                sql = sql.replace('<prop_key>',prop_key)
                                sql = sql.replace('<prop_val>',str(prop_val))
                                #print "!!debug2.5!! ",sql
                                cursor.execute(sql)
                            for prop_key,pkey in listprop_dic.iteritems():
                                prop_val = '\'%s\''%oligomer_i.properties[pkey]
                                prop_val = prop_val.replace(']','}')
                                prop_val = prop_val.replace('[','{')
                                logger.info("  UPDATE %s "%(prop_key))
                                sql = copy.deepcopy(sql_root)
                                sql = sql.replace('<prop_key>',prop_key)
                                sql = sql.replace('<prop_val>',str(prop_val))
                                #print "!!debug2!! ",sql
                                cursor.execute(sql)
                            connection.commit()
                            cursor.close()
                            connection.close()
                            # Add properties 
                            prop_table = "%s_calc_prop"%(table_ab)
                            logger.info("Adding properties for candidate_oligomer_id %d ctag %s to %s"%(candidate_oligomer_id,oligomer_i.properties['ctag'],prop_table))
                            #
                            command = 'UPDATE'
                            sql_root = "<command> <schema>.<table> SET property_value = <prop_val> WHERE candidate_oligomer_id = <candidate_oligomer_id> AND property_name = \'<property_name>\';"
                            sql_root = sql_root.replace('<command>',command)
                            sql_root = sql_root.replace('<schema>',schema)
                            sql_root = sql_root.replace('<table>',prop_table)
                            sql_root = sql_root.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
                            # 
                            connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
                            cursor = connection.cursor()       
                            for prop_key,prop_val in properties.iteritems():
                                name_i = prop_name[prop_key]
                                desc_i = prop_desc[prop_key]
                                if( len(name_i) and len(desc_i) ):
                                    # See if property is already in table 
                                    sql = "SELECT COUNT(1) FROM  <schema>.<table>  WHERE candidate_oligomer_id = <candidate_oligomer_id> AND property_name  = \'<property_name>\';"
                                    sql = sql.replace('<schema>',schema)
                                    sql = sql.replace('<table>',prop_table)                            
                                    sql = sql.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
                                    sql = sql.replace('<property_name>',name_i)
                                    cursor.execute(sql)
                                    entries = cursor.fetchone()[0]
                                    connection.commit()

                                    if( entries == 0 ):
                                        logger.info("   INSERT property %s "%(name_i))
                                        sql = "INSERT INTO <schema>.<table> (candidate_oligomer_id,property_name,property_method,property_value) VALUES (%s,%s,%s,%s)"
                                        sql = sql.replace('<schema>',schema)
                                        sql = sql.replace('<table>',prop_table)
                                        #print "!!debug 3 !!",sql,candidate_oligomer_id,name_i,desc_i,prop_val
                                        cursor.execute(sql,(candidate_oligomer_id,name_i,desc_i,prop_val))
                                        connection.commit()
                                    else:
                                        logger.info("  UPDATE %s"%(name_i))
                                        sql = copy.deepcopy(sql_root)
                                        sql = sql.replace('<property_name>',str(name_i))
                                        sql = sql.replace('<prop_val>',str(prop_val))
                                        #print "!!debug 4 !!",sql
                                        cursor.execute(sql)
                                        connection.commit()
                            cursor.close()
                            connection.close()

                else:
                    logger.info("Entry in %s.%s found skipping generation "%(schema,table_ab))
            else:
                logger.info(" Not including A (%s) == B (%s) "%( a_i.properties['ctag'],b_i.properties['ctag']))
        #
        p.barrier()
        # sys.exit("mpi test 2394872938")
        
    return

def DelCand(p,backbone_type_name):
    schema = 'opv_cand'
    table_ab = 'candidate_oligomer'
    prop_table = "%s_calc_prop"%(table_ab)
    #
    # Get mol_ids of B's
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()
    sql = "SELECT candidate_oligomer_id FROM  <schema>.<table>;"
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',table_ab)
    cursor.execute(sql)
    cand_id_return = cursor.fetchall()
    connection.commit()
    cand_ids = []
    for r in cand_id_return:
        (mol_id,) = r
        cand_ids.append(mol_id)
    p.barrier()
    for candidate_oligomer_id in cand_ids:
        sql = "DELETE FROM <schema>.<table> WHERE candidate_oligomer_id = <candidate_oligomer_id>;"
        sql = sql.replace('<schema>',schema)
        sql = sql.replace('<table>',prop_table)
        sql = sql.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
        #print sql
        cursor.execute(sql)
        connection.commit()
        sql = "DELETE FROM <schema>.<table> WHERE candidate_oligomer_id = <candidate_oligomer_id>;"
        sql = sql.replace('<schema>',schema)
        sql = sql.replace('<table>',table_ab)
        sql = sql.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
        print sql
        #cursor.execute(sql)
        connection.commit()
        
    
    cursor.close()
    connection.close()


def CreateDimers2(p,table_a,table_b,calc_tag, backbone_type_name):
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    #
    # MPI setup
    #
    rank = p.getRank()
    size = p.getCommSize()

    local = resource.Resource('local')
    local.load_json()
    peregrine = resource.Resource('peregrine')
    peregrine.load_json()
    bbdir ="/home/tkemper/BuildingBlocks-cand_opt" # options.bbdir  #
    bbdir = '/Users/tkemper/Projects/opv-project/BuildingBlocks-cand_opt'
    listprop_dic = {'term_list':'term_list','func_list':'func_list','sub_list':'sub_list','bond_list':'bond_list'}
    #
    schema = 'opv_cand'
    table_ab = 'candidate_oligomer'
    rotateA = True
    write_xyz = False
    regen = False 
    n_generated = 0 
    n_listinserts = 0 
    n_propinserts = 0 
    n_candinserts = 0 
    #
    # script_dir = '/Users/tkemper/Projects/opv-project/scripts'
    # os.chdir(script_dir)
    # 
    prop_name,prop_desc = set_prop_names_desc()
    #
    #
    #
    # Get mol_ids of A's
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov' port=5432 ")
    cursor = connection.cursor()
    sql = "SELECT mol_id FROM  <schema>.<table>;"
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',table_a)
    cursor.execute(sql)
    a_molreturn = cursor.fetchall()
    #connection.commit()
    #
    a_mol_ids = []
    for r in reversed(a_molreturn):
        (mol_id,) = r
        a_mol_ids.append(mol_id)
    p.barrier()
    #a_mol_ids = a_mol_ids[100:124]    
    logger.info("a_mol_ids %d "%(len(a_mol_ids)))
    #
    # Get mol_ids of B's
    #
    sql = "SELECT mol_id FROM  <schema>.<table>;"
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',table_b)
    cursor.execute(sql)
    b_molreturn = cursor.fetchall()
    #connection.commit()
    
    b_mol_ids = []
    for r in b_molreturn:
        (mol_id,) = r
        b_mol_ids.append(mol_id)
    p.barrier()
    #b_mol_ids = b_mol_ids[570:574]    
    #b_mol_ids = [len(b_mol_ids)-1,len(b_mol_ids)-2]
    logger.info("b_mol_ids %d "%(len(b_mol_ids)))

    
    angle_rad = 90.0*math.pi/180.0


    # Set position to attach building blocks 
    bbid_a = "T"
    n_a = 1
    Xn_a = 0  # Number of term atom in neighbor list of cap atom 
    bbid_b = "T"
    n_b = 0
    Xn_b = 0  # Number of term atom in neighbor list of cap atom
    align_a = 1
    align_b = -1

    if( write_xyz ):
        # INSERT needed oligomers
        dir_i = bbdir +"/"+backbone_type_name
        if ( not os.path.isdir(dir_i) ):
            print "Making %s "%(dir_i)
            os.mkdir(dir_i)
        os.chdir(dir_i)

    fout = open('bad_oligomers.csv','w')
    writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
    writer.writerow( ['type',"tag"] )
    fout.close()  

    p.barrier()
    fg_list_prep = []
    a_mol_ids_p = p.splitListOnProcs(a_mol_ids)        
    p.barrier()
    #
    for a_mol_id in a_mol_ids_p:
        
        a_i = convert_belmol_strucC2(cursor,schema,table_a,a_mol_id)
        
        logger.info("Processing building block A %s on proc %d "%(a_i.properties['ctag'],rank))
                    
        a_i.tag += "_"
        a_i.properties['common_tag'] += "_"
        a_i.properties['ctag'] += "_"
        a_i.properties['deptag'] += "_"
        a_i.properties['name']  += "_"
        a_i.properties['iupac']  += "-"
        # Find keys of attachment points 
        Rkey_a,Xkey_a = a_i.find_XR(bbid_a,n_a,Xn_a)
        a_i.align_bond(Rkey_a,Xkey_a)
        a_i.shift_pos(-1.0*a_i.positions[Xkey_a] )
        if( rotateA ):
            # Find attached heavy atoms
            # align fist heavy neighbor with y axis 
            for key_k in a_i.bonded_nblist.getnbs(Xkey_a):
                    particle_k = a_i.particles[key_k]
                    if( particle_k.properties['number'] != 1 ):                        
                            a_i.align_yaxis(key_k,angle_rad)
                            break
        #
        # Remove atoms at R in building block i
        a_i.del_particle(Rkey_a)
        # Set particle key with dangling bond 
        a_i.danglkey = a_i.keyupdate[Xkey_a]
        a_i.attach_p = a_i.keyupdate[Xkey_a]
        
        for b_mol_id in b_mol_ids:
            #
            b_i = convert_belmol_strucC2(cursor,schema,table_b,b_mol_id)
            
            b_i.tag += "_"
            b_i.properties['common_tag'] += "_"
            b_i.properties['ctag'] += "_"
            b_i.properties['deptag'] += "_"
            b_i.properties['name']  += "_"
            b_i.properties['iupac']  += "-"

            if( a_i.properties['ctag'] !=  b_i.properties['ctag']  ):
                
                ctag =  a_i.properties['ctag'] +b_i.properties['ctag'] 
                #
                logger.info("Checking for  %s in %s.%s "%(ctag,schema,table_ab))
                #
                sql = "SELECT COUNT(1) FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                sql = sql.replace('<schema>',schema)
                sql = sql.replace('<table>',table_ab)
                sql = sql.replace('<ctag>',ctag)
                cursor.execute(sql)
                cand_id_cnt = cursor.fetchone()[0]
                #connection.commit()
                logger.info("cand_id_cnt %d "%(cand_id_cnt))
                if( cand_id_cnt == 0 or regen ):
                    # Building  oligomer  
                    logger.info("Building oligomer %s on proc %d "%(ctag,rank))
                    # bb_prepped = copy.deepcopy(bb_i)
                    #
                    # Find keys of attachment points
                    #
                    Rkey_b,Xkey_b = b_i.find_XR(bbid_b,n_b,Xn_b)
                    # Align building blocks along bonds of attachment atoms 
                    b_i.align_bond(Xkey_b,Rkey_b)
                    # Remove atoms at R in building block i
                    b_i.del_particle(Rkey_b)
                    # Set particle key with dangling bond 
                    b_i.danglkey = b_i.keyupdate[Xkey_b]
                    b_i.attach_p = b_i.keyupdate[Xkey_b]
                    # 
                    oligomer_i = attachprepedAB(a_i,b_i)
                    # 
                    # Process oligomer 
                    oligomer_i.proc_bbid()
                    oligomer_i.properties['bond_list'] = []
                    for bkey,b_i in oligomer_i.bonds.iteritems():
                        oligomer_i.properties['bond_list'].append([b_i.pkey1,b_i.pkey2])
                    # Get mol str     
                    xyz_str = oligomer_i.write_xyz_str()
                    pybelmol = pybel.readstring('xyz',xyz_str)
                    #
                    molstr  = pybelmol.write('mol')
                    oligomer_i.properties['molstr'] = molstr
                    #
                    properties = ob_properties(pybelmol)
                    #
                    # Check properties from openbabel
                    #                
                    add_oligomer = True 
                    if( oligomer_i.n_particles != properties['atoms'] ):
                            print "%s has an Error in atoms "%(oligomer_i.properties['ctag'])
                            add_oligomer = False
                            oligomer_i.write_xyz()
                            fout = open('bad_oligomers.csv','a')
                            writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                            writer.writerow( ['atoms',oligomer_i.tag] )
                            fout.close()  
                            
                            #sys.exit(1)
                    if( oligomer_i.n_bonds != properties['bonds'] ):
                            print "%s has an Error in bonds "%(oligomer_i.properties['ctag'])
                            print oligomer_i.n_bonds , properties['bonds']
                            
                            oligomer_i.write_xyz()
                            fout = open('bad_oligomers.csv','a')
                            writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                            writer.writerow( ['bonds',oligomer_i.tag] )
                            fout.close()   
                            add_oligomer = False
                            #sys.exit(2)
                    if( add_oligomer and cand_id_cnt == 0  ):
                        # Create entry 
                        logger.info("Creating new entry in %s.%s  proc %d  "%(schema,table_ab,rank))
                        
                        sql = "INSERT INTO <schema>.<table> (tag,ctag,iupac,name,common_tag,mol,backbone_type_name) VALUES (%s,%s,%s,%s,%s,%s,%s)"
                        sql = sql.replace('<schema>',schema)
                        sql = sql.replace('<table>',table_ab)
                        #print "!!debug0!! ",sql,oligomer_i.properties['ctag'],oligomer_i.properties['iupac'],oligomer_i.properties['name'],oligomer_i.properties['common_tag'],oligomer_i.properties['molstr']
                        cursor.execute(sql,(oligomer_i.tag,oligomer_i.properties['ctag'],oligomer_i.properties['iupac'],oligomer_i.properties['name'],oligomer_i.properties['common_tag'],oligomer_i.properties['molstr'],backbone_type_name))
                        connection.commit()
                        n_candinserts += 1 
                        #
                    else:
                        logger.warning(" Oligomer %s not INSERTED ")
                        #
                    if( add_oligomer ):
                        #if( add_props ):
                        logger.info("Checking for  %s in %s.%s "%(oligomer_i.properties['ctag'],schema,table_ab))
                        #     
                        sql = "SELECT COUNT(1) FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                        sql = sql.replace('<schema>',schema)
                        sql = sql.replace('<table>',table_ab)
                        sql = sql.replace('<ctag>',oligomer_i.properties['ctag'])
                        cursor.execute(sql)
                        candidate_oligomer_id_cnt = cursor.fetchone()[0]
                        #connection.commit()
                        #
                        if( candidate_oligomer_id_cnt == 0 ):
                            logger.warning(" ctag %s not found in %s.%s "%(oligomer_i.properties['ctag'],schema,table_ab))
                        else:
                            # Get candidate_oligomer_id
                            #    
                            sql = "SELECT candidate_oligomer_id FROM <schema>.<table> WHERE ctag = \'<ctag>\';"
                            sql = sql.replace('<schema>',schema)
                            sql = sql.replace('<table>',table_ab)
                            sql = sql.replace('<ctag>',oligomer_i.properties['ctag'])
                            cursor.execute(sql)
                            candidate_oligomer_id = cursor.fetchone()[0]
                            #connection.commit()
                            #
                            logger.info("Found %s with candidate_oligomer_id %d "%(oligomer_i.properties['ctag'],candidate_oligomer_id))
                            # Add lists 
                            command = 'UPDATE'
                            sql_root = "<command> <schema>.<table> SET <prop_key> = <prop_val> WHERE candidate_oligomer_id = '<candidate_oligomer_id>';"
                            sql_root = sql_root.replace('<command>',command)
                            sql_root = sql_root.replace('<schema>',schema)
                            sql_root = sql_root.replace('<table>',table_ab)
                            sql_root = sql_root.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
                            #                   
                            if( regen ):
                                # Update mol string as other properties should not have changed
                                prop_key = 'mol'
                                pkey = 'molstr'
                                prop_val = '\'%s\''%oligomer_i.properties[pkey]
                                logger.info("  UPDATE %s "%(prop_key))
                                sql = copy.deepcopy(sql_root)
                                sql = sql.replace('<prop_key>',prop_key)
                                sql = sql.replace('<prop_val>',str(prop_val))
                                #print "!!debug2.5!! ",sql
                                cursor.execute(sql)
                            for prop_key,pkey in listprop_dic.iteritems():
                                prop_val = '\'%s\''%oligomer_i.properties[pkey]
                                prop_val = prop_val.replace(']','}')
                                prop_val = prop_val.replace('[','{')
                                logger.info("  UPDATE %s "%(prop_key))
                                sql = copy.deepcopy(sql_root)
                                sql = sql.replace('<prop_key>',prop_key)
                                sql = sql.replace('<prop_val>',str(prop_val))
                                #print "!!debug2!! ",sql
                                cursor.execute(sql)
                            #connection.commit()
                            n_listinserts += 1 
                            #
                            # Add properties 
                            prop_table = "%s_calc_prop"%(table_ab)
                            logger.info("Adding properties for candidate_oligomer_id %d ctag %s to %s"%(candidate_oligomer_id,oligomer_i.properties['ctag'],prop_table))
                            #
                            command = 'UPDATE'
                            sql_root = "<command> <schema>.<table> SET property_value = <prop_val> WHERE candidate_oligomer_id = <candidate_oligomer_id> AND property_name = \'<property_name>\';"
                            sql_root = sql_root.replace('<command>',command)
                            sql_root = sql_root.replace('<schema>',schema)
                            sql_root = sql_root.replace('<table>',prop_table)
                            sql_root = sql_root.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
                            # 
                            for prop_key,prop_val in properties.iteritems():
                                name_i = prop_name[prop_key]
                                desc_i = prop_desc[prop_key]
                                if( len(name_i) and len(desc_i) ):
                                    # See if property is already in table 
                                    sql = "SELECT COUNT(1) FROM  <schema>.<table>  WHERE candidate_oligomer_id = <candidate_oligomer_id> AND property_name  = \'<property_name>\';"
                                    sql = sql.replace('<schema>',schema)
                                    sql = sql.replace('<table>',prop_table)                            
                                    sql = sql.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
                                    sql = sql.replace('<property_name>',name_i)
                                    cursor.execute(sql)
                                    entries = cursor.fetchone()[0]
                                    #connection.commit()
                                    #
                                    if( entries == 0 ):
                                        logger.info("   INSERT property %s "%(name_i))
                                        sql = "INSERT INTO <schema>.<table> (candidate_oligomer_id,property_name,property_method,property_value) VALUES (%s,%s,%s,%s)"
                                        sql = sql.replace('<schema>',schema)
                                        sql = sql.replace('<table>',prop_table)
                                        #print "!!debug 3 !!",sql,candidate_oligomer_id,name_i,desc_i,prop_val
                                        cursor.execute(sql,(candidate_oligomer_id,name_i,desc_i,prop_val))
                                        connection.commit()
                                        n_propinserts +=1 
                                    else:
                                        logger.info("  UPDATE %s"%(name_i))
                                        sql = copy.deepcopy(sql_root)
                                        sql = sql.replace('<property_name>',str(name_i))
                                        sql = sql.replace('<prop_val>',str(prop_val))
                                        #print "!!debug 4 !!",sql
                                        cursor.execute(sql)
                                        #connection.commit()
                            connection.commit()

                else:
                    logger.info("Entry in %s.%s found skipping generation "%(schema,table_ab))
            else:
                logger.info(" Not including A (%s) == B (%s) "%( a_i.properties['ctag'],b_i.properties['ctag']))
        #
        p.barrier()
        # sys.exit("mpi test 2394872938")
    # Close Connection
    p.barrier()
    cursor.close()
    connection.close()
    p.barrier()
    #
    #if( rank == 0 ):
    logger.info("structures generated %d on proc %d "%(n_generated,rank))
    logger.info("lists inserted %d on proc %d "%(n_listinserts,rank))
    logger.info("structures inserted %d on proc %d "%(n_candinserts,rank))
    logger.info("properties inserted %d on proc %d "%(n_propinserts,rank))
    p.barrier() 
    #
    return 

def DelCand(p,backbone_type_name):
    schema = 'opv_cand'
    table_ab = 'candidate_oligomer'
    prop_table = "%s_calc_prop"%(table_ab)
    #
    # Get mol_ids of B's
    #
    connection = psycopg2.connect("dbname='opv' user='tkemper' password='M0l3cul316' host='kestrel.hpc.nrel.gov'")
    cursor = connection.cursor()
    sql = "SELECT candidate_oligomer_id FROM  <schema>.<table>;"
    sql = sql.replace('<schema>',schema)
    sql = sql.replace('<table>',table_ab)
    cursor.execute(sql)
    cand_id_return = cursor.fetchall()
    connection.commit()
    cand_ids = []
    for r in cand_id_return:
        (mol_id,) = r
        cand_ids.append(mol_id)
    p.barrier()
    for candidate_oligomer_id in cand_ids:
        sql = "DELETE FROM <schema>.<table> WHERE candidate_oligomer_id = <candidate_oligomer_id>;"
        sql = sql.replace('<schema>',schema)
        sql = sql.replace('<table>',prop_table)
        sql = sql.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
        #print sql
        cursor.execute(sql)
        connection.commit()
        sql = "DELETE FROM <schema>.<table> WHERE candidate_oligomer_id = <candidate_oligomer_id>;"
        sql = sql.replace('<schema>',schema)
        sql = sql.replace('<table>',table_ab)
        sql = sql.replace('<candidate_oligomer_id>',str(candidate_oligomer_id))
        print sql
        #cursor.execute(sql)
        connection.commit()
        
    
    cursor.close()
    connection.close()

        

    
if __name__=="__main__":

    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--bbdir", dest="bbdir", type="string", default="/home/tkemper/BuildingBlocks-cand_opt", help="BuildingBlocks dir")
    parser.add_option("--ppn", dest="ppn", type="int", default=24, help="Processors per node for nwchem calculations ")
    parser.add_option("-a","--table_a", dest="table_a", type="string", default="donors_func_fullsym", help="table of building block for A in A-B dimer")
    parser.add_option("-b","--table_b", dest="table_b", type="string", default="acceptors_func_fullsym", help="table of building block for B in A-B dimer")
    parser.add_option("--bbtype", dest="bbtype", type="string", default="DA", help="backbone_type_name")
    (options, args) = parser.parse_args()
     
    #comm = MPI.COMM_WORLD
    #size = comm.Get_size()
    #rank = comm.Get_rank()

    #print "I am %d in a world of %d" % (rank,size)     
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    #
    # MPI setup
    #    
    rank = p.getRank()
    size = p.getCommSize()
    print " running on proc %d/%d"%(rank,size)
    # sys.exit("04293ur023")

    if( len(args) < 1 ):
        calc_tag = 'BDT_tddft_ob'
    else:
        calc_tag =  args[0]

    hdlr = logging.FileHandler('%s.log'%(calc_tag),mode='w')
    hdlr.setLevel(logging.DEBUG)
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)

    start_time = datetime.now()
    if( rank == 0 ):
        logger.info('Started %s '%(start_time))
    # Set up resources
    # set_res(calc_tag,options)

    # add_struc(calc_tag,options)
    # cal_obprops(calc_tag,options)
    # func_all(calc_tag,options,p)
    backbone_type_name = options.bbtype
    CreateDimers2(p,options.table_a,options.table_b,calc_tag,backbone_type_name)
    p.barrier()
    # DelCand(p,backbone_type_name)

    finish_time = datetime.now()
    delt_t = finish_time - start_time
    if( rank == 0 ):
        logger.info('Finished %s in %f seconds '%(finish_time,delt_t.seconds))
    p.barrier()
    
