#! /usr/bin/env python
"""
Run the simulaitons in a project 
"""



__author__ = "Travis W. Kemper"
__version__ = "0.1"
__email__ = "travkemp@gmail.com"
__status__ = "Alpha"


import   os, os.path, sys , copy ,shutil, logging, math, json, csv 
import numpy as np
from datetime import datetime
from optparse import OptionParser

from streamm import project


def run_calc(proj_tag,options):
    #
        
    logging.info('Running project %s  on  %s '%(proj_tag,datetime.now()))        
        
    proj_i = project.Project(calc_tag)
    proj_i.load_json()
    for calc_key,calc_i in proj_i.calculations.iteritems():
        print "Running %s "%calc_i.tag
        os.chdir(calc_i.dir['scratch'])
        calc_i.check()
        if( calc_i.meta['status'] == 'running' ):
            calc_i.meta['status'] = 'written'
        calc_i.run()
        calc_i.check()
        os.chdir(calc_i.dir['home'])

    
if __name__=="__main__":
    
    usage = "usage: %prog tag \n"
    parser = OptionParser(usage=usage)
    
    (options, args) = parser.parse_args()

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # 
    # logging.basicConfig(level=logging.DEBUG,
    #                format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    #                datefmt='%m-%d %H:%M',
    #                filemode='w')

    if( len(args) < 1 ):
        calc_tag = 'chain_sp'
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
    
    run_calc(calc_tag,options)

    finish_time = datetime.now()
    delt_t = finish_time - start_time
    #if( rank == 0 ):
    logger.info('Finished %s in %f seconds '%(finish_time,delt_t.seconds))
