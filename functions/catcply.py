import sys 
from copy import deepcopy
import numpy as np


import mpiBase
from buildingblocks import Buildingblock
# from buildingblocks import segment
# from simulation import siminfo
# from material import materialinfo 

import time, datetime



def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-o","--output_id", dest="output_id", type="string", default="", help="id")
    parser.add_option("-r","--repeat", dest="repeat", type="int", default=1, help="number of repeat untis ")    
    parser.add_option("-b", "--bbdir", dest="bbdir",  type="string", default=".",help="where to look for the building blocks")
    
    (options, args) = parser.parse_args()
        
    return options, args

def tokenize(s):
    import re
    #  print "tokenizing ", s
    p1 = re.split("\(([^)]*)\)", s)
    # splits "D (R1 R2 ) A (R4)" into D, R1 R2, A, R4
    #    print "p1 = ", p1
    allp = []
    even = False
    for p in p1:
        if even:
            allp.append ('(')
        p2 = p.split()
        #        print "p = ", p
        #        print "p2 = ", p2
        for g in p2:
            allp.append(g)
        if even:
            allp.append(')')
        if even:
            even=False
        else:
            even=True
    return allp

def write_log_line(log_out,rank,log_line):
    """
    write log line 
    """
    if( rank == 0 ):
        log_out.write(log_line)


def main():
    """
    Read in data file and replicate it 
    """

    #        
    # Read options 
    #
    options, args = get_options()
    #
    # Formated ouput varables
    #
    sperator_line = "\n---------------------------------------------------------------------"
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    p.verbose = False 
    
    rank = p.getRank()
    size = p.getCommSize()
    #
    # Read in cply file into a structure container object
    #
    if( rank == 0 ):
        # Open log files 
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w") 

    # Notes from peter's code
    # main
    #   > gen_struct
    #      > build_from_str(bblocks, input_str, options)
    input_str = args[0]
    if( rank == 0 ):
        log_line = "setting up files for input string = {} ".format( input_str)
        log_out.write(log_line)
    
    bb_id_list = tokenize(input_str)
    
    bb_list = [] # list of buildingblock objects
    
    set_func = False
    bb_cnt = -1
    for bb_id in bb_id_list:
        if( bb_id == "(" ):
            set_func = True
            bb_unit_i = bb_list[bb_cnt]
            bb_unit_i.func_cnt = 0
            
            
        if( bb_id == ")" ):
            set_func = False 
        
        if( bb_id != ")" and  bb_id != "(" ):
            bb_cnt += 1 
            cply_file = "{}/{}.cply".format(options.bbdir,bb_id)
            if( rank == 0 ):
                log_line = "Reading {}".format(cply_file)
                log_out.write(log_line)
            bb_o = Buildingblock(name=bb_id,verbose=False)
            bb_o.read_cply(cply_file)

            
            # If bonds are in file
            if( len(bb_o.bondC) > 0 ):
                bb_o.bondC_nblist()
            else:
                bb_o.build_bonded_nblist(max_nn=12.0,radii_buffer=1.25)
                bb_o.nblist_bonds()

            bb_o.ptclC.guess_radii()

            # Set type bassed on number of connection points in cply file
            # or if read in inbetween brackets 
            bb_o.set_type()
            if(set_func ):
                bb_o.type = "func"
                bb_unit_i.func_cnt += 1                
            
            if( rank == 0 ):
                log_line = "  Type {}  ".format(bb_o.type)
                log_line += "  with {} connections ".format(bb_o.connect_cnt)
                log_line += "  with {} func connections ".format(bb_o.func_connect_cnt)
                log_out.write(log_line)
            #bb_o.align_termcaps()
            

            bb_o.verbose = options.verbose
            
            bb_list.append(deepcopy(bb_o))
                
            del bb_o
    #
    #  Set tags
    #
    tag =""
    for bb_i in bb_list:
        # segment_i = segment(bb_i.name)
        
        if( bb_i.type == "unit" and bb_i.func_cnt > 0 ):
            tag += "{}_".format(bb_i.name)
        else: 
            tag += "{}".format(bb_i.name)

    # matinfo = materialinfo(tag)
            
    # Add functinoal groups
            
    bb_units = []
    func_cnt = 0
    unit_found = False

    write_log_line(log_out,rank,"Adding functional groups to units ")
    
    for bb_i in bb_list:
        if( bb_i.type == "unit" ):
            bb_unit_i = bb_i
            # Initialize functional connection point
            #  This will be increased during iadd
            unit_found = True 
            bb_unit_i.connectionpoint = -1
            n_func = bb_unit_i.func_cnt #bb_unit_i.func_connect_cnt
            func_cnt = 0
            write_log_line(log_out,rank,"  Initializing repeat unit {} with {} functional points ".format(bb_i.name,n_func))
            segment_i = {} #segment(str(col[1]))
                           # segment_i = segment()
            segment_i['tag'] = str(bb_i.name)
            segment_i['unit'] = str(bb_i.name)
            segment_i['func'] = []
            #segment_i['segment'] = {}
            
        if( bb_i.type == "func" ):
            func_cnt += 1
            write_log_line(log_out,rank,"  Adding the {} functional group {} to repeat unit ".format(func_cnt,bb_i.name))
            bb_unit_i += bb_i
            segment_i['func'].append( str(bb_i.name) )
            
        if( unit_found ):
            if( func_cnt == n_func ):
                if( options.verbose ):
                    log_line = "\n  All {}  functinoal groups have been added to repeat unit ".format(func_cnt)
                    log_line +=  sperator_line
                    write_log_line(log_out,rank,log_line)
                bb_unit_i.segments.append(segment_i)
                bb_units.append(deepcopy(bb_unit_i))
                unit_found = False 

    
    # Build repeat unit
    write_log_line(log_out,rank,"Concatenating units to build repeat unit")
    bb_chain = Buildingblock(verbose=False)
    bb_chain.verbose = False 
    for bb_i in bb_units:
        bb_chain += bb_i

    # Replicate repeat unit 
    if( options.repeat > 1 ):
        tag += "_n{}".format(options.repeat)
        bb_i = deepcopy(bb_chain)
        for n in range(2,options.repeat+1):
            write_log_line(log_out,rank, " Adding repeat  unit {}".format(n))
            bb_chain += deepcopy(bb_i)
            
    # Add Terminals
    if(  bb_list[0].type == "term" ):
        write_log_line(log_out,rank, " Adding terminal in first position ")
        bb_f = deepcopy(bb_list[0])
        bb_f += bb_chain
        
    else:
        bb_f = bb_chain
        
    if( bb_list[-1].type == "term" ):
        write_log_line(log_out,rank, " Adding terminal in last position ")
        bb_f += bb_list[-1]

    bb_f.verbose = False 
    # Set chain to 1 since all the same chain
    bb_f.set_tag("chain",1)

    if( len(options.output_id) ):
        tag = options.output_id
    
    # Print final structure 
    write_log_line(log_out,rank, str(bb_f))
    print  bb_f
    
    comment = "Read in from {}  to output {} ".format(args,tag)
    append = False 
    bb_f.ptclC.write_xmol("{}.xyz".format(tag),comment,append)
    comment = "Read in from {}  to output {} ".format(args,tag)
    append = True
    for pid, ptclObj  in bb_f.ptclC :
         ptclObj.type = ptclObj.tagsDict["fftype"]
    bb_f.ptclC.write_xmol("{}.xyz".format(tag),comment,append)

    bb_f.write_cply("{}.cply".format(tag),write_ff=True,write_bonds=True)

    print tag
    
    # matinfo.write()
    
if __name__=="__main__":
    main()

    
