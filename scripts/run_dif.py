#!/opt/local/bin/python
"""
Run diffractometer
"""

import sys, datetime
from modeler_hoomd.analyze import diffractometer
from modeler_hoomd.builder import hoomd_xml
import numpy

from structureContainer import StructureContainer
import particles 
import mpiNREL,  file_io 


def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    # json files to act on
    parser.add_option("-j","--in_json", dest="in_json", type="string", default="", help="Input json file, which read in first then over writen by subsequent input files")
    parser.add_option("-o","--output_id", dest="output_id", default="trans",type="string",help=" prefix for output files  ")

    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--itp", dest="itp_file",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_xmol", dest="in_xmol", type="string", default="", help="Input xmol file")
    parser.add_option("--xmol_format", dest="xmol_format", type="string", default="atomic_symb", help="Format of xmol file: atomic_symb - first colum atomic symbols; lammps - first colum particle type atomic symbols will be set to carbon; ")
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")
    #
    # Filters
    #
    parser.add_option("--id", dest="id", type="string", default="", help=" select atoms of group by number  ")    
    parser.add_option("--symbol", dest="symbol", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--type", dest="type", type="string", default="", help=" select type  ")    
    parser.add_option("--chains", dest="chains", type="string", default="", help="select atoms of group by chain number  ")    
    parser.add_option("--ring", dest="ring", type="string", default="", help="select atoms of group by particlesn a ring   ")    
    parser.add_option("--resname", dest="resname", type="string", default="", help="select atoms of group by residue name  ")    
    parser.add_option("--residue", dest="residue", type="string", default="", help="select atoms of group by resudue number  ")    
    parser.add_option("--linkid", dest="linkid", type="string", default="", help="select atoms of group by  linkd ")    
    parser.add_option("--fftype", dest="fftype", type="string", default="", help="select atoms of group by force field type  ")
    parser.add_option("--gtype", dest="gtype", type="string", default="", help="select atoms of group by force field type  ")
    
    #parser.add_option("--filter_lmptype", dest="filter_lmptype", type="string", default="", help=" filter atoms by lammps type ")
 
    (options, args) = parser.parse_args()
        
    return options, args
   

def main():
    """
    Read in files, filter particles and produce an xrd pattern 

    Arguments

    Return
    
    """
    
    #
    # Formated ouput varables
    #
    sperator_line = "\n---------------------------------------------------------------------"
    

    # Initialize mpi
    p = mpiNREL.getMPIObject()

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
    if( options.verbose ):
        print " Reading in files"
    
    struc_o = StructureContainer()
    struc_o = file_io.getstrucC(struc_o, options.in_json, options.in_gro , options.in_top, options.in_data,options.in_xmol,options.xmol_format )

    
    #
    # Get paticle and bond structures
    #
    ptclC_o = struc_o.ptclC
    bondC_o  = struc_o.bondC
    
    #   
    # Filter particles
    #
    if( options.verbose ):
        print " Filter particles"

    search_o = particles.create_search(options.id,options.type,options.symbol,options.chains,options.ring,options.resname,options.residue,options.linkid,options.fftype,options.gtype)
    if( rank == 0 ):
        if( options.verbose ): print " Filter input by ",search_o,len(search_o)
    list_f = ptclC_o.getParticlesWithTags(search_o)
    sum_f = len(list_f)

    #
    # Put coordinates in a numpy array for diffraction
    #
    debug = False 
    coor_list = [] #numpy.array()
    for pid in list_f:
        coor_list.append( ptclC_o[pid].position )
        if( debug ): print "  adding ",pid, ptclC_o[pid].tagsDict["symbol"] , ptclC_o[pid].position

    if( debug): sys.exit(" debug ptcls ")
        
    coor_i = numpy.array( coor_list )
    




        #x = numpy.loadtxt('just_coords')
    #
    # Run diffractometer
    #
    if( options.verbose ):
        print " Running diffractometer"
        
    D = diffractometer.diffractometer()
    D.set(zoom=2, n_angles=10, length_scale=1.0)
    latvec = struc_o.get_latvec()
    box = [ latvec[0][0], latvec[1][1], latvec[2][2] ]
    D.load(coor_i,box)
    D.prep_matrices()
    D.average()


if __name__=="__main__":
    main()
   
