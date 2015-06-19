#! /usr/bin/env python
"""
Run Molecular Dynamics simulation at low temperature NVT to get a minimized structure
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/02/2014
# travis.kemper@nrel.gov

from structureContainer import StructureContainer
from parameters import ParameterContainer

import pbcs,topology
import mpiBase, file_io

import numpy as np 
import sys , datetime, random, math 

#const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

def get_options():
    """
    Set options
    """
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("-p","--ptime", dest="ptime", default=False,action="store_true", help="Print performance information  ")
    parser.add_option("-o","--output_id", dest="output_id", default="replicate",type="string",help=" prefix for output files  ")
    parser.add_option("-d","--dir_id", dest="dir_id", default=".",type="string",help=" path of output files   ")
    parser.add_option("--calc_overlap", dest="calc_overlap", default=1,type="int",help=" Turn on or off calculation of molecular overlap, handy to build a quick input file" )
    parser.add_option("--calc_type", dest="calc_type", help="Type of calculations output: gromacs or lammps  ", default="gromacs")

    # SWS: adding option so random stream is fixed for testing
    parser.add_option("--fixed_rnd_seed", dest="fixed_rnd_seed", action="store_true", default=False, help=" Fixes random number stream in replicate method in structureContainer (set to True for regression testing")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on ")
    parser.add_option("--sol_json", dest="sol_json", default="",type="string",help=" json files of solvents to add to act on ")
    # Gromacs files to act on 
    parser.add_option("--top", dest="top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--gro", dest="gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--sol_top", dest="sol_top", type="string", default="", help="Input gromacs topology file of solvents (.top) ")
    parser.add_option("--sol_gro", dest="sol_gro", type="string", default="", help="Input gromacs structure file of solvents (.gro) ")

    # Solv
    parser.add_option("--sol_buf", dest="sol_buf",  type=float, default=3.0, help=" intra solvent buffer " )
    # Cut-off
    #   should be replaced by Van der Waals radi  
    parser.add_option("--atomic_cut", dest="atomic_cut", type=float, default=7.5, help="Minimum distance between atoms of molecules ")

    # Replication 
    parser.add_option("--den_target", dest="den_target", type=float, default=0.01, help="Target density g/cm^3 ")
    parser.add_option("--atoms_target", dest="atoms_target", type=int, default=100000, help="Target number of atoms ")
    parser.add_option("--max_mol_place", dest="max_mol_place", type=float, default=50, help="Maximum attempts to place a molecule  ")
    parser.add_option("--max_sys", dest="max_sys", type=float, default=3, help="Maximum system recreations at a certain lattice constant ")
    parser.add_option("--lc_expand", dest="lc_expand", type=float, default=0.100, help="Fraction of the box size to increase system size after max_sys is excieded ")
    parser.add_option("--perc_sol", dest="perc_sol", type=float, default=0.0, help="Percent solvent by mass ")

    # Force filed stuff should be read in from json 
    parser.add_option("--itp_file", dest="itp_file", help=" itp file for all force field parameters  ", default="oplsaa_biaryl.itp")
    parser.add_option("--norm_dihparam", dest="norm_dihparam",default=False,action="store_true",help="Normalize dihedral potential terms if single dihedral is specified in itp file  ")

    (options, args) = parser.parse_args()
        
    return options, args


def main():
    """
    Replicate structures  

    Returns: None
    
    """
    import pbcs, gromacs, lammps , xmol

    options, args = get_options()

    
    # Initialize mpi
    p = mpiBase.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    # Read in oligomers
    #   from json files
    oligo_array = []
    #oligo_param_array = []
    paramC = ParameterContainer()

    strC_i = StructureContainer()  
    #    from gromacs  files
    if( len(options.gro) > 0 ):
        oligo_array,paramC = file_io.struc_array_gromacs(oligo_array,options.gro,options.top,paramC)
        for strC_j in oligo_array:
            strC_i += strC_j
    # Read in solvents
    #   from json files 
    sol_array = []
    sol_param_array = []
    #   from gromacs  files
    if( len(options.gro) > 0 ):
        sol_array,paramC = file_io.struc_array_gromacs(sol_array,options.sol_gro,options.sol_top,paramC)
        for strC_j in sol_array:
            strC_i += strC_j

    
    #
    # Check input
    #
    if( options.calc_type == "lammps" ):
        error_line = " data file output not working yet "
        sys.exit(error_line)

    # Check to make sure all the interactions have parameters  
    norm_dihparam = False
    # Don't need since not using norm_dihparam
    cov_nblist = []
    cov_nbindx = []
    paramC_f,strC_i  = topology.set_param(strC_i,paramC,norm_dihparam,cov_nblist, cov_nbindx)

    f_rep = StructureContainer()    
    f_new = pbcs.replicate(p,options,oligo_array,sol_array)

    # Label all the interactions according to their interaction type
    #   for lammps input
    #   this should be done pre replication then passed to each copy
    norm_dihparam = False 
    paramC_f,f_new  = topology.set_param(f_new,paramC_f,norm_dihparam,cov_nblist, cov_nbindx)


    #print " Sorted parameters "
    #print str(paramC_f)
    #sys.exit("debug 1 ")
    #f_rep.compressPtclIDs()
    #print " f_new prop "
    #print f_new.ptclC
    #print f_new.bondC
    
    if( rank == 0 ):    
        #  Write xmol file 
        xmol_file = options.dir_id +"/" + options.output_id + ".xmol"
        n_part = f_new.getPtclNum()
        n_chains = f_new.getchainnumb()
        comment = " structure with %d particles and %d chains "%(n_part,n_chains)
        append = False
        xmol.write(f_new.ptclC,xmol_file,comment,append)
        #f_new.write_xmol(xmol_file,comment,append)

        # Write json file
	#file_io.write_json(f_new,paramC,options.dir_id,options.output_id )

        # Write gro file 
        #if( options.calc_type == "gromacs" ):
        path_gro_file = options.dir_id +"/" + options.output_id + ".gro"
        if( options.verbose ): print  "     - Writing  ",path_gro_file
        gromacs.print_gro(f_new,path_gro_file)
            
        
        #  Write Lammps input file
        #if( options.calc_type == "lammps" ):
        path_data_file = options.dir_id +"/" + options.output_id + ".data"
        if( options.verbose ): print  "     - Writing  ",path_data_file
        lammps.write_data(f_new,paramC_f,path_data_file)

if __name__=="__main__":
    main()
