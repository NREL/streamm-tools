#! /usr/bin/env python
"""
Run Molecular Dynamics simulation at low temperature NVT to get a minimized structure
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/02/2014
# travis.kemper@nrel.gov

from structureContainer import StructureContainer
import pbcs
import mpiNREL

import numpy as np 
import sys , datetime, random, math 

const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

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
    parser.add_option("--atomic_cut", dest="atomic_cut", type=float, default=2.5, help="Minimum distance between atoms of molecules ")

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
    import pbcs

    options, args = get_options()

    
    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    f_rep = StructureContainer()
    f_new = f_rep.replicate(p,options)
    #f_rep.compressPtclIDs()

    #print " f_new prop "
    #print f_new.ptclC
    #print f_new.bondC
    
    if( rank == 0 ):    
        #  Write xmol file 
        xmol_file = options.dir_id +"/" + options.output_id + ".xmol"
        n_part = f_new.getpartnumb()
        n_chains = f_new.getchainnumb()
        comment = " structure with %d particles and %d chains "%(n_part,n_chains)
        append = False 
        f_new.write_xmol(xmol_file,comment,append)

        # Write json file
	f_new.write_json(options.dir_id,options.output_id )
        
        #  Write Lammps input file
        #path_data_file = options.dir_id +"/" + options.output_id + ".data"
        #print " Writint file ",path_data_file
        #f_new.lmp_writedata(path_data_file,options.norm_dihparam,options.itp_file)

        # Write gromacs input files python replicate.py   --gro SOL.gro --top SOL.top   --sol_gro SOL.gro   --sol_top SOL.top   --den_target 0.1  --atoms_target 100    --perc_sol 90
        
        f_new.write_gro(options.dir_id,options.output_id )

if __name__=="__main__":
    main()
