#! /usr/bin/env python
"""
Run Molecular Dynamics simulation at low temperature NVT to get a minimized structure
"""


# Dr. Travis Kemper
# NREL
# Initial Date 7/02/2014
# travis.kemper@nrel.gov

from structureContainer import StructureContainer
import numpy as np 
import mpiNREL

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
    parser.add_option("-o","--output_id", dest="output_id", default="replicate",type="string",help=" prefix for output files  ")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on ")
    parser.add_option("--sol_json", dest="sol_json", default="",type="string",help=" json files of solvents to add to act on ")

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
    
    (options, args) = parser.parse_args()
        
    return options, args

def struc_array_json(json_list):
    """
    Loop over list of json files and place each structure container into a list

    Arguments
      json_list (list) list of json files
      
    Returns
      struc_array (list) of structure objects
      
    """

    # Read in structure containers from json files 
    struc_array = []
    if( len(json_list)):
        # Loop over json files 
        json_files = json_list.split(',')
        if( len(json_files) > 0 ):
            # Read index files from args
            for json_file in json_files:
                struc_i = StructureContainer()
                struc_i.getsys_json(json_file)
                struc_array.append(struc_i)
                
    return struc_array


def main():
    """
    Replicate structures  

    Returns: None
    """

    # Initialize mpi
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    options, args = get_options()

    # Open log file 
    if( rank == 0  ):
        log_file = options.output_id + ".log"
        log_out = open(log_file,"w")
        
    # Read in structure containers from json files 
    struc_array = struc_array_json(options.json)

    # Read in solvent containers from json files 
    sol_array = struc_array_json(options.sol_json)

    struc_cnt   = 0   # Number of structures to replicate 
    struc_nprt    = 0   # Total number of particles in all the structure 
    struc_mass  = 0.0    # Total mass of particles in all the structure 
    for struc_i in struc_array:
        struc_cnt += 1
        nprt = len( struc_i.ptclC )
        struc_nprt += nprt 
        tot_mass = struc_i.getTotMass()
        struc_mass += tot_mass
        nbonds = len( struc_i.bondC )
        if( options.verbose ):
            log_lines = ""
            log_lines += "  Structure %d \n"%struc_cnt 
            log_lines +=  "    Particles  %d \n"%nprt
            log_lines +=  "    Total mass   %f \n"%tot_mass
            log_lines +=  "    Bonds  %d \n"%nbonds
            print log_lines
            log_out.write(log_lines)
            

    sol_cnt   = 0   # Number of solvent to replicate 
    sol_nprt    = 0   # Total number of solvent in all the solture 
    sol_mass  = 0.0    # Total mass of solvent in all the solture 
    for sol_i in sol_array:
        sol_cnt += 1
        nprt = len( sol_i.ptclC )
        sol_nprt += nprt 
        tot_mass = sol_i.getTotMass()
        sol_mass += tot_mass
        nbonds = len( sol_i.bondC )
        if( options.verbose ):
            log_lines = ""
            log_lines += "  Solvents %d \n"%sol_cnt 
            log_lines +=  "    Particles  %d \n"%nprt
            log_lines +=  "    Total mass   %f \n"%tot_mass
            log_lines +=  "    Bonds  %d \n"%nbonds
            print log_lines
            log_out.write(log_lines)

    # Calculate the number of structures and  solvent molecules
    if( options.perc_sol > 0.0 ):
        # perc_sol = n_solperstruc*sol_mass/( n_solperstruc*sol_mass + struc_mass )
        n_solperstruc =   struc_mass /( sol_mass/perc_sol - sol_mass )
        print " %d sets of solvent per set of structures "%n_solperstruc
    else:
        ntot_sol = 0
        ntot_struc = int(float(options.atoms_target)/float(struc_nprt))
        


    #
    # Calculate the box size for a target density 
    #
    target_density_amuang = options.den_target*const_avo/10.0 # densit in AMU/Angstrom^3
    total_mass = struc_mass*ntot_struc
    volume_target_ang = total_mass/target_density_amuang    
    len_target_ang = volume_target_ang**(1.0/3.0)
    cut_ij_sq = options.atomic_cut* options.atomic_cut

    replic_struc = StructureContainer()  # Output replicated structure 
    # Set lattice vector to new box size
    latvec_list = [np.array([len_target_ang,0.0,0.0]),np.array( [0.0,len_target_ang,0.0]),np.array( [0.0,0.0,len_target_ang]) ]    
    replic_struc.setLatVec(latvec_list)
    
    # Print script information and settings 
    if( rank == 0  ):
        
        log_lines =  "  Replication settings \n"
        log_lines += "   - Molecules \n"
        log_lines += "       Set of structures will replicated %d times \n"%(ntot_struc)
        log_lines += "   - Solvents \n"
        log_lines += "       Set of solvent structures will replicated %d times \n"%(ntot_sol)
        log_lines += "   - Tragets \n"
        log_lines += "       Solvent mass percentage %f \n"%(options.perc_sol)
        log_lines += "       Density %f g/cm^3 %f AMU/Angstrom^3 \n"%(options.den_target,target_density_amuang)
        log_lines += "       Total atoms %d \n"%(options.atoms_target)
        log_lines += "       Giving a cubic cell with length of %f Angstroms \n"%(len_target_ang)
        log_lines += "   - Placement \n"
        log_lines += "       Maximum structure placements %d \n"%(options.max_mol_place)
        log_lines += "       Maximum number of restarts before box is expanded %d \n"%(options.max_sys)
        log_lines += "       Percent of box size to add during expantion %8.2f \n"%(100*options.lc_expand)
        print log_lines
        log_out.write(log_lines)

    p.barrier()
    
    # Record initial time
    if( rank == 0  ): 
	t_i = datetime.datetime.now()

    sys_mol_n = 0    # Number of molecules add to the system
    sys_attempts = 0  # Number of times the system has been reset
    #
    # Start adding molecules to the system
    #
    add_mol = 1
    while ( add_mol and options.mult_i > 0 ):
        #
        # Initialize 
        #
        add_mol = 1
        overlap_sum = 1
        moladd_atempts = 0
        
        # Record intial time for pereformance testing 
        if( rank == 0  ): 
            tadd_i = datetime.datetime.now()

        
        for struc_i in struc_array:

            
            
    if( rank == 0 ):
        log_out.close()
        
if __name__=="__main__":
    main()
   

