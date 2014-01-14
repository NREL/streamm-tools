#! /usr/bin/env python
# Read in group ref file and pull data from qm calc

# Dr. Travis Kemper
# NREL
# Initial Date 12/18/2013
# Email travis.kemper@nrel.gov
# Version 2.00 

def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] "
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False, help="Verbose output ")

    parser.add_option("--nproc", dest="nproc", type=int, default=24, help=" Number of processors for mpi run ")

    parser.add_option("--user", dest="user", type="string",  default="tkemper", help=" User name  ")
    parser.add_option("--pbcs",  dest="pbcs", default= True , help=" Use periodic boundry conditions ")


    parser.add_option("--cluster_host", dest="cluster_host",type="string", default="peregrine",help=" name of cluster ")

    # System properties
    parser.add_option("--prop_dim", dest="prop_dim", type="int",default="3", help="Spacial dimension of the system")

    # Input output files     
    parser.add_option("-r","--ref_file", dest="ref_file", type="string", default="sim.ref", help=" Reference file used for supesquent reading of data ")
 
    # Gromacs
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file ")
 
    # FF settings
    parser.add_option("--ff_software", dest="ff_software", type="string", default="gromacs", help=" Software package to use for Force Field calculations  ")
    parser.add_option("--atom_types", dest="atom_types", type="string", default="", help="Read atom types that will replace default elements ")
    parser.add_option("--nlist_bonds", dest="nlist_bonds", default=True, help="Build neighbor list from bonds")
       
    # Gromacs
    parser.add_option("--gromacs_sufix", dest="gromacs_sufix", type="string", default="", help=" Sufix for gromacs calculations such as _mpi ")
    parser.add_option("--gromacs_dir", dest="gromacs_dir", type="string", default="", help=" gromacs dir ")
    parser.add_option("--g_center", dest="g_center",  default=False, help=" center minimization ")
    parser.add_option("--load_gromacs", dest="load_gromacs", type="string",  default="", help=" module comand to load gromacs ")

    # Groups
    parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
    parser.add_option("--rm_VS", dest="rm_VS", default=False, help=" remove virtual sites ")
    parser.add_option("--groups_com_cut", dest="groups_com_cut", type="float", default="10.0", help="Cut off for neighbors in QM output ")
        
    (options, args) = parser.parse_args()

    return options, args

def main():
    import os, sys, numpy 
    import gromacs, xmol , nwchem 

    options, args = get_options()
    
    if( options.verbose ):
        print " Collecting qm data speciiefied in ",args
        

    # If input file specified, read it in a reset options accordingly
    Lines = []
    # Read in lines of the input files 
    for input_file in args:
        F = open(input_file,'r' )
        Lines_i = F.readlines()
        F.close()
        Lines = Lines + Lines_i    
    
    # Get lattice vectors from .gro file 
    if( len(options.in_gro) ):
        if( options.verbose ): print "  Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options,options.in_gro)
	
    # Find groups
    group_cnt = 0
    for line in Lines :
        col = line.split()
        if ( len(col) > 1 ):
            if( col[0] == "group_ij" ):
		g_indx = int( col[1])
		if( g_indx > group_cnt ): group_cnt = g_indx
    print "  ",group_cnt," groups found in reference files "

    
    f_et = open("et.dat","w")
    f_et.write("# neutral # ; cation #; N-N (A) ; Vij ")
    
    
    f_err = open("et.err","w")
    f_err.write("# neutral # ; cation #; N-N (A) ; Vij ")
    
    # Record group information
    for line in Lines :
        col = line.split()
        if ( len(col) > 1 ):
            if( col[0] == "grp_et" ):
		g_i = ( int( col[1]) )
		g_j = (  int( col[2]) )
		q_i = (  int( col[3]) ) 
		q_j = (  int( col[4]) )
		mag_dr = (  float( col[5]) )
		g_id = ( col[6] )
		nw_log = g_id +"/" + g_id +".log"
		
		if( nwchem.check_log(nw_log)):
		    
		    geom_name  = "GEOMI"
		    ASYMB_i,R_i = nwchem.read_log(nw_log,geom_name)
		    
		    
		    geom_name  = "GEOMJ"
		    ASYMB_j,R_j = nwchem.read_log(nw_log,geom_name)
		    
		    
		    # Get N position
		    for atom_i in range(len(ASYMB_i)):
			if( ASYMB_i[atom_i].strip() == "N" ):
			    n_i = numpy.array( R_i[atom_i] )
		    
		    
		    # Get N position
		    for atom_j in range(len(ASYMB_j)):
			if( ASYMB_j[atom_j].strip() == "N" ):
			    n_j = numpy.array( R_j[atom_j] )
			    
		    dr_nn = n_j - n_i
		    mag_dr_nn = numpy.linalg.norm(dr_nn)
    
		    geom_name  = "GEOMIJ"
		    ASYMB_ij,R_ij = nwchem.read_log(nw_log,geom_name)
		    
		    
		    out_xyz = "p_"+str(g_i)+"_"+str(g_j) +".xyz"
		    xmol.write_xyz(ASYMB_ij,R_ij,out_xyz)
    
		    vij = nwchem.read_etlog(nw_log)
		    
		    if( q_i == 1 ):
			cation = g_i
			neut = g_j
		    else:
			cation = g_j
			neut = g_i
			
		    
		    f_et.write("\n %d %d  %f %f " % (neut,cation,mag_dr_nn,vij))
		    
		else:
		    f_err.write("\n %s " % (nw_log))
		    
		    

                    
		    
		
    f_et.close()
    f_err.close()
    
if __name__=="__main__":
    main()
