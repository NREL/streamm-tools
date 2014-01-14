#! /usr/bin/env python
# Minimize indvidual groups 

# Dr. Travis Kemper
# NREL
# Initial Date 10/28/2013
# travis.kemper@nrel.gov





def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog <file to convert> [options] "
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    parser.add_option("--cluster_host", dest="cluster_host",type="string", default="",help=" name of cluster ")

    parser.add_option("--pbcs",  dest="pbcs", default= True , help=" Use periodic boundry conditions ")
    parser.add_option("--npros", dest="npros", type="int",default="1", help="Number of processors to run calculations")

    # Group 
    parser.add_option("--h_term", dest="h_term", default=True, help=" Hydrogen terminate dangeling bonds of group")

    parser.add_option("--software", dest="software", type="string", default="gromacs", help=" Software package to use for Mininimization calculations  ")
    
    # FF settings 
    parser.add_option("--atom_types", dest="atom_types", type="string", default="", help="Read atom types that will replace default elements ")    
    parser.add_option("--nlist_bonds", dest="nlist_bonds", default=True, help="Build neighbor list from bonds")
    parser.add_option("--find_rings", dest="find_rings", default=False, help="Find rings ")
    #parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
    parser.add_option("--ff_file",dest="ff_file",  type="string", default="ff-new.itp", help=" itp file to read in for ff parameters ")
    parser.add_option("--zero_q", dest="zero_q", default=False, help=" Zero atomic charges ")

    # Gromacs
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file ")
    parser.add_option("--out_gro", dest="out_gro", type="string", default="out.gro", help="Output gromacs structure file ")

    parser.add_option("--gromacs_sufix", dest="gromacs_sufix", type="string", default="", help=" Sufix for gromacs calculations such as _mpi ")
    parser.add_option("--gromacs_dir", dest="gromacs_dir", type="string", default="", help=" gromacs dir ")
    parser.add_option("--g_center", dest="g_center",  default=False, help=" center minimization ")
    parser.add_option("--load_gromacs", dest="load_gromacs", type="string",  default="", help=" module comand to load gromacs ")

    #Lammps
    parser.add_option("--lammp_dir", dest="lammp_dir", type="string",  default="", help=" Lammps executable directory ")
    
    # Groups
    parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
    parser.add_option("--opt_groups", dest="opt_groups", default=True, help=" Optimize groups individually ")
    parser.add_option("--rm_VS", dest="rm_VS", default=False, help=" remove virtual sites ")
    

    # QM Settings
    parser.add_option("--qm_method", dest="qm_method", type="string",default="HF", help="Method of QM calculation ")
    parser.add_option("--qm_basis", dest="qm_basis", type="string",default="6-31G", help="Basis set of QM calculation ")
    parser.add_option("--qm_kywd", dest="qm_kywd", type="string",default="", help="Key words for QM calculation ")
    parser.add_option("--qm_mult", dest="qm_mult", type="int",default="0", help=" Shift in default spin multiplicity ( singlet,doublet) QM calculation, allows for triplets ")
    parser.add_option("--qm_npros", dest="qm_npros", type="int",default="1", help="Number of processors to run secondary calculations")
    parser.add_option("--qm_exclude", dest="qm_exclude", type="string",default="1", help="Exclude ")
    parser.add_option("--qm_charge", type="int",action="append",default="0", help="Input gaussain log file ")
    parser.add_option("--qm_triplet", dest="qm_triplet", default=False,  help="Output triplets ")
    
    (options, args) = parser.parse_args()

    return options, args
    
def main():
    import os, sys
    import gromacs, elements, top, groups

    options, args = get_options()

    # Print options
    #simpy.print_options(options)
    
    if( options.verbose ):
        print "  software ",options.software 
        print "  pbcs ",options.pbcs
        print "  "

    if( len(options.in_gro) ):
        if( options.verbose ): print "  Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options,options.in_gro)
        
    # Read in gro file
    if( len(options.in_top) ):
        if( options.verbose ): print "  Reading in ",options.in_top
        ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST  = gromacs.read_top(options,options.in_top)
        ASYMB , ELN  = elements.mass_asymb(AMASS) 
        
    # Retype special atom types and replace element
    if( len( options.atom_types)): 
        if( options.verbose ): print "  Reading in ",options.atom_types
        ASYMB , ELN  = top.special_types(ATYPE,ASYMB , ELN , options.atom_types)
        
    # Print system information
    if( options.verbose ):
	print " prop "
        #prop.print_prop( AMASS,ELN,LV,CHARGES )
        ##top.print_prop( BONDS,ANGLES,DIH )
        
    # Create neighbor list
    if( options.nlist_bonds ):
        if( options.verbose ): print "  Creating neighbor list with bonds"
        NBLIST,NBINDEX  = groups.build_nablist_bonds(ELN,BONDS)

    #Find rings
    #if( options.find_rings ):
    #   RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)

    # Find groups
    if( options.group_ptma ):
        if( options.verbose ): print "  Find groups of TEMPO "
        group_index_i,group_list_i,group_numb = groups.tempo(  ATYPE,ELN,NBLIST,NBINDEX, options )
    
    if( options.opt_groups ):
        if( options.verbose ): print "  Print groups into structure files"
        R_opt = groups.ff_opt_groups(group_index_i,group_list_i, ELN,ASYMB,R,ATYPE,GTYPE,CHARGES,CHARN,AMASS,RESID,RESN,BONDS,ANGLES,DIH ,LV,NBLIST,NBINDEX,options)
        
    ##os.chdir(work_dir)
    #options.out_gro = "min_groups.gro"
    if( options.out_gro ):
        if( options.verbose ): print "  Writing gro file ",options.out_gro
	gromacs.print_gro(options.out_gro,GTYPE,RESID,RESN,R_opt,LV)
        
        #pbc_whole = "echo -e \" 0 \n \" | trjconv -f " + options.out_gro -s ${RNID}.tpr -trans 1 1 1 -o  ${ID}-MIN.gro -pbc whole
        

if __name__=="__main__":
    main()