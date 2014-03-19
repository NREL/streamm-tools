#! /usr/bin/env python
# Find volume for a set of molecules to give a specified density 

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
    parser.add_option("--atom_types", dest="atom_types", type="string", default="", help="Read atom types that will replace default elements ")    

    (options, args) = parser.parse_args()

    return options, args

  
def main():
    import sys, gromacs, prop, top ,  elements

    options, args = get_options()

    gro_file = sys.argv[1]
    top_infile = sys.argv[2]
     
    GTYPE, R, VEL, LV = gromacs.read_gro(options,gro_file)
    ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST= gromacs.read_top(options,top_infile)
    ASYMB , ELN  = elements.mass_asymb(AMASS) 
    
    
    # Retype special atom types and replace element
    if( len( options.atom_types)): 
        if( options.verbose ): print "  Reading in ",options.atom_types
    
        # Zero virtual site mass
        ASYMB , ELN  = top.special_types(ATYPE,ASYMB , ELN , options.atom_types)
        for atom_i in range(len(ELN)):
            if( ELN[atom_i] ==  0 ):
                AMASS[atom_i] =  0.0
                        
    
    DEN = prop.vecdensity(  AMASS,LV )
    
    print "   Density ",DEN
    #
    #mol_mass_amu = prop.total_mass( AMASS )
    #
    #
    #mass_amu = mol_mass_amu*mol_mult
    #
    #volume_target_nm = mass_amu/target_density_amunm
    #len_target_nm = volume_target_nm**(1.0/3.0)
    #
    #print mass_amu ,' AMU ',volume_target_nm,' volume nm^3',len_target_nm,' nm cubed '
    #

if __name__=="__main__":
    main()
