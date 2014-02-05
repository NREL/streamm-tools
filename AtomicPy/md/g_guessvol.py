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

    (options, args) = parser.parse_args()

    return options, args

  
def main():
    import sys, gromacs, prop

    options, args = get_options()

    N_a_mag = 6.022141129 

    # Declare gromacs atom id's array
    GTYPE = []
    # Declare coordinate array
    R = []
    # Declare atom type array
    ATYPE = []
    # Declare residue number
    RESNR = []
    # Declare residue id
    RESIDU = []
    # Declare gromacs id
    GTYPE = []
    # Declare charge group number
    CGNR = []
    # Declare atomic charge
    CHARGE = []
    # Declare atomic mass array
    AMASS = []

    gro_file = sys.argv[1]
    top_infile = sys.argv[2]
    target_density_gcm = float(sys.argv[3])  #'g/cm^3'
    mol_mult = float(sys.argv[4])  # number of molecules 
    target_density_amunm = target_density_gcm*N_a_mag*100.0

    ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST= gromacs.read_top(options,top_infile)
    mol_mass_amu = prop.total_mass( AMASS )
    mass_amu = mol_mass_amu*mol_mult
    
    volume_target_nm = mass_amu/target_density_amunm
    len_target_nm = volume_target_nm**(1.0/3.0)
    
    print mass_amu ,' AMU ',volume_target_nm,' volume nm^3',len_target_nm,' nm cubed '
    

if __name__=="__main__":
    main()
