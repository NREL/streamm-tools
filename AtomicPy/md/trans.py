#! /usr/bin/env python
"""
Create a supercell of randomly placed and rotated molecules 
"""

# Dr. Travis Kemper
# Initial Date April 2014
# travis.kemper@nrel.gov

# length - angstroms
# mass   - AMU
# volume - angstroms^3


const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--itp", dest="itp_file",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")

    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help=" Output xyz file ")
    parser.add_option("--out_gro", dest="out_gro", type="string", default="", help=" Output gromacs gro file ")
    parser.add_option("--out_data", dest="out_data",type="string",default="",help=" Output Lammps data file ")
    parser.add_option("--out_xmol", dest="out_xmol", type="string", default="", help=" Output xmol file ")
    
    (options, args) = parser.parse_args()
        
    return options, args
   

def main():
    #
    # Caclulate rdf
    #
    
    import os, sys, numpy , math , random
    import datetime
    import time    
    import gromacs, elements, xmol, prop, file_io, groups,lammps , top #, vectors 

    # Load information onto all processors 
    #
    options, args = get_options()
    prop_dim = 3
    

    #
    # Get lammps xyz file 
    #
    if( len(options.in_lammpsxyz) ):
        if( options.verbose ): print  "     - Reading in ",options.in_lammpsxyz

        lammpsxyz_F = open(options.in_lammpsxyz , 'r' )
        lammpsxyz_lines = lammpsxyz_F.readlines()
        lammpsxyz_F.close()        
        
        
    #
    # Read in top file
    #
    if( len(options.in_top) ):
        if( options.verbose ): print  "     - Reading in ",options.in_top
        ATYPE_i,RESN_i,RESID_i,GTYPE_i,CHARN_i,CHARGES_i,AMASS_i,BONDS_i,ANGLES_i,DIH_i,MOLNUMB_i,MOLPNT_i,MOLLIST_i = gromacs.read_top(options,options.in_top)
        ASYMB_i,ELN_i  = elements.mass_asymb(AMASS_i)
    #
    # Get gro file 
    #
    if( len(options.in_gro) ):
        if( options.verbose ): print  "     - Reading in ",options.in_gro
        GTYPE_i,R_i,VEL_i,LV_i = gromacs.read_gro(options,options.in_gro)        
    #
    # Read in parameter file 
    #


    #
    # Read in ff file
    #
    if( len(options.itp_file) ):
        if( options.verbose ): print  "     - Reading in ",options.itp_file    
        FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(options.itp_file)
        if(  options.ff_software == "lammps"  ):
            # Identify total number of atom types for lammps output 
            ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND_i , BTYPE_REF, ANGTYPE_IND_i , ANGTYPE_REF, DTYPE_IND_i , DTYPE_REF = lammps.lmp_types(ELN_i,ATYPE_i,AMASS_i,BONDS_i,ANGLES_i,DIH_i)

            #   Build covalent nieghbor list for bonded information 
            NBLIST, NBINDEX = top.build_covnablist(ELN_i,R_i)
            
            # Check atom types to be sure each atom of the same type has the same number of neighbors 
            ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE_i,NBLIST,NBINDEX)
        
            ATYPE_EP, ATYPE_SIG = top.atom_parameters(options.itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
            BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(options.itp_file,BTYPE_IND_i , BTYPE_REF,FF_BONDTYPES)
            ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(options.itp_file,ANGTYPE_IND_i , ANGTYPE_REF,FF_ANGLETYPES)
            DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(options.itp_file, options.norm_dihparam, DTYPE_IND_i , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )
            
            IMPTYPE_F  = top.imp_parameters(options.itp_file)
            
            
    elif(len(options.in_data) == 0 ):
        if(   len(options.out_data) > 0  ):
            print " An itp file specified with the --itp_file option is needed to create a lammps input file "
            sys.exit(" Read in error")
    #
    # Get lammps data file 
    #
    if( len(options.in_data) ):
        if( options.verbose ): print  "     - Reading in ",options.in_data
        ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,DIH_i,DTYPE_IND_i,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,RESN_i,ATYPE_IND_i,CHARGES_i,R_i , ATYPE_i, BONDS_i ,BTYPE_IND_i, ANGLES_i ,ANGTYPE_IND_i, LV_i = lammps.read_data(options.in_data)
        
        AMASS_i = []
        for atom_i in range(len(ATYPE_IND_i)):
            type_ind  = ATYPE_IND_i[atom_i]
            AMASS_i.append( ATYPE_MASS[type_ind])
        ASYMB_i ,ELN_i = elements.mass_asymb(AMASS_i)
        
        #if(  options.ff_software == "gromacs"  ):
        GTYPE_i = []
        RESID_i = []
        CHARN_i = []
        for i in range( len(ELN_i) ):
            GTYPE_i.append(ASYMB_i[i])
            RESID_i.append("MOL")
            CHARN_i.append(RESN_i[i])
            
        #CHARN_i = top.set_chargegroups(options,verbose,CG_SET,CHARN,ATYPE_i,ASYMB_i,ELN,R,NBLIST,NBINDEX, RING_NUMB,LAT_CONST)
        
    #
    # Test that geometry was read in
    #
    try:
        R_i
    except NameError:
        sys.exit("Geometry read in error ")
    

    
    #
    # Initialize system
    #
    ASYMB_sys = []
    ATYPE_IND_sys = []
    ELN_sys  = []
    ATYPE_sys = []
    RESN_sys = []
    RESID_sys = []
    GTYPE_sys = []
    CHARN_sys = []
    CHARGES_sys = []
    AMASS_sys = []
    GTYPE_sys = []
    R_sys = []
    VEL_sys = []
    BONDS_sys = []
    BTYPE_IND_sys = []
    ANGLES_sys = []
    ANGTYPE_IND_sys = []
    DIH_sys = []
    DTYPE_IND_sys = []
        

    for atom_i in range( len(ELN_i) ):
        ASYMB_sys.append( ASYMB_i[atom_i])
        ELN_sys .append( ELN_i[atom_i])
        ATYPE_sys.append( ATYPE_i[atom_i])
        RESN_sys.append( RESN_i[atom_i])
        RESID_sys.append( RESID_i[atom_i])
        GTYPE_sys.append( GTYPE_i[atom_i])
        CHARN_sys.append( CHARN_i[atom_i])
        CHARGES_sys.append( CHARGES_i[atom_i])
        AMASS_sys.append( AMASS_i[atom_i])
        R_sys.append( R_i[atom_i])
        ATYPE_IND_sys.append( ATYPE_IND_i[atom_i])
        
    
    #
    # Output files 
    #
    if( len(options.out_xyz) ):
        if( options.verbose ): print  "     - Writing  ",options.out_xyz
        xmol.write_xyz(ASYMB_sys,R_sys,options.out_xyz)
    

    if( len(options.out_gro) ):
        if( options.verbose ): print  "     - Writing  ",options.out_gro
        gromacs.print_gro(options.out_gro,GTYPE_sys,RESID_sys,RESN_sys,R_sys,LV)

    if( len(options.out_data) ):
        if( options.verbose ): print  "     - Writing  ",options.out_data
        
        # Find topology for entire system
        #   Repeat molecular topology for all molecules added to the system
        mol_mult = 1
        BONDS_sys,BTYPE_IND_sys = replicate_bonds(mol_mult,ELN_i,BONDS_i,BTYPE_IND_i) 
        ANGLES_sys,ANGTYPE_IND_sys = replicate_angles(mol_mult,ELN_i,ANGLES_i,ANGTYPE_IND_i) 
        DIH_sys,DTYPE_IND_sys = replicate_dih(mol_mult,ELN_i,DIH_i,DTYPE_IND_i) 

        debug = 0
        if( debug ):
            print " BONDS ", len(BONDS_sys)
            print "     ",BONDS_sys[0]
            for bond_indx in range(len(BONDS_sys)):
                print BONDS_sys[bond_indx][0]+1, BONDS_sys[bond_indx][1]+1
                
            print " ANGLES_sys ", len(ANGLES_sys)
            print "     ",ANGLES_sys[0]
            print " DIH_sys ", len(DIH_sys)
            print "     ",DIH_sys[0]
            
            sys.exit("debug 6 ")

        lammps.print_lmp(options.out_data,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
          BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
          ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
          DIH_sys,DTYPE_IND_sys,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
          RESN_sys,ATYPE_IND_sys,CHARGES_sys,R_sys , ATYPE_sys,
          BONDS_sys ,BTYPE_IND_sys, ANGLES_sys ,ANGTYPE_IND_sys, LV)

    if( len(options.out_xmol) ):
        if( options.verbose ): print  "     - Writing  ",options.out_xmol


        str_file = open( options.out_xmol, 'w' )
        
        if( len(options.in_lammpsxyz) ):
            if( options.verbose ): print    "       - Reprocessing  ",options.out_xmol
            if(  len(options.in_data)  ):
                if( options.verbose ): print    "       -  with data from ",options.in_data
            else:
                print " no atype type reference "
                sys.exit(" missing topology information ")
            
            #
            # Find the atomic symbol and element number for each atom type
            #
            ATYPE_SYMB , ATYPE_ELN = elements.mass_asymb(ATYPE_MASS)
            n_frames = int( float( len(lammpsxyz_lines) )/ float( len(ASYMB_sys) + 2) )
            if( options.verbose ):
                print " Translating atom types "
                for type_i in range(len(ATYPE_MASS)):
                    print type_i + 1," -> ",ATYPE_SYMB[type_i] #, ATYPE_ELN[type_i],ATYPE_MASS[type_i]
                print " print ",n_frames," frames "
                
            line_cnt = -1
            
            for frame_i in range(n_frames):
                
                line_cnt += 1
                str_file.write( lammpsxyz_lines[line_cnt] )
                line_cnt += 1
                
                if( options.verbose ):
                    print " reading frame ",frame_i," starting at line ",line_cnt-1," with comment ",lammpsxyz_lines[line_cnt] 
                
                str_file.write( lammpsxyz_lines[line_cnt] )
                
                #sys.exit("debug 4")
                # R_frame_i = []
                for atom_i in range(len(ASYMB_sys) ):   
                    line_cnt += 1
                    
                    if( line_cnt > len(lammpsxyz_lines)-1):
                        print " frame is missing some atoms ",atom_i," not found "
                        # sys.exit("read in e)
                        
                    col =  lammpsxyz_lines[line_cnt].split()
                    if( len(col) >= 4 ):
                        type_i = int(col[0]) - 1
                        r_x = float(col[1])
                        r_y = float(col[2])
                        r_z = float(col[3])
                        str_file.write( "%s %f %f %f \n" % (ATYPE_SYMB[type_i],r_x,r_y,r_z) )
                        
                        # R_frame_i.append( numpy.array( [r_x,r_y,r_z] ) )
                
            # xmol.print_xmol(ASYMB_sys,R,file_xmol)
        str_file.close()
        
            
                    
                    
                    
        
        

if __name__=="__main__":
    main()
   
