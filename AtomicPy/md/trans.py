#! /usr/bin/env python
"""
Translate between various file types

 length - angstroms
 mass   - AMU
 volume - angstroms^3



ELN_i       element number
ASYMB_i
CTYPE_i
CHARGES_i
UNITNUMB_i
UNITTYPE_i
R_i

"""

# Dr. Travis Kemper
# Initial Date April 2014
# travis.kemper@nrel.gov


const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    # json files to act on
    #parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")    
    parser.add_option("-j","--in_json", dest="in_json", type="string", default="", help="Input json file")
    
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file (.top) ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file (.gro) ")
    parser.add_option("--itp", dest="itp_file",  type="string", default="",help="gromacs force field parameter file")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")

    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help=" Output xyz file ")
    parser.add_option("--out_data", dest="out_data",type="string",default="",help=" Output Lammps data file ")
    parser.add_option("--out_xmol", dest="out_xmol", type="string", default="", help=" Output xmol file ")
    #
    #
    #
    parser.add_option("--out_top", dest="out_top", type="string", default="", help=" Output gromacs topology file ")
    parser.add_option("--out_gro", dest="out_gro", type="string", default="", help=" Output gromacs structure file ")
    parser.add_option("--out_itp", dest="out_itp", type="string", default="", help=" Output gromacs parameter file ")
    #
    # 
    #
    parser.add_option("--out_json", dest="out_json", type="string", default="", help=" Output json file ")
    #
    # Filters
    #
    
    parser.add_option("--filter_rings", dest="filter_rings", default=False,action="store_true", help="Only print atoms in rings ")
    #
    # Modify input system file 
    #
    parser.add_option("--replicate_n", dest="replicate_n", default=1,type=int, help="Replicate system n times ")
    
    (options, args) = parser.parse_args()
        
    return options, args
   

def main():
    #
    # Caclulate rdf
    #
    
    import os, sys, numpy , math , random, json
    import datetime
    import time    
    import gromacs, elements, xmol, prop, file_io, groups,lammps , top, jsonapy

    # Load information onto all processors 
    #
    options, args = get_options()
    prop_dim = 3

 
    #
    # Read in json file
    #
    if( len(options.in_json) ):
        if( options.verbose ):
            print  "     - Reading in ",options.in_json
	json_data,json_success = jsonapy.read_jsondata(options.in_json)
	if(  json_success ):
	    #
            mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
            #
            # Need meta data to proceed 
            #      		    
            if( metadata_found ):
                ELN_i,ASYMB_i,CTYPE_i,CHARGES_i,UNITNUMB_i,UNITTYPE_i,R_i,VEL_i,ATYPE_i,AMASS_i,MOLNUMB_i,RING_NUMB_i,RESID_i,RESN_i,CHARN_i,LV_i,json_atomicdata  = jsonapy.read_atomic(json_data)
                BONDS_i  = jsonapy.read_connections(json_data)
                NBLIST_i, NBINDEX_i  = groups.build_nablist_bonds(ELN_i,BONDS_i)
                NA_i = len(ELN_i)
                #BONDS_i = top.nblist_bonds(NA,NBLIST, NBINDEX)
                ANGLES_i = top.nblist_angles(NA_i,NBLIST_i, NBINDEX_i)
                #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)

                # New options that need to be passed 
                limdih =  0
                limitdih_n = 1
                DIH_i = top.nblist_dih(NA_i,NBLIST_i, NBINDEX_i,limdih,limitdih_n)
                IMPS_i = top.nblist_imp(NA_i,NBLIST_i, NBINDEX_i,ELN_i)
			                
                #
                # 
                #
                GTYPE_i = []
                for i in range( len(ELN_i) ):
                    GTYPE_i.append(ASYMB_i[i])

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
        ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,DIH_i,DTYPE_IND_i,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,MOLNUMB_i,ATYPE_IND_i,CHARGES_i,R_i , ATYPE_i, BONDS_i ,BTYPE_IND_i, ANGLES_i ,ANGTYPE_IND_i, LV_i = lammps.read_data(options.in_data)

        
        #
        # Find the atomic symbol and element number for each atom type
        #
        ATYPE_SYMB , ATYPE_ELN = elements.mass_asymb(ATYPE_MASS)
        if( options.verbose ):
            print " Translating atom types "
            for type_i in range(len(ATYPE_MASS)):
                print type_i + 1," -> ",ATYPE_SYMB[type_i] #, ATYPE_ELN[type_i],ATYPE_MASS[type_i]
                
        AMASS_i = []
        for atom_i in range(len(ATYPE_IND_i)):
            type_ind  = ATYPE_IND_i[atom_i]
            AMASS_i.append( ATYPE_MASS[type_ind])
        ASYMB_i ,ELN_i = elements.mass_asymb(AMASS_i)
        
        #if(  options.ff_software == "gromacs"  ):
        GTYPE_i = []
        RESID_i = []
        RESN_i = []
        CHARN_i = []
        VEL_i = []
        CTYPE_i = []
        UNITNUMB_i = []
        UNITTYPE_i = []
        RING_NUMB_i  = []
        for i in range( len(ELN_i) ):
            GTYPE_i.append(ASYMB_i[i])
            RESID_i.append("MOL")
            CHARN_i.append(MOLNUMB_i[i])
            RESN_i.append(MOLNUMB_i[i])
            VEL_i.append( numpy.array( [0.0 ,0.0 ,0.0]) )
            CTYPE_i.append(MOLNUMB_i[i])
            UNITNUMB_i.append(MOLNUMB_i[i])
            UNITTYPE_i.append(MOLNUMB_i[i])
            RING_NUMB_i.append(MOLNUMB_i[i])
            
        #CHARN_i = top.set_chargegroups(options,verbose,CG_SET,CHARN,ATYPE_i,ASYMB_i,ELN,R,NBLIST,NBINDEX, RING_NUMB,LAT_CONST)
        N_MOL_i,MOLPNT_i,MOLLIST_i = groups.molecule_list(MOLNUMB_i)
        
    #
    # Test that geometry was read in
    #
    try:
        R_i
    except NameError:
        sys.exit("Geometry read in error ")
    
    #
    # Calculate properties
    #
    if( options.filter_rings ):
        # NBLIST_i, NBINDEX_i = top.build_covnablist(ELN_i,R_i)
        NBLIST_i, NBINDEX_i  = groups.build_nablist_bonds(ELN_i,BONDS_i)
        RINGLIST_i, RINGINDEX_i , RING_NUMB_i = top.find_rings(ELN_i,NBLIST_i,NBINDEX_i)
        #
        #
		
    #
    # Test that bonds have been set 
    #
    try:
        BONDS_i
    except NameError:
        if( options.verbose):
            print " No bonds found from read in will generated with a covalent nieghbor list"
            
        #   Build covalent nieghbor list for bonded information 
        NBLIST, NBINDEX = top.build_covnablist(ELN_i,R_i)

        NA = len(ELN_i)

        BONDS_i = top.nblist_bonds(NA,NBLIST, NBINDEX)
        ANGLES_i = top.nblist_angles(NA,NBLIST, NBINDEX)
        #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)

        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        DIH_i = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS_i = top.nblist_imp(NA,NBLIST, NBINDEX,ELN_i)
			
    #
    # Test that lammps identifiers have been set 
    #
    try:
        ATYPE_IND_i
    except NameError:
        # if no lammps data read in

        if( options.verbose):
            print " No lammps types found in read in will generate"
                    
        # Identify total number of atom types for lammps output 
        ATYPE_IND_i , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN_i,ATYPE_i,AMASS_i,BONDS_i,ANGLES_i,DIH_i)

    #
    # Print system properties
    #

    
    NA_i = len(ELN_i)
    mass_amu = prop.total_mass( AMASS_i )
    volume_i = prop.volume( LV_i )
    density_i = prop.vecdensity(  AMASS_i,LV_i )
    rank = 0
    if( rank == 0  ):
        print "   - Input system "
        print "       Atoms ",NA_i
        print "       Mass ",mass_amu
        print "       Box size ",LV_i #[0][0]
        print "       Volume ",volume_i
        print "       Density ",density_i
        print "       Molecules ",max(MOLNUMB_i) + 1 
        #print "       Conjugated rings ",max(RING_NUMB_i) + 1 
        #print "       Residues ",max(RESN_i) + 1 
        #print "       Units ",max(UNITNUMB_i) + 1 
        #print "       Charge groups ",max(CHARN_i) + 1 
        print "       "
        print "       "
        print "       "
        print "       "

        print "       BONDS ",len(BONDS_i)
        print "       "

    if( options.replicate_n > 1 ):
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


        CTYPE_sys  = []
        UNITNUMB_sys  = []
        UNITTYPE_sys  = []
        MOLNUMB_sys  = []
        RING_NUMB_sys  = []


        if( options.verbose):
            print " Replicating atomic information ", options.replicate_n," times "

        for unit_n in range( options.replicate_n ):
            # Loop over specifed number of molecules
            #   default is 1

            if( options.verbose):
                print " Replicating unit  ", unit_n," with ",NA_i," atoms "

            for atom_i in range(NA_i):

                #add_atom = 0
                #if( options.filter_rings ):
                #    if( RING_NUMB_i[atom_i] > 0 ):
                #        add_atom = 1
                #else:
                #    add_atom = 1

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
                VEL_sys.append( VEL_i[atom_i])
                ATYPE_IND_sys.append( ATYPE_IND_i[atom_i])
                CTYPE_sys.append( CTYPE_i[atom_i])
                UNITTYPE_sys.append( [atom_i])
                # Shift numbered groups 
                UNITNUMB_sys.append(  UNITNUMB_i[atom_i] + unit_n*max(UNITNUMB_i) )
                MOLNUMB_sys.append( MOLNUMB_i[atom_i] + unit_n*max(MOLNUMB_i)  )
                RING_NUMB_sys.append( RING_NUMB_i[atom_i] + unit_n*max(RING_NUMB_i)  )
    else:
        
        ASYMB_sys = ASYMB_i
        ATYPE_IND_sys = ATYPE_IND_i
        ELN_sys  = ELN_i
        ATYPE_sys = ELN_i
        RESN_sys = RESN_i
        RESID_sys = RESID_i
        GTYPE_sys = GTYPE_i
        CHARN_sys = CHARN_i
        CHARGES_sys = CHARGES_i
        AMASS_sys = AMASS_i
        R_sys = R_i
        VEL_sys = VEL_i


        CTYPE_sys  = CTYPE_i
        UNITNUMB_sys  = UNITNUMB_i
        UNITTYPE_sys  = UNITTYPE_i
        MOLNUMB_sys  = MOLNUMB_i
        RING_NUMB_sys  = RING_NUMB_i


    #
    # System connection information
    #
    if( options.verbose):
        print " Replicating connection information "
                        
    BONDS_sys = []
    BTYPE_IND_sys = []
    ANGLES_sys = []
    ANGTYPE_IND_sys = []
    DIH_sys = []
    DTYPE_IND_sys = []

    LV_sys = LV_i
    BONDS_sys,BTYPE_IND_sys = top.replicate_bonds(options.replicate_n,ELN_i,BONDS_i,BTYPE_IND_i) 
    ANGLES_sys,ANGTYPE_IND_sys = top.replicate_angles(options.replicate_n,ELN_i,ANGLES_i,ANGTYPE_IND_i) 
    DIH_sys,DTYPE_IND_sys = top.replicate_dih(options.replicate_n,ELN_i,DIH_i,DTYPE_IND_i) 


    
    NA_sys = len(ELN_sys)
    mass_amu = prop.total_mass( AMASS_sys )
    volume_sys = prop.volume( LV_sys )
    density_sys = prop.vecdensity(  AMASS_sys,LV_sys )
    rank = 0
    if( rank == 0  ):
        print "   - Output system "
        print "       Atoms ",NA_sys
        print "       Mass ",mass_amu
        print "       Box size ",LV_sys[0][0]
        print "       Volume ",volume_sys
        print "       Density ",density_sys
        print "       Molecules ",max(MOLNUMB_sys)
        print "       Conjugated rings ",max(RING_NUMB_sys)
        print "       Residues ",max(RESN_sys)
        print "       Units ",max(UNITNUMB_sys)
        print "       Charge groups ",max(CHARN_sys)
        print "       "
        print "       "
        print "       "
        print "       "

        print "       BONDS ",len(BONDS_sys)
        print "       "
        
    #
    # Get lammps xyz file 
    #
    if( len(options.in_lammpsxyz) ):
        if( options.verbose ): print  "     - Reading in ",options.in_lammpsxyz

        #if(  len(options.in_data)  ):
            #if( options.verbose ): print    "       -  with data from ",options.in_data
        #else:
            #print " no atype type reference "
            #sys.exit(" missing topology information ")
            

        lammpsxyz_F = open(options.in_lammpsxyz , 'r' )
        lammpsxyz_lines = lammpsxyz_F.readlines()
        lammpsxyz_F.close()        

        if( len(options.out_xmol) ):
            if( options.verbose ): print  "     - Writing  ",options.out_xmol
            str_file = open( options.out_xmol, 'w' )
        
        if( len(options.in_lammpsxyz) ):
            
            line_cnt = -1
            n_frames = int( float( len(lammpsxyz_lines) )/ float( len(ASYMB_sys) + 2) )
            for frame_i in range(n_frames):
                
                line_cnt += 1
                if( len(options.out_xmol) ):
                    str_file.write( lammpsxyz_lines[line_cnt] )
                line_cnt += 1
                
                if( options.verbose ):
                    print " reading frame ",frame_i," starting at line ",line_cnt-1," with comment ",lammpsxyz_lines[line_cnt] 
                
                if( len(options.out_xmol) ):
                    str_file.write( lammpsxyz_lines[line_cnt] )
                
                #sys.exit("debug 4")
                R_sys = []
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
                        if( len(options.out_xmol) ):
                            str_file.write( "%s %f %f %f \n" % (ATYPE_SYMB[type_i],r_x,r_y,r_z) )
                        
                        R_sys.append( numpy.array( [r_x,r_y,r_z] ) )
                
            # xmol.print_xmol(ASYMB_sys,R,file_xmol)
        if( len(options.out_xmol) ):
            str_file.close()
        
            
                    
            
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


    if( len(options.out_json) ):
        if( options.verbose ):
            print  "     - Writing  ",options.out_json
        json_data = jsonapy.write_atomic(json_data,ELN_sys,ASYMB_sys,CTYPE_sys,CHARGES_sys,UNITNUMB_sys,UNITTYPE_sys,R_sys,VEL_sys,ATYPE_sys,AMASS_sys,MOLNUMB_sys,RING_NUMB_sys,RESID_sys,RESN_sys,CHARN_sys,LV_sys)

        json_data = jsonapy.write_connections(json_data,BONDS_sys)
                                         
        f = open(options.out_json, 'w')
        json.dump(json_data, f, indent=2)
        f.close()
	
				
        

if __name__=="__main__":
    main()
   
