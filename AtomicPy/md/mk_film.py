#! /usr/bin/env python
# Subroutines for calculating structural properties

# Dr. Travis Kemper
# Gatech
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

    parser.add_option("--atomic_cut", dest="atomic_cut", type=float, default=2.5, help="Minimum distance between atoms of molecules ")

    parser.add_option("--den_target", dest="den_target", type=float, default=1.0, help="Target density g/cm^3 ")
    parser.add_option("--atoms_target", dest="atoms_target", type=float, default=100000.0, help="Target number of atoms ")
    parser.add_option("--max_mol_place", dest="max_mol_place", type=float, default=10, help="Maximum attempts to place a molecule  ")
    parser.add_option("--max_sys", dest="max_sys", type=float, default=10, help="Maximum system recreations at a certian lattice constant ")
    parser.add_option("--lc_expand", dest="lc_expand", type=float, default=2.5, help="Distance (angstroms) to increase system size after max_sys is excieded ")
    # parser.add_option("--nb_list",dest="nb_list", default=False,action="store_true", help="Use neighbor list ")

    parser.add_option("--out_gro", dest="out_gro", type="string", default="out.gro", help="gromacs output file ")
    
    (options, args) = parser.parse_args()
        
    return options, args
   
   
def main():
    #
    # Caclulate rdf
    #
    
    import os, sys, numpy , math 
    import datetime
    import gromacs, elements, xmol, prop, file_io, groups  #, vectors 
    import random
    
    debug = 0 
    
    options, args = get_options()
    
    
    p_time = 0
    
    if( options.verbose ): print "   Reading in files to establish intial conditions and atom types "
    #
    # Read in top file
    #
    if( len(options.in_top) ):
        if( options.verbose ): print "      Reading in ",options.in_top
        ATYPE_i,RESN_i,RESID_i,GTYPE_i,CHARN_i,CHARGES_i,AMASS_i,BONDS_i,ANGLES_i,DIH_i,MOLNUMB_i,MOLPNT_i,MOLLIST_i = gromacs.read_top(options,options.in_top)
        ASYMB_i,ELN_i  = elements.mass_asymb(AMASS_i)
    #
    # Get coord
    #
    if( len(options.in_gro) ):
        if( options.verbose ): print "     Reading in ",options.in_gro
        GTYPE_i,R_i,VEL_i,LV_i = gromacs.read_gro(options,options.in_gro)

    NA = len(GTYPE_i)
    #
    # Calculate the target density
    #
    mol_mult = int( options.atoms_target/float(NA) )
    target_density_amuang = options.den_target*const_avo/10.0
    mol_mass_amu = prop.total_mass( AMASS_i )
    mass_amu = mol_mass_amu*mol_mult
    
    volume_target_ang = mass_amu/target_density_amuang
    len_target_ang = volume_target_ang**(1.0/3.0)
    
    LV = numpy.zeros([3,3])
    
    LV[0,0] = len_target_ang
    LV[1,1] = len_target_ang
    LV[2,2] = len_target_ang

    #
    # Find target volume for single molecule
    #
    mol_vol = float(mol_mult)/volume_target_ang
    mol_unit_l =  mol_vol**(-1.0/3.0)
        
    if( options.verbose ):
        print "  Input molecule has ",NA," atoms and mass of ",mol_mass_amu," AMU "
        print "  To achieve ",options.atoms_target," atom systems it will be multiplied ",mol_mult
        print "  giving ",mol_vol," mol/Angstrom^3 and a molecular unit length of ",mol_unit_l," Angstorms "
        print "  For the target density of ",options.den_target," g/cnm^3 ",target_density_amuang," AMU Angstrom^-3"
        print "    a target cubic unit cell of ",len_target_ang," will be needed "
        print "Options "
        print "  in_top",options.in_top
        print "  in_gro",options.in_gro
        print "  atomic_cut",options.atomic_cut
        print "  den_target",options.den_target
        print "  atoms_target",options.atoms_target
        print "  max_mol_place",options.max_mol_place
        print "  max_sys",options.max_sys
        print "  lc_expand",options.lc_expand
        print "  out_gro",options.out_gro
    
    #
    # Shift molecule to have center of mass at origin 
    # 
    r_shift = numpy.array( [0.0,0.0,0.0] )
    R_moli_c = prop.shift_cent_mass(AMASS_i,R_i,r_shift)
  
  
    out_xyz ="R_cent.xyz"
    file_xmol ="R_cent.xmol"
    xmol.write_xyz(ASYMB_i,R_moli_c,out_xyz)
    
    lv_i = len_target_ang
    
    # sys.exit(" debug 2  ")
    
    add_mol = 1
    prop_dim = len(R_i[0] )
    
    ang_acc = 1000
    
    #
    # Initialize system
    #
    ASYMB_sys = []
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
    # LV_sys = []
    
    sys_mol_n = 0
    sys_attempts = 0
    
    cut_ij_sq = options.atomic_cut* options.atomic_cut
    
    while ( add_mol ):
        add_mol = 1
        overlap = 1
        moladd_atempts = 0
        
        while ( overlap ):
            moladd_atempts += 1    
            overlap = 0
            
            #
            # Rotate molecule randomly 
            #
            rot_angle_i = float(random.randrange(0,ang_acc))*numpy.pi/float(ang_acc)
            rot_angle_j = float(random.randrange(0,ang_acc))*numpy.pi/float(ang_acc)
            
            if( debug ): print " rotating molecule by ",rot_angle_i,rot_angle_j
            
            cy = math.cos(rot_angle_i)
            sy = math.sin(rot_angle_i)
            cz = math.cos(rot_angle_j)
            sz = math.sin(rot_angle_j)
            
            R_rot = []
            for atom_i in range( len(ELN_i) ):
                xd = R_moli_c[atom_i][0]
                yd = R_moli_c[atom_i][1]
                zd = R_moli_c[atom_i][2]
                r_x =  cy*cz*xd - sz*cy*yd + sy*zd 
                r_y =  sz*xd    + cz*yd            
                r_z = -sy*cz*xd + sy*sz*yd + cy*zd 
                r_i =  numpy.array( [r_x,r_y,r_z] )
                R_rot.append(  r_i )
            
            #
            # Shift molecule to random point 
            #
            mol_origin = numpy.zeros(prop_dim)
            for d in range(prop_dim):
                mol_origin[d] = random.randrange(0,int(lv_i))
                
            R_shift = []
            
            if( debug ): print "  Shifting molecule to random point ",mol_origin
            
            for atom_i in range( len(ELN_i) ):
                R_shift.append(  R_rot[atom_i] +  mol_origin  )
            
            if( debug ): print  " new center of mass ",prop.cent_mass(AMASS_i,R_shift)," should be mol origin ",mol_origin
            
            
            #
            # Calculate overlap with existing molecules in the system 
            #
            
            if( len(ELN_sys) > 0 ):
                
                for atom_i in range( len(ELN_i) ):
                    r_i = R_shift[atom_i]
                    for sys_atom in range( len(ELN_sys) ):
                        r_j = R_sys[sys_atom]
                        r_ij_sq = prop.sq_drij_c(r_i,r_j,LV)
                        if( r_ij_sq < cut_ij_sq ): overlap = 1
                
            else:
                overlap = 0
                
            if( overlap ==  0 ):
                sys_mol_n += 1
                
                if( debug ):  print " Adding molecule to the system ",sys_mol_n
                
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
                    R_sys.append( R_shift[atom_i])
                    #VEL_sys.append( VEL_i[atom_i])
            
            if( moladd_atempts >= options.max_mol_place ):
                if( options.verbose ):
                    print " Attempts to add molecule ",sys_mol_n," is ",moladd_atempts," which has exceeded max attempts ",options.max_mol_place," system will be reset "
                    
                sys_attempts += 1 
                ASYMB_sys = []
                
                ELN_sys  = []
                ATYPE_sys = []
                RESN_sys = []
                RESID_sys = []
                GTYPE_sys = []
                CHARN_sys = []
                CHARGES_sys = []
                AMASS_sys = []
                R_sys = []
                VEL_sys = []
                
                
                if( sys_attempts >= options.max_sys  ):
                    
                    if( options.verbose ):
                        print ' Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                        print ' Lattice vectors will be expanded by (option lc_expand)',options.lc_expand
                        
                    LV[0,0] = LV[0,0] + options.lc_expand
                    LV[1,1] = LV[1,1] + options.lc_expand
                    LV[2,2] = LV[2,2] + options.lc_expand
                
                
                sys_mol_n = 0
                
        if( sys_mol_n ==  mol_mult ): add_mol = 0
        
    
    gromacs.print_gro(options.out_gro,GTYPE_sys,RESID_sys,RESN_sys,R_sys,LV)

if __name__=="__main__":
    main()
   
