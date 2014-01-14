#! /usr/bin/env python
# Convert gromacs structure file into NWCHEM files for VAB calculation

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

    parser.add_option("--debug", dest="debug", default=False, help=" Debug only produce sim.ref ")
        
    (options, args) = parser.parse_args()

    return options, args

def main():
    import os, sys, numpy 
    import gromacs,elements , top    , groups , prop, xmol , nwchem 

    options, args = get_options()
    
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
        group_index_i,group_list_i,group_numb,group_cnt = groups.tempo(  ATYPE,ELN,NBLIST,NBINDEX, options )

    # Find center of mass of groups
    if( len(group_index_i) > 0 ):
        group_cent = groups.cent_mass( group_index_i,group_list_i, R, AMASS,options )
        
    # Create neighbor list 
    g_nb_ind,g_nb_list,g_r_list = groups.nablist(group_cnt,group_cent,LV,options)

    #Initialize reference file
    f_ref = file(options.ref_file, "w")
    f_ref.write( "\n #  Groups:" )
    f_ref.write( "\n #     group ; index number " )
    f_ref.write( "\n #     g_qm_calc ; number of groups ; file id ; index number " )
    f_ref.write( "\n #     gset_qm_calc ; number of groups ; file id ; inter group center of mass ('$\AA$') ; index number " )
    f_ref.close()
    
    f_ref = file(options.ref_file, "a")
    f_ref.write( "#  group_ij ; grp_i ; grp_j ; r_ij ")
    f_ref.write( "\n#  calc id ; grp_i ; grp_j ; q_i ; q_j ")
    f_ref.close
    
    for g_i in range(group_cnt):
        # Write reference file
        f_ref = file(options.ref_file, "a")
        f_ref.write( " \n group %d " % g_i )
	f_ref.close()

        N_o = group_index_i[g_i]
        N_f = group_index_i[g_i + 1 ] - 1
        NA_g = N_f - N_o + 1
        if( options.verbose ): print "        Printing gaussian input for group_i ",g_i," with atoms ",NA_g
        F_id = "g_" + str(g_i) 

        # Record id
        group_id = F_id
        calc_id = F_id


        # Cent group i at origin and apply pbc's so pairs are whole and centered
        r_cm = []
	for d in range(options.prop_dim):
            r_cm.append( -1.0*group_cent[g_i][d] )
        # Shift all atoms to center group i
        
        
        R_cent = prop.shift_r(r_cm,R,options)
        R_pbc = prop.atom_pbc(R_cent,LV,options)
        
        print " shifted to center ",r_cm
    
        R_i = []
        ASYMB_i = []
        ATYPE_i = []
        CHARGES_i = []
        REFNUMB_i = []
	FIX_i = []
	GHOST_i = []
        
        # Place group into local array
        ELECTRONS_i = 0 
        for a_indx in range( N_o,N_f+1):
            i = group_list_i[a_indx]
            if( int(ELN[i]) > 0 ):

                ASYMB_i.append( ASYMB[i] )
                R_i.append( R_pbc[i] )
                ATYPE_i.append( ATYPE[i] )
                CHARGES_i.append( CHARGES[i] )
                ELECTRONS_i += int( ELN[i] )
                REFNUMB_i.append( i )
                FIX_i.append( 0 )
                GHOST_i.append( 0 )

        # H-term
        ASYMB_hterm , R_hterm , ATYPE_hterm , CHARGES_hterm, REFNUMB_hterm,ELECTRONS_hterm ,FIX_hterm ,GHOST_hterm ,FIX_i = prop.hterm_seg(ELN,ASYMB,R_pbc,ATYPE,CHARGES,NBLIST,NBINDEX,ASYMB_i,REFNUMB_i,FIX_i,calc_id,options)

	ASYMB_ih = ASYMB_i + ASYMB_hterm
	R_ih = R_i + R_hterm
	ATYPE_ih = ATYPE_i + ATYPE_hterm
	CHARGES_ih = CHARGES_i + CHARGES_hterm
	REFNUMB_ih = REFNUMB_i + REFNUMB_hterm
	FIX_ih = FIX_i + FIX_hterm
	GHOST_ih = GHOST_i + GHOST_hterm

	g_o = g_nb_ind[g_i]
	g_f = g_nb_ind[g_i+1]
	#n_id_i = n_calcid[g_i]
	#c_id_i = c_calcid[g_i]
	
	#print " group i has nb grps ",g_o,g_f
	
	for j_indx in range(g_o,g_f):
	    g_j = g_nb_list[ j_indx]
            
            # to prevent double counting with symmetric neutral geometries
            if( g_j > g_i ):
                    
                if( options.verbose ):
                    print "          group_j ",g_j
                mag_dr = g_r_list[ j_indx]
                #n_id_j = n_calcid[g_j]
                #c_id_j = c_calcid[g_j]
    
                calc_id_j = "g_" + str(g_j) 
    
    
                #print " group i  j ",g_i,g_j,mag_dr
    
                Nj_o = group_index_i[g_j]
                Nj_f = group_index_i[g_j + 1 ] - 1
            
                R_j = []
                ASYMB_j = []
                ATYPE_j = []
                CHARGES_j = []
                REFNUMB_j = []
                FIX_j = []
                GHOST_j = []
    
                # Place group into local array
                ELECTRONS_j = 0 
                for a_indx in range( Nj_o,Nj_f+1):
                    i = group_list_i[a_indx]
                    if( int(ELN[i]) > 0 ):
    
                        ASYMB_j.append( ASYMB[i] )
                        R_j.append( R_pbc[i] )
                        ATYPE_j.append( ATYPE[i] )
                        CHARGES_j.append( CHARGES[i] )
                        ELECTRONS_j += int( ELN[i] )
                        REFNUMB_j.append( i )
                        FIX_j.append( 0 )
                        GHOST_j.append( 0 )
    
    
                # H-term
                #ASYMB_hterm,R_hterm,ATYPE_hterm,CHARGES_hterm ,FIX_hterm,GHOST_hterm = structure.hterm_seg(ELN,ASYMB,R,ATYPE,CHARGES,REFNUMB_j,ASYMB_j,R_j,ATYPE_j,CHARGES_j,FIX_j,GHOST_j,NBLIST,NBINDEX,options)
    
                #ASYMB_hterm , R_hterm , ATYPE_hterm , CHARGES_hterm, REFNUMB_hterm,ELECTRONS_hterm ,FIX_hterm ,GHOST_hterm ,FIX_j = structure.hterm_seg2(ELN,ASYMB,R,ATYPE,CHARGES,REFNUMB,NBLIST,NBINDEX,ASYMB_j,REFNUMB_j,FIX_j,calc_id,options)
    
                
                ASYMB_hterm , R_hterm , ATYPE_hterm , CHARGES_hterm, REFNUMB_hterm,ELECTRONS_hterm ,FIX_hterm ,GHOST_hterm ,FIX_j = prop.hterm_seg(ELN,ASYMB,R_pbc,ATYPE,CHARGES,NBLIST,NBINDEX,ASYMB_j,REFNUMB_j,FIX_j,calc_id_j,options)
    
    
                ASYMB_jh = ASYMB_j + ASYMB_hterm
                R_jh = R_j + R_hterm
                ATYPE_jh = ATYPE_j + ATYPE_hterm
                CHARGES_jh = CHARGES_j + CHARGES_hterm
                REFNUMB_jh = REFNUMB_j + REFNUMB_hterm
                FIX_jh = FIX_j + FIX_hterm
                GHOST_jh = GHOST_j + GHOST_hterm
                
                    
                # Aproximate cation geom with neutral geom
                ASYMB_n_i = ASYMB_ih
                R_n_i = R_ih
                ASYMB_c_i = ASYMB_ih
                R_c_i = R_ih
                
                ASYMB_n_j = ASYMB_jh
                R_n_j = R_jh
                ASYMB_c_j = ASYMB_jh
                R_c_j = R_jh
                
    
                #nw_log = n_id_i+"/"+ n_id_i+ ".log"
                #ASYMB_n_i, R_n_i = nwchem.read_log(nw_log)
                #nw_log = c_id_i+"/"+ c_id_i+ ".log"
                #ASYMB_c_i, R_c_i = nwchem.read_log(nw_log)
                #nw_log = n_id_j+"/"+ n_id_j+ ".log"
                #ASYMB_n_j, R_n_j = nwchem.read_log(nw_log)
                #nw_log = c_id_j+"/"+ c_id_j+ ".log"
                #ASYMB_c_j, R_c_j = nwchem.read_log(nw_log)
    
                # Calculate distance
                cm_n_i = prop.sys_cent_mass(ASYMB_n_i, R_n_i , options)
                cm_c_i = prop.sys_cent_mass(ASYMB_c_i, R_c_i , options)
                cm_n_j = prop.sys_cent_mass(ASYMB_n_j, R_n_j , options)
                cm_c_j = prop.sys_cent_mass(ASYMB_c_j, R_c_j , options)
                
                # Print xyz
                id_name = "p_"+str(g_i)+"n_"+str(g_j)+"n" 
                out_xyz = id_name +".xyz"
                ASYMB_ij = ASYMB_n_i  + ASYMB_c_j
                R_ij = R_n_i  + R_c_j
                ASYMB_ij.append( "Ar" )
                R_ij.append(cm_n_i)
                ASYMB_ij.append( "Kr" )
                R_ij.append(cm_c_j)
                xmol.write_xyz(ASYMB_ij,R_ij,out_xyz)
                #
                ### Print xyz
                #id_name = "p_"+str(g_i)+"n_"+str(g_j)+"c" 
                #out_xyz = id_name +".xyz"
                #ASYMB_ij = ASYMB_n_i  + ASYMB_c_j
                #R_ij = R_n_i  + R_c_j
                #ASYMB_ij.append( "Ar" )
                #R_ij.append(cm_n_i)
                #ASYMB_ij.append( "Kr" )
                #R_ij.append(cm_c_j)
                #xmol.write_xyz(ASYMB_ij,R_ij,out_xyz)
                #        
                #id_name = "p_"+str(g_j)+"n_"+str(g_i)+"c" 
                #out_xyz = id_name +".xyz"
                #ASYMB_ij = ASYMB_c_i  + ASYMB_n_j
                #R_ij = R_c_i  + R_n_j
                #ASYMB_ij.append( "Ar" )
                #R_ij.append(cm_c_i)
                #ASYMB_ij.append( "Kr" )
                #R_ij.append(cm_n_j)
                #xmol.write_xyz(ASYMB_ij,R_ij,out_xyz)
                #
                f_n_i = prop.r_frac(cm_n_i,LV)
                f_c_i = prop.r_frac(cm_c_i,LV)
                f_n_j = prop.r_frac(cm_n_j,LV)
                f_c_j = prop.r_frac(cm_c_j,LV)
                
                    
                #f_ref = file(options.ref_file, "a")
                
                # Charge state of optimized geom
                Q_i = 0
                Q_j = 1
                
                f_ij = f_n_i - f_c_j
                df_pbc = prop.f_pbc(f_ij,options)
                dr_pbc = prop.f_real(df_pbc,LV)
                mag_dr =  numpy.linalg.norm(dr_pbc)
                
                f_ref = file(options.ref_file, "a")
                f_ref.write( " \n group_ij  %d %d %f " % (g_i,g_j,mag_dr )  )
    
                id_name = "et_"+str(g_i)+"n_"+str(g_j)+"n" 
                out_et = id_name +".nw"
                calc_id = "grp_et"
                calc_info = str(g_i) + " " + str(g_j) + " " + str(Q_i)+" "+str(Q_j) +" "+ str(mag_dr)
                f_ref.write( " \n %s %s %s " % (calc_id,calc_info,id_name )  )            
                
                nwchem.write_et( ASYMB_n_i, R_n_i,Q_i,ASYMB_c_j, R_c_j,Q_j,out_et )
                ASYMB_ij = ASYMB_n_i + ASYMB_c_j
                #
                ## Charge state of optimized geom
                #Q_i = 0
                #Q_j = 1
                #
                #f_ij = f_n_i - f_c_j
                #df_pbc = structure.f_pbc(f_ij,options)
                #dr_pbc = structure.f_real(df_pbc,LV)
                #mag_dr =  numpy.linalg.norm(dr_pbc)
                #
                #f_ref = file(options.ref_file, "a")
                #f_ref.write( " \n group_ij  %d %d %f " % (g_i,g_j,mag_dr )  )
                #
                #id_name = "et_"+str(g_i)+"n_"+str(g_j)+"c" 
                #out_et = id_name +".nw"
                #calc_id = "grp_et"
                #calc_info = str(g_i) + " " + str(g_j) + " " + str(Q_i)+" "+str(Q_j) +" "+ str(mag_dr)
                #f_ref.write( " \n %s %s %s " % (calc_id,calc_info,id_name )  )            
                #
                #nwchem.write_et( ASYMB_n_i, R_n_i,Q_i,ASYMB_c_j, R_c_j,Q_j,out_et )
                #ASYMB_ij = ASYMB_n_i + ASYMB_c_j
                #
                #
                ## Charge state of optimized geom
                #Q_i = 1
                #Q_j = 0
                #
                #f_ij = f_c_i - f_n_j
                #df_pbc = structure.f_pbc(f_ij,options)
                #dr_pbc = structure.f_real(df_pbc,LV)
                #mag_dr =  numpy.linalg.norm(dr_pbc)
                #
                #
                #id_name = "et_"+str(g_j)+"n_"+str(g_i)+"c" 
                #out_et = id_name +".nw"
                #calc_id = "grp_et"
                #calc_info = str(g_i) + " " + str(g_j) + " " + str(Q_i)+" "+str(Q_j) +" "+ str(mag_dr)
                #f_ref.write( " \n %s %s %s " % (calc_id,calc_info,id_name )  )
                #
                #nwchem.write_et( ASYMB_c_i, R_c_i,Q_i,ASYMB_n_j, R_n_i,Q_j,out_et )
                #
                f_ref.close()
	    
	    
if __name__=="__main__":
    main()
