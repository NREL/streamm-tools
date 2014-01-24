#! /usr/bin/env python
# Read in group ref file and pull data from qm calc

# Dr. Travis Kemper
# NREL
# Initial Date 12/18/2013
# Email travis.kemper@nrel.gov
# Version 2.00 
#


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] "
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False, help="Verbose output ")

    parser.add_option("--in_et", dest="in_et", type="string", default="et.dat", help="Iinput file of V_AB values  ")

    # Gromacs
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file ")
 
    # Groups
    parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
 
    parser.add_option("--nlist_bonds", dest="nlist_bonds", default=True, help="Build neighbor list from bonds")

    parser.add_option("--atom_types", dest="atom_types", type="string", default="", help="Read atom types that will replace default elements ")

    (options, args) = parser.parse_args()

    return options, args

def ket(vij):
    import math 

    dg=0
    lamb=1.02
    kb=8.6173324*10**-5
    hbar=6.58211928*10**-16
    
    T=298.0

    ket = vij**2/hbar*math.sqrt( math.pi/(kb*T*lamb))*math.e**( -(lamb+dg)**2/(4*lamb*kb*T))
    
    return ket


def main():
    import os, sys, numpy , math 
    import file_io, gromacs , top, elements, groups, prop

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

    et_file = options.in_et
    if( file_io.file_exists(et_file)):
        f_et = open(et_file,"r")
        Lines = f_et.readlines()
        f_et.close()

    bin_size = 0.1 #(angstroms)
    v_min = 0.0
    v_max = 5000.0
    v_floor = int(v_min/bin_size)*bin_size
    val_range = v_max - v_floor
    n_bins = int(val_range/bin_size) + 1

    bins = numpy.zeros(n_bins)
    et_sq_bins  = numpy.zeros(n_bins)
    bin_cnt = numpy.zeros(n_bins)

    et_sq_bins_intra  = numpy.zeros(n_bins)
    bin_cnt_intra = numpy.zeros(n_bins)
    et_sq_bins_inter  = numpy.zeros(n_bins)
    bin_cnt_inter = numpy.zeros(n_bins)
    
    grp_et = numpy.empty((3000,3000),dtype="float")
    
    print len(grp_et)
    print len(grp_et[0])
    print grp_et[0][0]

    volume_i =   prop.volume( LV ) 

    for line in Lines:
        col = line.split()
        if( col[0] != "#" ):
            n_i = int( col[0] )
            c_i = int( col[1] )
            dr =  float(col[2] )
            et_value = float(col[3] ) 
            bin_v = int( dr/ bin_size )
            if( bin_v < n_bins ):
                bins[bin_v] += et_value
                et_sq = et_value*et_value 
                et_sq_bins[bin_v] += et_sq
                bin_cnt[bin_v] += 1
                grp_et[n_i][c_i] = et_value
                
                # Find molecule numbers 
                Ni_o = group_index_i[n_i]
                atom_i = group_list_i[Ni_o]
                mol_i = MOLNUMB[atom_i]
                
                # Find molecule numbers 
                Nj_o = group_index_i[c_i]
                atom_j = group_list_i[Nj_o]
                mol_j = MOLNUMB[atom_j]
                
                
                if(  mol_i == mol_j):
                    et_sq_bins_intra[bin_v] += et_sq
                    bin_cnt_intra[bin_v] += 1
                else:
                    et_sq_bins_inter[bin_v] += et_sq
                    bin_cnt_inter[bin_v] += 1
                    

    et_cut = 0.0000001
    for g_i in range( 500):
        for g_j in range( 500):
            if( grp_et[g_i][g_j] > et_cut and   grp_et[g_j][g_i] > et_cut  ):
                print g_i,g_j,grp_et[g_i][g_j] ,  grp_et[g_j][g_i] 
    

    SUM_rhog_v = 0.0
    SUM_rsq_rhog_v = 0.0
    SUM_rqd_rhog_v = 0.0
    SUM_rhog_v_intra = 0.0
    SUM_rsq_rhog_v_intra = 0.0
    SUM_rqd_rhog_v_intra = 0.0
    SUM_rhog_v_inter = 0.0
    SUM_rsq_rhog_v_inter = 0.0
    SUM_rqd_rhog_v_inter = 0.0
    
    f_et = open("et.hist","w")
    f_et_h = open("et_h.hist","w")
    f_et.write("# r_ij ; cnts ;  sum V_{AB} ; < V_{AB} >  ; sum V_{AB}^2 ; < V_{AB}^2 > ; g_v ; < V_{AB}^2 >_intra ; g_v_intra ; < V_{AB}^2 >_inter ; g_v_inter ; 4 pi r^2 g_v")
    #r_ij,bin_cnt[b_indx],b_sum, b_average, et_sq ,et_sq_ave, rhog_v ,et_sq_intra_ave,rhog_v_intra ,et_sq_inter_ave, rhog_v_inter 
    for b_indx in range(n_bins):
        r_ij = b_indx*bin_size + v_floor + bin_size*0.5
        r_in = r_ij - bin_size*0.5
        r_out = r_ij + bin_size*0.5
        
        if(  bin_cnt[b_indx] > 0 ):
            b_sum =  bins[b_indx] 
            b_average = b_sum/float( bin_cnt[b_indx] )

            et_sq = et_sq_bins[b_indx]
            et_sq_ave = et_sq/float( bin_cnt[b_indx] )
            
            dr_vol = 4.0*math.pi/3.0*( r_out**3 - r_in**3 )
            
            #g_v = et_sq/dr_vol/float( 2424 )
            dVg_v = 2.0*et_sq/float( group_cnt )
            rhog_v = dVg_v/dr_vol
            
            SUM_rhog_v += rhog_v
            SUM_rsq_rhog_v += r_ij*r_ij*rhog_v
            SUM_rqd_rhog_v += r_ij*r_ij*r_ij*r_ij*rhog_v

            if(  et_sq_bins_intra[b_indx] > 0 ):
                    
                et_sq_intra = et_sq_bins_intra[b_indx]
                et_sq_intra_ave = et_sq_intra/float( bin_cnt_intra[b_indx] )
                dVg_v_intra = 2.0*et_sq_intra/float( group_cnt )
                rhog_v_intra = dVg_v_intra/dr_vol
                SUM_rhog_v_intra += rhog_v_intra
                SUM_rsq_rhog_v_intra += r_ij*r_ij*rhog_v_intra
                SUM_rqd_rhog_v_intra += r_ij*r_ij*r_ij*r_ij*rhog_v_intra
                
            else:
                et_sq_intra_ave = 0.0
                rhog_v_intra = 0.0 
                dVg_v_intra = 0.0
            

            if(  bin_cnt_inter[b_indx] > 0 ):
                et_sq_inter = et_sq_bins_inter[b_indx]
                et_sq_inter_ave = et_sq_inter/float( bin_cnt_inter[b_indx] )
                dVg_v_inter = 2.0*et_sq_inter/float( group_cnt )
                rhog_v_inter = dVg_v_inter/dr_vol
                SUM_rhog_v_inter += rhog_v_inter
                SUM_rsq_rhog_v_inter += r_ij*r_ij*rhog_v_inter
                SUM_rqd_rhog_v_inter += r_ij*r_ij*r_ij*r_ij*rhog_v_inter
                
            else:
                et_sq_inter_ave = 0.0
                rhog_v_inter = 0.0 
                dVg_v_inter = 0.0 
                
            #ket_value = ket(b_average)q
        else:
            b_sum = 0.0 
            b_average = 0.0 
            et_sq  = 0.0 
            et_sq_ave = 0.0 
            rhog_v  = 0.0
            et_sq_intra_ave = 0.0 
            rhog_v_intra  = 0.0 
            et_sq_inter_ave = 0.0 
            rhog_v_inter= 0.0
            dVg_v = 0.0
            dVg_v_inter = 0.0 
            dVg_v_intra = 0.0
            
            
            
        f_et.write( "\n %f %d %f %f %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f  %16.12f  " %  (r_ij , bin_cnt[b_indx],b_sum, b_average, et_sq ,et_sq_ave, rhog_v ,et_sq_intra_ave,rhog_v_intra ,et_sq_inter_ave, rhog_v_inter ))
        f_et_h.write( "\n %f %d %f %f %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f  %16.12f %16.12f %16.12f %16.12f " %  (r_in,bin_cnt[b_indx],b_sum, b_average, et_sq ,et_sq_ave, rhog_v , et_sq_intra_ave,rhog_v_intra ,et_sq_inter_ave, rhog_v_inter,dVg_v,dVg_v_intra,dVg_v_inter ))
        f_et_h.write( "\n %f %d %f %f %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f  %16.12f %16.12f %16.12f %16.12f " %  (r_out,bin_cnt[b_indx],b_sum, b_average, et_sq ,et_sq_ave, rhog_v , et_sq_intra_ave,rhog_v_intra ,et_sq_inter_ave, rhog_v_inter,dVg_v,dVg_v_intra,dVg_v_inter ))

    f_et.close()
    f_et_h.close()
    
    
    N_pairs = sum(bin_cnt )
    N_Aq = N_pairs/volume_i
    
    print "   total pairs ",N_pairs,N_Aq,"  x10^-{21} "
    print "   Total: "
    print "     sum rho g_v ",SUM_rhog_v
    print "     sum r^2 rho g_v ",SUM_rsq_rhog_v
    print "     sum r^4 rho g_v ",SUM_rqd_rhog_v
    ave_l_sq = SUM_rqd_rhog_v/SUM_rsq_rhog_v
    ave_l = math.sqrt(ave_l_sq )
    print "     <l^2> ",SUM_rsq_rhog_v/SUM_rhog_v,ave_l
    
    #
    # Inter mol 
    #
    
    N_pairs_inter = sum(bin_cnt_inter )
    N_Aq_inter = N_pairs_inter/volume_i*100
    
    print "   Inter molecular pair : ", N_pairs_inter,N_Aq_inter,"  x10^-{19}/cm^3  ",100.0*N_pairs_inter/N_pairs
    print "     sum rho g_v ",SUM_rhog_v_inter
    print "     sum r^2 rho g_v ",SUM_rsq_rhog_v_inter
    print "     sum r^4 rho g_v ",SUM_rqd_rhog_v_inter
    ave_l_sq = SUM_rqd_rhog_v_inter/SUM_rsq_rhog_v_inter
    ave_l = math.sqrt(ave_l_sq )
    print "     <l^2> ",ave_l_sq,ave_l
    
    #
    # Intra mol 
    #
    
    N_pairs_intra = sum(bin_cnt_intra )
    N_Aq_intra = N_pairs_intra/volume_i*100
    
    print "   Intra molecular pair : ", N_pairs_intra,N_Aq_intra,"  x10^-{19}/cm^3 ",100.0*N_pairs_intra/N_pairs
    print "     sum rho g_v ",SUM_rhog_v_intra
    print "     sum r^2 rho g_v ",SUM_rsq_rhog_v_intra
    print "     sum r^4 rho g_v ",SUM_rqd_rhog_v_intra
    ave_l_sq = SUM_rqd_rhog_v_intra/SUM_rsq_rhog_v_intra
    ave_l = math.sqrt(ave_l_sq )
    print "     <l^2> ",ave_l_sq,ave_l
    

if __name__=="__main__":
    main()