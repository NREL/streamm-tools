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

def sigma_m(N,ave,ave_sq):
    """
    Calculate the standard deviation of the mean for a confidence
    interval of 95%. Will return zero for single valued data sets 
    """
    import numpy
    
    # Some website that probably does not exist
    #   http://mathworld.wolfram.com/Studentst-Distribution.html
    
    v = N - 1 #  Degrees of freedom
    
    # Set appropriate Student t prefix
    if( v > 100 ):
	Coefint_pre = 1.64487
    elif( v > 30 ):
	Coefint_pre = 1.66023
    elif( v > 10 ):
	Coefint_pre = 1.69726
    elif( v > 5 ):
	Coefint_pre = 1.81246
    elif( v == 5 ):
	Coefint_pre = 2.01505
    elif( v == 4 ):
	Coefint_pre = 2.13185
    elif( v == 3 ):
	Coefint_pre = 2.35336
    elif( v == 2 ):
	Coefint_pre = 2.91999
    elif( v == 1 ):
	Coefint_pre = 6.31375  
	
    if( N > 1 ):
	v_sqrt = numpy.sqrt(  N - 1 )
	sigma = numpy.sqrt(  ( ave_sq ) - (ave)**2 ) # Standard deviation
	sigma_m = Coefint_pre*sigma/v_sqrt
    else:
	sigma_m = 0.0  # Set zero for unknow error 
    
    return sigma_m

def main():
    import os, sys, numpy , math 
    import file_io, gromacs , top, elements, groups, prop
    
    debug = 0
    
    # Get options
    options, args = get_options()
    
    # Get reference structure file 
    if( len(options.in_gro) ):
        if( options.verbose ): print "  Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options,options.in_gro)
        
    # Read in topology file
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

    # Read in values for pairs of groups 
    et_file = options.in_et
    if( file_io.file_exists(et_file)):
        f_et_in = open(et_file,"r")
        Lines = f_et_in.readlines()
        f_et_in.close()

    # Initialize histogram settings
    bin_size = 0.1 #(angstroms)
    v_min = 0.0
    v_max = 20.0
    v_floor = int(v_min/bin_size)*bin_size 
    val_range = v_max - v_floor
    n_bins = int(val_range/bin_size) + 1

    bins = numpy.zeros(n_bins)
    et_sq_bins  = numpy.zeros(n_bins)
    et_qu_bins  = numpy.zeros(n_bins)
    bin_cnt = numpy.zeros(n_bins)

    et_sq_bins_intra  = numpy.zeros(n_bins)
    et_qu_bins_intra  = numpy.zeros(n_bins)
    bin_cnt_intra = numpy.zeros(n_bins)
    et_sq_bins_inter  = numpy.zeros(n_bins)
    et_qu_bins_inter  = numpy.zeros(n_bins)
    bin_cnt_inter = numpy.zeros(n_bins)
    
    
    grp_et = numpy.empty((3000,3000),dtype="float")
    
    print len(grp_et)
    print len(grp_et[0])
    print grp_et[0][0]

    # Find volume
    volume_i =   prop.volume( LV ) 

    for line in Lines:
        col = line.split()
        if( col[0] != "#" ):
            n_i = int( col[0] )   # Nuetral group #
            c_i = int( col[1] )   # Cation group #
            dr =  float(col[2] )  # inter nitrogen spacing 
            et_value = float(col[3] )     # value of V_{AB} from NWCHEM
            bin_v = int( dr/ bin_size )   # Bin value 
            if( bin_v < n_bins ):         
                bin_cnt[bin_v] += 1
		#
                bins[bin_v] += et_value
                et_sq = et_value*et_value 
                et_qu = et_sq*et_sq 
                et_sq_bins[bin_v] += et_sq
                et_qu_bins[bin_v] += et_qu
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
                    et_qu_bins_intra[bin_v] += et_qu
                    bin_cnt_intra[bin_v] += 1
                else:
                    et_sq_bins_inter[bin_v] += et_sq
                    et_qu_bins_inter[bin_v] += et_qu
                    bin_cnt_inter[bin_v] += 1
                    
		if( debug ):
		    if( bin_v == 50 ):
			print dr , et_value

    et_cut = 0.0000001
    for g_i in range( 500):
        for g_j in range( 500):
            if( grp_et[g_i][g_j] > et_cut and   grp_et[g_j][g_i] > et_cut  ):
                print g_i,g_j,grp_et[g_i][g_j] ,  grp_et[g_j][g_i] 
    
    #
    # Calaculate averages, g_v, 4 pi 
    #
    #  b - bin number
    #  N - number of groups
    # dV - volume element  (A^3)
    #
    #  < V_{AB} > (b)   = sum(VAB)(b) / cnt(b) 
    #  < V_{AB}^2 > (b) = sum(VAB^2)(b) / cnt(b) 
    #
    #             dVg_v = 2 * sum(VAB^2) / N 
    #            dr_vol = dV 
    #            rhog_v = dVg_v/dr_vol            # (rho * g_v )
    #   
    #  
    
    
    FOURPI =  4.0*math.pi
    VOL_CONST = FOURPI/3.0
    
    Coef_interval =  1.00 # 96
    
    SUM_rhog_v = 0.0
    SUM_rsq_rhog_v = 0.0
    SUM_rsq_rhog_v_sig_mean = 0.0
    SUM_rqd_rhog_v = 0.0
    SUM_rqd_rhog_v_sig_mean = 0.0
    SUM_rhog_v_intra = 0.0
    SUM_rsq_rhog_v_intra = 0.0
    SUM_rsq_rhog_v_intra_sig_mean = 0.0
    SUM_rqd_rhog_v_intra = 0.0
    SUM_rqd_rhog_v_intra_sig_mean = 0.0
    SUM_rhog_v_inter = 0.0
    SUM_rsq_rhog_v_inter = 0.0
    SUM_rsq_rhog_v_inter_sig_mean = 0.0
    SUM_rqd_rhog_v_inter = 0.0
    SUM_rqd_rhog_v_inter_sig_mean = 0.0
    
    f_VAB = open("VAB.hist","w")
    f_VABsq = open("VABsq.hist","w")
    f_gv = open("rhogv.hist","w")
    f_rsqgv = open("rsqgv.hist","w")
    f_rqugv = open("rqugv.hist","w")
    f_dVgv = open("dVgv.hist","w")
    #f_et_h = open("et_h.hist","w")
    #f_et.write("# r_ij ; cnts ;  sum V_{AB} ; < V_{AB} >  ; sum V_{AB}^2 ; < V_{AB}^2 > ; g_v ; < V_{AB}^2 >_intra ; g_v_intra ; < V_{AB}^2 >_inter ; g_v_inter ; 4 pi r^2 g_v")
    f_VAB.write("# r_ij ; cnts ; <V_{AB}>(eV);  mstd( <V_{AB}>(eV)  )  ")
    f_VABsq.write("# r_ij ; cnts ; < V_{AB}^2 >  (eV^2); mstd( <V_{AB}^2>(eV^2)  ),   intra  < V_{AB}^2 >  (eV^2); intra mstd( <V_{AB}^2>(eV^2)  ), inter  < V_{AB}^2 >  (eV^2); inter mstd( <V_{AB}^2>(eV^2)  ) ")
    f_gv.write("# r_ij ; cnts ; <rhog_v> (eV^2/A^3);  mstd( <rhog_v>(eV^2/A^3)  ); intra <rhog_v> (eV^2/A^3); intra  mstd( <rhog_v>(eV^2/A^3)  ); inter <rhog_v> (eV^2/A^3); inter mstd( <rhog_v>(eV^2/A^3) )  ")
    f_rsqgv.write("# r_ij ; cnts ; <rhog_v> (eV^2/A);  mstd( <rhog_v>(eV^2/A)  ); intra <rhog_v> (eV^2/A); intra  mstd( <rhog_v>(eV^2/A)  ); inter <rhog_v> (eV^2/A); inter mstd( <rhog_v>(eV^2/A) )  ")
    f_rqugv.write("# r_ij ; cnts ; <rhog_v> (eV^2 A);  mstd( <rhog_v>(eV^2 A)  ); intra <rhog_v> (eV^2 A); intra  mstd( <rhog_v>(eV^2 A)  ); inter <rhog_v> (eV^2 A); inter mstd( <rhog_v>(eV^2 A) )  ")
    f_dVgv.write("# r_ij ; cnts ; <dVrhog_v> (eV^2);  mstd( <dVrhog_v> (eV^2)  ); intra <dVrhog_v> (eV^2); intra  mstd(<dVrhog_v> (eV^2)); inter <dVrhog_v> (eV^2); inter mstd( <dVrhog_v> (eV^2) )  ")
    #r_ij,bin_cnt[b_indx],b_sum, b_average, et_sq ,et_sq_ave, rhog_v ,et_sq_intra_ave,rhog_v_intra ,et_sq_inter_ave, rhog_v_inter 
    for b_indx in range(n_bins):
        r_ij = b_indx*bin_size + v_floor + bin_size*0.5
        r_in = r_ij - bin_size*0.5
        r_out = r_ij + bin_size*0.5
        N_cnts = bin_cnt[b_indx]
	
        if(  N_cnts > 0 ):
	    
	    dr_sq = r_ij*r_ij
	    dr_qu = dr_sq*dr_sq
	    N_sqrt =  numpy.sqrt(  N_cnts)
	    
            et_sum =  bins[b_indx] 
            et_average = et_sum/float(N_cnts)

            et_sq = et_sq_bins[b_indx]
            et_sq_ave = et_sq/float(N_cnts)
            
            et_qu = et_qu_bins[b_indx]
            et_qu_ave = et_qu/float(N_cnts)
            
	    et_sig_mean = sigma_m(N_cnts,et_average,et_sq_ave)
	    et_sq_sig_mean = sigma_m(N_cnts,et_sq_ave,et_qu_ave)
		
            dr_vol = VOL_CONST*( r_out**3 - r_in**3 )
            
            #g_v = et_sq/dr_vol/float( 2424 )
            dVg_v = 2.0*et_sq/float( group_cnt )
	    dVg_v_sig_mean = 2.0*et_sq_sig_mean*float(N_cnts)/float( group_cnt )
	    
            rhog_v = dVg_v/dr_vol
	    rhog_v_sig_mean = dVg_v_sig_mean/dr_vol
	    
	    fourpi_rsqr_rhog_v = FOURPI*dr_sq*rhog_v
	    fourpi_rsqr_rhog_v_sig_mean = FOURPI*dr_sq*rhog_v_sig_mean 
	    
            SUM_rhog_v += rhog_v
            SUM_rsq_rhog_v += dr_sq*rhog_v
	    SUM_rsq_rhog_v_sig_mean  += dr_sq*rhog_v_sig_mean*dr_sq*rhog_v_sig_mean
            SUM_rqd_rhog_v += dr_qu*rhog_v
	    SUM_rqd_rhog_v_sig_mean += dr_qu*rhog_v_sig_mean*dr_qu*rhog_v_sig_mean

	    N_cnt_intra = float( bin_cnt_intra[b_indx] )
	    N_cnt_inter = float( bin_cnt_inter[b_indx] )
            if(  N_cnt_intra> 0 ):
                    
                et_sq_intra = et_sq_bins_intra[b_indx]
                et_sq_intra_ave = et_sq_intra/N_cnt_intra
                et_qu_intra = et_qu_bins_intra[b_indx]
                et_qu_intra_ave = et_qu_intra/N_cnt_intra
		
		et_sq_sig_mean_intra = sigma_m(N_cnt_intra,et_sq_intra_ave,et_qu_intra_ave)
		    
                dVg_v_intra = 2.0*et_sq_intra/float( group_cnt )
		dVg_v_sig_mean_intra = 2.0*et_sq_sig_mean_intra*N_cnt_intra/float( group_cnt )
                rhog_v_intra = dVg_v_intra/dr_vol
		rhog_v_sig_mean_intra = dVg_v_sig_mean_intra/dr_vol
		
		
		fourpi_rsqr_rhog_v_intra = FOURPI*dr_sq*rhog_v_intra
		fourpi_rsqr_rhog_v_sig_mean_intra = FOURPI*dr_sq*rhog_v_sig_mean_intra 
	    
                SUM_rhog_v_intra += rhog_v_intra
                SUM_rsq_rhog_v_intra += dr_sq*rhog_v_intra
		SUM_rsq_rhog_v_intra_sig_mean += dr_sq*rhog_v_sig_mean_intra*dr_sq*rhog_v_sig_mean_intra
                SUM_rqd_rhog_v_intra += dr_qu*rhog_v_intra
                SUM_rqd_rhog_v_intra_sig_mean += dr_qu*rhog_v_sig_mean_intra*dr_qu*rhog_v_sig_mean_intra
		
                
            else:
                et_sq_intra_ave = 0.0
		et_sq_sig_mean_intra = 0.0 
                rhog_v_intra = 0.0
		fourpi_rsqr_rhog_v_intra = 0.0
		fourpi_rsqr_rhog_v_sig_mean_intra = 0.0
		
                dVg_v_intra = 0.0
		dVg_v_sig_mean_intra = 0.0 

            if(  N_cnt_inter > 0 ):
                et_sq_inter = et_sq_bins_inter[b_indx]
                et_sq_inter_ave = et_sq_inter/N_cnt_inter
                et_qu_inter = et_qu_bins_inter[b_indx]
                et_qu_inter_ave = et_qu_inter/N_cnt_inter
		
		et_sq_sig_mean_inter = sigma_m(N_cnt_inter,et_sq_inter_ave,et_qu_inter_ave)
		
                dVg_v_inter = 2.0*et_sq_inter/float( group_cnt )
		dVg_v_sig_mean_inter = 2.0*et_sq_sig_mean_inter*N_cnt_inter/float( group_cnt )
                rhog_v_inter = dVg_v_inter/dr_vol
		rhog_v_sig_mean_inter = dVg_v_sig_mean_inter/dr_vol
		
		fourpi_rsqr_rhog_v_inter = FOURPI*dr_sq*rhog_v_inter
		fourpi_rsqr_rhog_v_sig_mean_inter = FOURPI*dr_sq*rhog_v_sig_mean_inter 
	    
		
                SUM_rhog_v_inter += rhog_v_inter
                SUM_rsq_rhog_v_inter += dr_sq*rhog_v_inter
                SUM_rsq_rhog_v_inter_sig_mean += dr_sq*rhog_v_sig_mean_inter*dr_sq*rhog_v_sig_mean_inter
                SUM_rqd_rhog_v_inter += dr_qu*rhog_v_inter
                SUM_rqd_rhog_v_inter_sig_mean += dr_qu*rhog_v_sig_mean_inter*dr_qu*rhog_v_sig_mean_inter
                
            else:
                et_sq_inter_ave = 0.0
		et_sq_sig_mean_inter = 0.0 
                rhog_v_inter = 0.0
		rhog_v_sig_mean_inter = 0.0 
		fourpi_rsqr_rhog_v_inter = 0.0
		fourpi_rsqr_rhog_v_sig_mean_inter = 0.0 
                dVg_v_inter = 0.0
		dVg_v_sig_mean_inter = 0.0  
		
            #ket_value = ket(b_average)q
        else:
            et_sum = 0.0 
            et_average = 0.0
	    et_sig_mean = 0.0
	    #
            et_sq  = 0.0 
            et_sq_ave = 0.0
	    et_sq_sigma = 0.0 
            rhog_v  = 0.0
	    rhog_v_sig_mean = 0.0
	    fourpi_rsqr_rhog_v = 0.0
	    fourpi_rsqr_rhog_v_sig_mean = 0.0
            dVg_v = 0.0
	    dVg_v_sig_mean = 0.0 
	    #
	    # Set intra values
	    #
            et_sq_intra_ave = 0.0
	    et_sq_sig_mean_intra = 0.0
	    
            dVg_v_intra = 0.0
            dVg_v_sig_mean_intra = 0.0 
            	    
            rhog_v_intra  = 0.0
	    rhog_v_sig_mean_intra = 0.0 
	    fourpi_rsqr_rhog_v_intra = 0.0
	    fourpi_rsqr_rhog_v_sig_mean_intra = 0.0
		
	    # Set inter values 
            et_sq_inter_ave = 0.0
	    et_sq_sig_mean_inter = 0.0 
            rhog_v_inter= 0.0
	    rhog_v_sig_mean_inter = 0.0
	    fourpi_rsqr_rhog_v_inter = 0.0
	    fourpi_rsqr_rhog_v_sig_mean_inter = 0.0 
            dVg_v_inter = 0.0
	    dVg_v_sig_mean_inter = 0.0
	    
        f_VAB.write( "\n %f %d %16.12f  %16.12f  " % (r_ij , bin_cnt[b_indx], et_average, et_sig_mean ) )
        f_VABsq.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , bin_cnt[b_indx], et_sq_ave, et_sq_sigma, et_sq_intra_ave, et_sq_sig_mean_intra, et_sq_inter_ave, et_sq_sig_mean_inter  ) )

        f_gv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , bin_cnt[b_indx], rhog_v, rhog_v_sig_mean, rhog_v_intra, rhog_v_sig_mean_intra, rhog_v_sig_mean_inter, rhog_v_sig_mean_inter  ) )
        f_rsqgv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , bin_cnt[b_indx], fourpi_rsqr_rhog_v, fourpi_rsqr_rhog_v_sig_mean, fourpi_rsqr_rhog_v_intra, fourpi_rsqr_rhog_v_sig_mean_intra, fourpi_rsqr_rhog_v_inter, fourpi_rsqr_rhog_v_sig_mean_inter  ) )
        f_dVgv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , bin_cnt[b_indx], dVg_v, dVg_v_sig_mean, fourpi_rsqr_rhog_v_intra, fourpi_rsqr_rhog_v_sig_mean_intra, fourpi_rsqr_rhog_v_inter, fourpi_rsqr_rhog_v_sig_mean_inter  ) )

	if( debug ):
	    if( b_indx == 50 ):
		print "N_cnts",bin_cnt[b_indx]
		print " et_sum ",et_sum
		print " et_average ",et_average
		print " et_sq ",et_sq
		print " et_sq_ave ",et_sq_ave
		print " et_qu_ave ",et_qu_ave
		print " et_sq_sig_mean ",et_sq_sig_mean
		print " r_out, r_in ",r_out, r_in
		print " dr_vol ",dr_vol
		print " group_cnt ",group_cnt
		print " dVg_v ",dVg_v
		print " dVg_v_sig_mean ",dVg_v_sig_mean
		print " rhog_v ",rhog_v
		print " rhog_v_sig_mean ",rhog_v_sig_mean
		print " fourpi_rsqr_rhog_v ",fourpi_rsqr_rhog_v
		print " fourpi_rsqr_rhog_v_sig_mean ",fourpi_rsqr_rhog_v_sig_mean
		
		print et_sq_ave, et_sq_sigma
		print rhog_v, rhog_v_sig_mean
		print fourpi_rsqr_rhog_v, fourpi_rsqr_rhog_v_sig_mean
		print dVg_v, dVg_v_sig_mean
		sys.exit(" debug mode! ")

    f_VAB.close()
    f_VABsq.close()
    f_gv.close()
    f_rsqgv.close()
    f_rqugv.close()
    f_dVgv.close()
    
    N_pairs = sum(bin_cnt )
    N_Aq = N_pairs/volume_i
    
    print "   total pairs ",N_pairs,N_Aq,"  x10^-{21} "
    print "   Total: "
    print "     sum rho g_v ",SUM_rhog_v
    print "     sum r^2 rho g_v ",SUM_rsq_rhog_v
    print "     sum r^4 rho g_v ",SUM_rqd_rhog_v
    ave_l_sq = SUM_rqd_rhog_v/SUM_rsq_rhog_v
    ave_l_sq_sig_mean = numpy.sqrt( SUM_rqd_rhog_v_sig_mean) / numpy.sqrt( SUM_rsq_rhog_v_sig_mean )
    
    ave_l = numpy.sqrt(ave_l_sq )
    ave_l_sig_mean = numpy.sqrt(ave_l_sq_sig_mean )
    print "     <l^2> ",ave_l_sq, "+-",ave_l_sq_sig_mean,"     l ",ave_l, "+-",ave_l_sig_mean
    
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
    
    ave_l_sq_sig_mean = numpy.sqrt( SUM_rqd_rhog_v_intra_sig_mean) / numpy.sqrt( SUM_rsq_rhog_v_intra_sig_mean )
    
    ave_l = numpy.sqrt(ave_l_sq )
    ave_l_sig_mean = numpy.sqrt(ave_l_sq_sig_mean )
    print "    Intra-molecular <l^2> ",ave_l_sq, "+-",ave_l_sq_sig_mean,"     l ",ave_l, "+-",ave_l_sig_mean
        
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
    ave_l_sq_sig_mean = numpy.sqrt( SUM_rqd_rhog_v_inter_sig_mean) / numpy.sqrt( SUM_rsq_rhog_v_inter_sig_mean )
    
    ave_l = numpy.sqrt(ave_l_sq )
    ave_l_sig_mean = numpy.sqrt(ave_l_sq_sig_mean )
    print "    Inter-molecular <l^2> ",ave_l_sq, "+-",ave_l_sq_sig_mean,"     l ",ave_l, "+-",ave_l_sig_mean
        
    
if __name__=="__main__":
    main()
