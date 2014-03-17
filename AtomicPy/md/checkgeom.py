#! /usr/bin/env python

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# Email travis.kemper@nrel.gov
# Version 2.00 

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    
    parser.add_option("-a","--accuracy", dest="accuracy", default=0.00001, help="Verbose output ")

    (options, args) = parser.parse_args()
        
    return options, args


def checkredundant(R_1, GT_1, AT_1, AMASS_1, BONDS_1,ANGLES_1,DIHS_1,R_2, GT_2, AT_2, AMASS_2, BONDS_2,ANGLES_2,DIHS_2):
    import sys
    
    for i in range( len( AT_1)):
        if ( AT_1[i] != AT_2[i] ):
            print ' atom types differ '
            sys.exit('invalad geometry comparsion ')
    print ' bonds ',len( BONDS_1)
    for i in range( len( BONDS_1)):
        if ( BONDS_1[i][0] != BONDS_2[i][0] and  BONDS_1[i][1] != BONDS_2[i][1] ):
            print ' bonds differ '
            sys.exit('invalad geometry comparsion ')
    for i in range( len( ANGLES_1)):
        if ( ANGLES_1[i][0] != ANGLES_2[i][0] and  ANGLES_1[i][1] != ANGLES_2[i][1] and  ANGLES_1[i][2] != ANGLES_2[i][2] ):
            print ' angles differ '
            sys.exit('invalad geometry comparsion ')
    for i in range( len( DIHS_1)):
        if ( DIHS_1[i][0] != DIHS_2[i][0] and  DIHS_1[i][1] != DIHS_2[i][1] and  DIHS_1[i][2] != DIHS_2[i][2] and  DIHS_1[i][3] != DIHS_2[i][3]):
            print ' dihedrals differ '
            sys.exit('invalad geometry comparsion ')


def check_bonds(R_1,BONDS_1, R_2,BONDS_2):
    import sys, numpy 
    
    print ' bonds ',len( BONDS_1)
    for i in range( len( BONDS_1)):
        i_1 = BONDS_1[i][0]
        j_1 = BONDS_1[i][1]
        
        r_i = R_1[i_1]
        r_j = R_1[j_1]
        r_ij_1 = r_j - r_i
        R_ij_1.append(  numpy.linalg.norm(  r_ij_1 ))
        
        i_2 = BONDS_2[i][0]
        j_2 = BONDS_2[i][1] 

        r_i = R_2[i_2]
        r_j = R_2[j_2]
        r_ij_2 = r_j - r_i
        R_ij_2.append( numpy.linalg.norm(  r_ij_2 ))
        
        
    return  R_ij_1,R_ij_2

def check_angles(R_1, ANGLES_1 ,R_2, ANGLES_2 ):
    import sys, math , prop
    #
    #  k - i - j
    #
    debug = 0
    
    print ' bonds ',len( ANGLES_1)
    for angle_indx in range( len( ANGLES_1)):
        a_k = ANGLES_1[angle_indx][0]
        a_i = ANGLES_1[angle_indx][1]
        a_j = ANGLES_1[angle_indx][2]
        r_k = R_1[a_k]
        r_i = R_1[a_i]
        r_j = R_1[a_j]
        r_ij = r_j - r_i
        r_ik = r_k - r_i
        angle_1 = prop.getAngle(r_ik,r_ij)
        ANGLE_kij_1.append( angle_1 )
        if( debug) :
            print  ' angle ', angle_indx, a_k,a_i,a_j,math.acos( cos_1 )*180.0/math.pi
            #print ' r_i ', r_i
            #print ' r_j ', r_j
            ##print ' r_k ', r_k
            #print ' r_ij ',r_ij
            #print ' r_ik ',r_ik

        a_k = ANGLES_2[angle_indx][0]
        a_i = ANGLES_2[angle_indx][1]
        a_j = ANGLES_2[angle_indx][2]
        r_k = R_2[a_k]
        r_i = R_2[a_i]
        r_j = R_2[a_j]
        r_ij = r_j - r_i
        r_ik = r_k - r_i
        angle_2 = prop.getAngle(r_ik,r_ij)
        ANGLE_kij_2.append( angle_2 )

    if(debug): sys.exit( 'check_angles')
    
    return ANGLE_kij_1,ANGLE_kij_2


def check_dihedrals(R_1, DIHS_1, R_2, DIHS_2):
    import sys, numpy , math, prop 
    #
    #  k
    #    \
    #      i  -  j
    #             \
    #              l


    debug = 0
    for dih_indx in range( len( DIHS_1)):
        a_k = DIHS_1[dih_indx][0]
        a_i = DIHS_1[dih_indx][1]
        a_j = DIHS_1[dih_indx][2]
        a_l = DIHS_1[dih_indx][3]
        
        r_k = R_1[a_k]
        r_i = R_1[a_i]
        r_j = R_1[a_j]
        r_l = R_1[a_l]
        
	#r_ik = r_k - r_i
	#r_ij = r_j - r_i
	#r_il = r_l - r_j
		
	dih_angle = prop.getDihedral(r_k,r_i,r_j,r_l)
                
        DIH_kijl_1.append(  dih_angle )

        if( debug):
            print 'dih_angle',dih_angle
            print ' x', x_il_ij
            print '  r_ij', r_jk 
            print ' x r_ij', vector_inner( x_il_ij, r_jk ) 
            print " index ",a_l," or index ",a_i," or index ",a_j," or index ",a_k,dih_angle
        
        a_k = DIHS_2[dih_indx][0]
        a_i = DIHS_2[dih_indx][1]
        a_j = DIHS_2[dih_indx][2]
        a_l = DIHS_2[dih_indx][3]
        
        r_k = R_2[a_k]
        r_i = R_2[a_i]
        r_j = R_2[a_j]
        r_l = R_2[a_l]
        
	dih_angle = prop.getDihedral(r_k,r_i,r_j,r_l)
                
        DIH_kijl_2.append(  dih_angle )

    if(debug): sys.exit( 'check_dihedrals')

    return DIH_kijl_1 , DIH_kijl_2

def print_check(sys_name,AT_1,GT_1,BONDS_1,R_ij_1,ANGLES_1, ANGLE_kij_1,R_2,ANGLES_2,ANGLE_kij_2,R_ij_2,DIHS_1, DIH_kijl_1, DIHS_2 , DIH_kijl_2):
    import math,sys

    d_sq = 0.0
    cnt = 0
    # print data to latex table file 
    Ftab = open( 'rmsd.tab', 'w')
    Ftab_mult = open( 'rmsd_mult.tab', 'w')
    #Ftab.write( " \\begin{tabular}{c c c c} \n")
    Ftab.write( " $ %10s $ &  &  & \\\  \n" % (sys_name))
    Ftab.write( "  type & RMSD & maximum & atoms \\\  \n" )
    Ftab.write( " \hline \n" )
    #
    F = open( 'bonds.dif' , 'w' )
    F.write( " Bonds %d \n" % ( len(BONDS_1) )  )
    delta_max = 0.0 
    delta_sq_max = 0.0 
    for i in range( len( BONDS_1)):
        i_1 = BONDS_1[i][0]
        j_1 = BONDS_1[i][1]
        i_max = i_1
        j_max = j_1
        no_h = 1
        include_i = 1
        if( no_h ):
            if( AT_1[i_1][:1] == 'H' ) :  include_i = 0 
            if( AT_1[j_1][:1] == 'H' ) :  include_i = 0 
            if( AT_1[i_1] == 'LP' ) :  include_i = 0 
            if( AT_1[j_1] == 'LP' ) :  include_i = 0 
        if( include_i ):
            print ' Calculating ', AT_1[i_1][0:],AT_1[j_1][:1]
            delta = R_ij_2[i] - R_ij_1[i]
            err = delta
            delta_sq =  delta*delta
            d_sq = d_sq + delta_sq 
            cnt = cnt + 1 
            F.write( " index %7d or index  %7d %12.4f %12.4f %12.4f \n" % (i_1,j_1,R_ij_1[i],R_ij_2[i],delta) )
            if ( delta_sq > delta_sq_max  ):
                delta_sq_max =  delta_sq
                delta_max = delta
                i_max = i_1
                j_max = j_1
    d_ave = d_sq / float(cnt)
    rmsd = math.sqrt( d_ave )
    print ' Bond RMSD ',rmsd
    Ftab.write( " bond length ($\AA$) & %12.4f  & %12.4f & %s - %s  \\\ \n" % ( rmsd , delta_max , GT_1[i_max] , GT_1[j_max]  ))
    Ftab_mult.write( " %12.4f  & %12.4f & %s - %s  &  \n" % ( rmsd , delta_max , GT_1[i_max] , GT_1[j_max]  ))
    F.write( '\n' )
    
    d_sq = 0.0
    cnt = 0
    F.write( " Angles %d \n" % ( len(ANGLES_1) )  )
    delta_max = 0.0 
    delta_sq_max = 0.0 
    for angle_indx in range( len( ANGLES_1)):
        a_k = ANGLES_1[angle_indx][0]
        a_i = ANGLES_1[angle_indx][1]
        a_j = ANGLES_1[angle_indx][2]
        i_max = a_k
        j_max = a_i
        k_max = a_j
        no_h = 1
        include_i = 1
        if( no_h ):
            if( AT_1[a_k][:1] == 'H' ) :  include_i = 0 
            if( AT_1[a_i][:1] == 'H' ) :  include_i = 0 
            if( AT_1[a_j][:1] == 'H' ) :  include_i = 0 
            if( AT_1[a_k] == 'LP' ) :  include_i = 0 
            if( AT_1[a_i] == 'LP' ) :  include_i = 0 
            if( AT_1[a_j] == 'LP' ) :  include_i = 0 
        if( include_i ):
            delta = ANGLE_kij_2[angle_indx] - ANGLE_kij_1[angle_indx]
            err = delta
            delta_sq =  delta*delta
            d_sq = d_sq + delta_sq 
            cnt = cnt + 1 
            F.write( "index %7d or index %7d or index %7d %12.4f %12.4f %12.4f \n" % (a_k,a_i,a_j,ANGLE_kij_1[angle_indx] , ANGLE_kij_2[angle_indx],delta) )
            if ( delta_sq > delta_sq_max  ):
                delta_sq_max =  delta_sq
                delta_max = delta
                i_max = a_k
                j_max = a_i
                k_max = a_j
    d_ave = d_sq / float(cnt)
    rmsd = math.sqrt( d_ave )
    print ' Angle RMSD ',rmsd
    Ftab.write( " angle rmsd (deg) & %12.4f  & %12.4f & %s - %s - %s \\\ \n" % ( rmsd, delta_max , GT_1[i_max] , GT_1[j_max]  , GT_1[k_max]  ))
    Ftab_mult.write( " %12.4f  & %12.4f & %s - %s - %s & \n" % ( rmsd, delta_max , GT_1[i_max] , GT_1[j_max]  , GT_1[k_max]  ))
    
    F.write( '\n' )
    if( len(DIHS_1) > 0 ):
        d_sq = 0.0
        cnt = 0
        F.write( " Dihedrals %d \n" % ( len(DIHS_1) )  )
        delta_max = 0.0 
        delta_sq_max = 0.0 
        for angle_indx in range( len( DIHS_1)):
            a_l = DIHS_1[angle_indx][0]
            a_i = DIHS_1[angle_indx][1]
            a_j = DIHS_1[angle_indx][2]
            a_k = DIHS_1[angle_indx][3]
            i_max = a_k
            j_max = a_i
            k_max = a_j
            l_max = a_l
            no_h = 1
            include_i = 1
            if( no_h ):
                if( AT_1[a_k][:1] == 'H' ) :  include_i = 0 
                if( AT_1[a_i][:1] == 'H' ) :  include_i = 0 
                if( AT_1[a_j][:1] == 'H' ) :  include_i = 0 
                if( AT_1[a_l][:1] == 'H' ) :  include_i = 0 
                if( AT_1[a_k] == 'LP' ) :  include_i = 0 
                if( AT_1[a_i] == 'LP' ) :  include_i = 0 
                if( AT_1[a_j] == 'LP' ) :  include_i = 0 
                if( AT_1[a_l] == 'LP' ) :  include_i = 0 
            if( include_i ):
                delta = DIH_kijl_2[angle_indx] - DIH_kijl_1[angle_indx]
                if ( delta < -180 ): delta = delta + 360 
                if ( delta > 180 ): delta = delta - 360 
                err = delta
                delta_sq =  delta*delta
                d_sq = d_sq + delta_sq 
                cnt = cnt + 1 
                F.write( "index %7d or index %7d or index %7d or index %7d %12.4f %12.4f %12.4f \n" % (a_k,a_i,a_j,a_l,DIH_kijl_1[angle_indx] , DIH_kijl_2[angle_indx],delta) )
                if ( delta_sq > delta_sq_max ):
                    delta_sq_max =  delta_sq
                    delta_max = delta
                    i_max = a_k
                    j_max = a_i
                    k_max = a_j
                    l_max = a_l
        d_ave = d_sq / float(cnt)
        rmsd = math.sqrt( d_ave )
        print ' Dihedral RMSD ',rmsd
        Ftab.write( " dihedral rmsd (deg) & %12.4f  & %12.4f & %s - %s - %s  - %s \\\ \n" % ( rmsd, delta_max , GT_1[i_max] , GT_1[j_max] , GT_1[k_max] , GT_1[l_max]  ))
        Ftab_mult.write( "  %12.4f  & %12.4f & %s - %s - %s  - %s & \n" % ( rmsd, delta_max , GT_1[i_max] , GT_1[j_max] , GT_1[k_max] , GT_1[l_max]  ))
    
    F.write( '\n' )
    Ftab.write( " \hline   \hline  \n" )

    #Ftab.write( "\end{tabular} ")

    F.close()
    Ftab.close()
    Ftab_mult.close()
  
      
def main():
    #
    # Read in two geometries and compare differences in bonds, angles and dihedrals 
    #
    
    import sys
    #global R_1, GT_1, AT_1, AMASS_1, BONDS_1,ANGLES_1,DIHS_1
    #global R_2, GT_2, AT_2, AMASS_2, BONDS_2,ANGLES_2,DIHS_2
    # Read in first geometry
    sys_var = 1 
    REF_FILE = sys.argv[sys_var]
    print "reading in", REF_FILE
    REF_SUF = REF_FILE[-4:]
    Lines = files.getlines(REF_FILE)

    if( REF_SUF == '.gro' ):
        # read in gro
        GT_1, R_1, VEL_1, LV_1 =   read_gro(options,REF_FILE)

        # red in top
        sys_var = sys_var + 1
        REF_FILE = sys.argv[sys_var]
        REF_SUF = REF_FILE[-4:]
        if( REF_SUF == '.top' ):
            AT_1 , RESN_1 , RESID_1 , GTYPE_1 ,CHARN_1 , CHARGES_1 ,AMASS_1,BONDS_1,ANGLES_1,DIH_1, MOLNUMB_1, MOLPNT_1, MOLLIST_1 = read_top(options,REF_FILE)
            
        else :
            print ' second file should be a top file'
            sys.exit('input eror')
            
    # Read in second geometry 

    sys_var = sys_var + 1
    REF_FILE = sys.argv[sys_var]
    print "reading in", REF_FILE
    REF_SUF = REF_FILE[-4:]
    Lines = files.getlines(REF_FILE)

    if( REF_SUF == '.gro' ):
        # read in gro
        GT_2, R_2, VEL_2, LV_2 =   read_gro(options,REF_FILE)

        # red in top
        sys_var = sys_var + 1
        REF_FILE = sys.argv[sys_var]
        REF_SUF = REF_FILE[-4:]
        if( REF_SUF == '.top' ):
            AT_2 , RESN_2 , RESID_2 , GTYPE_2 ,CHARN_2 , CHARGES_2 ,AMASS_2,BONDS_2,ANGLES_2,DIH_2, MOLNUMB_2, MOLPNT_2, MOLLIST_2 = read_top(options,REF_FILE)
            
        else :
            print ' ',sys_var,'  file should be a top file'
            sys.exit('input eror')
            

    sys_var = sys_var + 1
    sys_name  = sys.argv[sys_var]
    
    # Check to make sure redundant geometrys
    print " Check structures "
    checkredundant(R_1, GT_1, AT_1, AMASS_1, BONDS_1,ANGLES_1,DIHS_1,R_2, GT_2, AT_2, AMASS_2, BONDS_2,ANGLES_2,DIHS_2)
    
    print " Checking bonds"
    R_ij_1,R_ij_2 = check_bonds(R_1, BONDS_1, R_2, BONDS_2)
    ANGLE_kij_1,ANGLE_kij_2 = check_angles(R_1, ANGLES_1 ,R_2, ANGLES_2 )
    DIH_kijl_1 , DIH_kijl_2 = check_dihedrals(R_1, DIHS_1, R_2, DIHS_2)
    print_check(sys_name,AT_1,GT_1,BONDS_1,R_ij_1,ANGLES_1, ANGLE_kij_1,R_2,ANGLES_2,ANGLE_kij_2,R_ij_2,DIHS_1, DIH_kijl_1, DIHS_2 , DIH_kijl_2)
    
    debug = 0
    if(debug):
        for i in range( len(GT_1) ):
            print AT_1[i],G
    
    
	
if __name__=="__main__":
    main()
   

    