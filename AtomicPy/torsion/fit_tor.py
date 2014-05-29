#! /usr/bin/env python
# run ab initio torsional potential 

# Dr. Travis Kemper
# NREL
# 12/09/2013
# travis.kemper@nrel.gov

EVTOKCAL = 23.0605

def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    usage = usage + "Input files \n"
    usage = usage + "  specify the destination name of an option followed by the value"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")
    
    parser.add_option("--tor_paramo", type="float", dest="tor_paramo", default=0.09, help="Value to initialize torsional potential coefficients ")
    parser.add_option("--out_ftor", dest="out_ftor", type="string", default="tor_fit.dat", help="Output fitted torsional potential data ")
    parser.add_option("--plot_tor", dest="plot_tor", type="string", default="tor_fit.pdf", help="Output fitted torsional potential graph ")
    parser.add_option("--plot_tor_comp", dest="plot_tor_comp", type="string", default="tor_comp.pdf", help="Output for components of torsional potential graph ")

    parser.add_option("--qm_tor", dest="qm_tor", type="string", default="",  help=" Data file of quantum target values ")
    parser.add_option("--ff_tor", dest="ff_tor", type="string", default="",  help=" Data file of force field values ")
    
    (options, args) = parser.parse_args()

    return options, args
    
def read_energy( energy_file ):
    import sys
    import file_io

    if ( file_io.file_exists(energy_file) ): 
        F = open(energy_file,'r')
        Lines = F.readlines()
        F.close()
    else:
        sys.exit(' data file does not exist')

    energy_data = []
    
    for line in Lines :
        col = line.split()
        if( len(col) > 0 ):
            if ( col[0][:1] != '#' ):
                energy_data.append(col)

    return energy_data

def plot_tor(options,DIH_PARAM,ANG_IND,ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten ):
    import matplotlib.pyplot as plt

    
    # plt.ylabel('Energy (kcal/mol)')
    plt.ylabel('Energy (eV)')
    plt.xlabel('dihedral angle (deg)')
    
    plt.plot( qm_angle,qm_en_s,'kx', label="mp2" )
    plt.plot( qm_angle,ff_en_s,"g+", label="$V_{tor}$=0" )
    plt.plot( qm_angle,target_en_s,"r^", label="mp2 - ($V_{tor}$=0)" )
    plt.plot( qm_angle,fitted_ptor,"b-", label="$V_{tor}$ fit" )
    plt.plot( qm_angle, fit_toten,"b--", label="$ff_{en}$ fit")
    plt.legend(loc=(0.67,0.72),prop={'size':10})
    
    #plt.show()
    
    plt.savefig(options.plot_tor,format='pdf')


def plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp ):
    import matplotlib.pyplot as plt

    plt.cla()
    # plt.ylabel('Energy (kcal/mol)')
    plt.ylabel('Energy (eV)')
    plt.xlabel('dihedral angle (deg)')
    
    plt.plot( qm_angle,target_en_s,"r^", label="mp2 - ($V_{tor}$=0)" )
    plt.plot( qm_angle,fitted_ptor,"b-", label="$V_{tor}$ fit" )
    
    print FF_DIHTYPES
    
    lstyle_list = ["b--","r--","g--","y--"]
    
    print " size ",len(fitted_comp[0])
    for param_ind in range(len(fitted_comp[0])):
        comp_i = []
        
        # Atom_types = " %s %s %s %s " % (FF_DIHTYPES[param_ind][0],FF_DIHTYPES[param_ind][1],FF_DIHTYPES[param_ind][2],FF_DIHTYPES[param_ind][3])
        
        for comp_subindx in range(len(fitted_comp)):
            comp_i.append( fitted_comp[comp_subindx][param_ind] )
        
        #print " compent i ",comp_i
        line_style =     lstyle_list[param_ind]
        plt.plot( qm_angle , comp_i,line_style, label="$V( theta)_{"+str(param_ind)+"}$") #+Atom_types)
        
    plt.legend(loc=(0.67,0.72),prop={'size':10})
    
    #plt.show()
    
    plt.savefig(options.plot_tor_comp,format='pdf')



def print_tor(options,DIH_PARAM,ANG_IND,ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param ):
    
    debug = 0
    
    fitted_ptor = []
    fit_toten = []
    fitted_comp = []
    
    # Calculate new torsional potential based on previous t(0)
    #
    for cent_indx in range(len(target_en_s)):
        ang_val = ff_angles[cent_indx]
        ff_en = ff_en_s[cent_indx]      
        # Sum torsional energies for each component
        fit_en,comp_en  =  sum_tor_comp(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        fitted_comp.append(comp_en)
        fitted_ptor.append( fit_en )
        fit_toten.append( fit_en + ff_en )

    fit_en_max = max( fit_toten ) 
    fit_en_min = min( fit_toten )
    
    
    
    # Write data file 
    f_fit = open(options.out_ftor,"w")
    f_fit.write(" angle ; mp2 ; V_tor ; mp2-Vtor ; V_tor fit ; ff_en fit; fit en components ")
    fit_toten_s = []
    for cent_indx in range(len(target_en_s)):
        qm_x = qm_angle[cent_indx]        # angle
        qm_en = qm_en_s[cent_indx]        # mp2
        ff_en = ff_en_s[cent_indx]        # $V_{tor}$=0
        t_en = target_en_s[cent_indx]     # 
        fit_en = fitted_ptor[cent_indx]   # 
        fit_ten = fit_toten[cent_indx]    # 
        fit_ten_s = fit_ten -fit_en_min   #
        fit_toten_s.append( fit_ten_s )
        
        fit_comps = ""
        for comp_i in fitted_comp[cent_indx] :
            fit_comps = fit_comps + " " + str(comp_i)
        
        f_fit.write( "\n %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f  %s " % (qm_x,qm_en,ff_en,t_en, fit_en,fit_ten_s,fit_comps))        
    
    f_fit.close()

    return (fitted_ptor , fit_toten_s, fitted_comp )



def print_tor_v2(options,DIH_PARAM,ANG_IND,ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param ):
    import math
    
    debug = 0
    
    fitted_ptor = []
    fit_toten = []
    fitted_comp = []
    
    # Calculate new torsional potential based on previous t(0)
    #
    for cent_indx in range(len(target_en_s)):
        ang_val = ff_angles[cent_indx]
        ff_en = ff_en_s[cent_indx]      
        # Sum torsional energies for each component
        #fit_en,comp_en  =  sum_tor_comp(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        
        fit_en = 0.0
        comp_en = []
        
        for angle_indx in range( len(ang_val)):
            theta = math.radians( ang_val[angle_indx] )
            param_ind = ANG_IND[angle_indx]
            F_coef =  DIH_PARAM[param_ind]
            
            if(   tor_type == "four" ):
                V_tor =  V_four(F_coef,theta)
            elif( tor_type == "RB" ):
                V_tor =  V_RB(F_coef,theta)
                
            elif( tor_type == "harmonic" ):
                V_tor =  V_k(F_coef,theta)
                
            elif(   tor_type == "four_F2" ):
                V_tor =  V_fourF2(F_coef,theta)
                
            else:
                print " unknown torsional function sum_tor_comp ",tor_type
                sys.exit(" option set incorrectly ")
            comp_en.append(V_tor)
            fit_en += V_tor
            
        fitted_comp.append(comp_en)
        fitted_ptor.append( fit_en )
        fit_toten.append( fit_en + ff_en )

    fit_en_max = max( fit_toten ) 
    fit_en_min = min( fit_toten )
    
    
    
    # Write data file 
    f_fit = open(options.out_ftor,"w")
    f_fit.write(" angle ; mp2 ; V_tor ; mp2-Vtor ; V_tor fit (eV); ff_en fit (eV); fit en components (eV) ")
    fit_toten_s = []
    for cent_indx in range(len(target_en_s)):
        qm_x = qm_angle[cent_indx]        # angle
        qm_en = qm_en_s[cent_indx]        # mp2
        ff_en = ff_en_s[cent_indx]        # $V_{tor}$=0
        t_en = target_en_s[cent_indx]     # 
        fit_en = fitted_ptor[cent_indx]   # 
        fit_ten = fit_toten[cent_indx]    # 
        fit_ten_s = fit_ten -fit_en_min   #
        fit_toten_s.append( fit_ten_s )
        
        fit_comps = ""
        for comp_i in fitted_comp[cent_indx] :
            fit_comps = fit_comps + " " + str(comp_i)
        
        f_fit.write( "\n %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f  %s " % (qm_x,qm_en,ff_en,t_en, fit_en,fit_ten_s,fit_comps))        
    
    f_fit.close()

    return (fitted_ptor , fit_toten_s, fitted_comp )


def V_four(p,theta_rad):
    import math
    # theta_rad dihedral angle in radians
    
    F1 = p[0]
    F2 = p[1]
    F3 = p[2]
    F4 = p[3]
    
    a1 =  theta_rad*1.0
    a2 =  theta_rad*2.0
    a3 =  theta_rad*3.0
    a4 =  theta_rad*4.0

    Vtor_f1 =  0.5*( F1*( 1 + math.cos( a1  ))  + F2*(1 -   math.cos( a2 )) + F3*( 1+ math.cos( a3 )) + F4*( 1- math.cos( a4 )))

    debug = 0
    if( debug ):
        print theta_rad,F1,F2,F3,F4,Vtor_f1,F2*(1 -   math.cos( a2 )) ,  math.cos( a2 )

    return Vtor_f1

def V_fourF2(p,theta_rad):
    import math
    # theta_rad dihedral angle in radians
    
    F2 = p[0]    
    
    a2 =  theta_rad*2.0

    Vtor_f1 =  F2*(1 -   math.cos( a2 )) 

    debug = 0
    if( debug ):
        print theta_rad,F2,a2,Vtor_f1,F2*(1 -   math.cos( a2 )) ,  math.cos( a2 )

    return Vtor_f1

def V_RB(p,theta_rad):
    import math
    import sys
    # theta_rad dihedral angle in radians
    # Ryckaert-Bellemans 
    
    
    debug = 0
    Vtor_f1 = -1.0*( sum(p) )
    if( debug ):
        print " p ",p
        print " Vtor_f1 ",Vtor_f1
        
    nmax = 4
    for n in range(1,nmax+1):
        Vtor_f1 += p[n-1]*( math.cos( theta_rad ) )**float(n)
        
        print "p[n]*( math.cos( theta_rad ) )**float(n)",p[n], ( math.cos( theta_rad ) ), float(n)
        
    if( debug ):
        print theta_rad,math.cos( theta_rad )
        print p
        sys.exit("  debug V_RB ")

    return Vtor_f1


def V_k(p,theta_rad):
    import math
    # theta_rad dihedral angle in radians
    # Ryckaert-Bellemans 
    
    
    debug = 0
    mult = 2.0 #float( int( p[1] ))
    phase = math.pi # float(int(p[2]))*math.pi # phase is usually 0 or pi
    ko =  p[0] #/ 10000
    Vtor_f1 = ko*( 1 +  (math.cos( mult*theta_rad - phase ) )  )
    
    if( debug ):
        print theta_rad,math.cos( theta_rad )
        print p

    return Vtor_f1


def sum_tor_test(DIH_PARAM,ang_val,ANG_IND):
    import math
    
    debug = 1
    
    n_param = len(DIH_PARAM)
    
    tor_lijk = 0.0

    theta = math.radians( ang_val )
    param_ind = ANG_IND[0]
    F_coef = []
        
    for p_indx in range( param_ind + n_param):
        F_coef.append(  DIH_PARAM[p_indx] )
        
    if( debug ):
        print " calc pot ",theta,F_coef
        
    V_tor =  V_four(F_coef,theta)
    tor_lijk += V_tor
    if( debug ):
        print " V_tor ",V_tor
    
    debug = 0
    if( debug ):
        sys.exit("sum_tor")
    
    print "tor_lijk ",tor_lijk
    
    return tor_lijk

def sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param):
    import math, sys 
    
    debug = 0
    
    #tor_type = "four"  # Fourier parameters
    #n_param = 4 #len(DIH_PARAM)
    
    #print len(DIH_PARAM[0])
    #sys.exit(" DIH_PARAM ")
    
    
    tor_lijk = 0.0
    
    for angle_indx in range( len(ang_val)):
        theta = math.radians( ang_val[angle_indx] )
        param_ind = ANG_IND[angle_indx]
        F_coef = []
        
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):            
            F_coef.append(  DIH_PARAM[p_indx] )
        
        if( debug ):
            print " calc pot ",theta,param_ind,param_ind*n_param,n_param,F_coef
            
        if(   tor_type == "four" ):
            V_tor =  V_four(F_coef,theta)
        elif( tor_type == "RB" ):
            V_tor =  V_RB(F_coef,theta)
            
        elif( tor_type == "harmonic" ):
            V_tor =  V_k(F_coef,theta)
            
        elif(   tor_type == "four_F2" ):
            V_tor =  V_fourF2(F_coef,theta)
        else:
            print " unknown torsional function in sum_tor",tor_type
            sys.exit(" option set incorrectly ")
        
        tor_lijk += V_tor
        if( debug ):
            print " V_tor ",V_tor
    
    debug = 0
    if( debug ):
        print "tor_lijk ",tor_lijk
        sys.exit("sum_tor")
    
    return tor_lijk

def sum_tor_comp(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param):
    import math, sys 
    
    debug = 0
    
    #tor_type = "four"  # Fourier parameters
    #n_param = 4 #len(DIH_PARAM)
    
    #print len(DIH_PARAM[0])
    #sys.exit(" DIH_PARAM ")
    
    
    tor_lijk = 0.0
    tor_comp = []
    
    for angle_indx in range( len(ang_val)):
        theta = math.radians( ang_val[angle_indx] )
        param_ind = ANG_IND[angle_indx]
        F_coef = []
        
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            F_coef.append(  DIH_PARAM[p_indx] )
        
        if( debug ):
            print " calc pot ",theta,param_ind,param_ind*n_param,n_param,F_coef
            
        if(   tor_type == "four" ):
            V_tor =  V_four(F_coef,theta)
        elif( tor_type == "RB" ):
            V_tor =  V_RB(F_coef,theta)
            
        elif( tor_type == "harmonic" ):
            V_tor =  V_k(F_coef,theta)
            
        elif(   tor_type == "four_F2" ):
            V_tor =  V_fourF2(F_coef,theta)
            
        else:
            print " unknown torsional function sum_tor_comp ",tor_type
            sys.exit(" option set incorrectly ")
        tor_comp.append(V_tor)
        tor_lijk += V_tor
        if( debug ):
            print " V_tor ",V_tor
    
    debug = 0
    if( debug ):
        print "tor_lijk ",tor_lijk
        sys.exit("sum_tor")
    
    return tor_lijk, tor_comp 

def residuals(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param):
    
    debug = 0
 
    #n_param = 4 #len(DIH_PARAM)
    
    #print "wt_coef",wt_coef
    #wt_coef = 0.0
    #DIH_PARAM[0] = 1.0
    #DIH_PARAM[4] = 2.0
    #DIH_PARAM[8] = 3.0
    #DIH_PARAM[12] = 4.0
    #print "DIH_PARAM",DIH_PARAM
    
    
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
    
    debug = 0
    if( debug ):
        sys.exit("residuals")
    
    debug = 0
    NUMB_DIHTYPES = max( ANG_IND )
    if( debug ): print " NUMB_DIHTYPES ",NUMB_DIHTYPES

    for param_ind in range( NUMB_DIHTYPES ):
        F_coef = []
        sq_coef = 0.0
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            F_coef.append(  DIH_PARAM[p_indx] )
            
            
            sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
            
        if( debug ):
            print param_ind
            print wt_coef,F_coef[0],F_coef[1],F_coef[2],F_coef[3]
            print sq_coef

        resid.append(sq_coef)

    return resid
    
def residuals_v2(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param):
    
    debug = 0
 
    #n_param = 4 #len(DIH_PARAM)
    
    #print "wt_coef",wt_coef
    #wt_coef = 0.0
    #DIH_PARAM[0] = 1.0
    #DIH_PARAM[4] = 2.0
    #DIH_PARAM[8] = 3.0
    #DIH_PARAM[12] = 4.0
    #print "DIH_PARAM",DIH_PARAM
    
    
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
    
    debug = 0
    if( debug ):
        sys.exit("residuals")
    
    debug = 1
    NUMB_DIHTYPES = max( ANG_IND )
    if( debug ): print " NUMB_DIHTYPES ",NUMB_DIHTYPES

    for param_ind in range( NUMB_DIHTYPES ):
        F_coef = []
        sq_coef = 0.0
        p_cnt = -1
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            p_cnt += 1
            F_coef.append(  DIH_PARAM[p_indx] )
            
            sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
            
            for param_j in range( NUMB_DIHTYPES ):
                #for p_jndx in range( param_j*n_param,param_j*n_param + n_param):
                p_jndx = param_j*n_param + p_cnt
                delta_c = DIH_PARAM[p_indx] - DIH_PARAM[p_jndx]
                
                if( debug): print " delta coef ",DIH_PARAM[p_indx] , DIH_PARAM[p_jndx] , delta_c
                
                sq_coef += wt_coef*( delta_c*delta_c )
                
        if( debug ):
            print param_ind
            print wt_coef,F_coef[0],F_coef[1],F_coef[2],F_coef[3]
            print sq_coef

        resid.append(sq_coef)

    return resid
    
def residuals_v3(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param):
    
    debug = 0
    wt_coef = 0.1 
 
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(DIH_PARAM) ):
        p_eV = DIH_PARAM[param_i] #*EVTOKCAL
        #p_kcalmol_rnd = round(p_kcalmol,4)
        #DIH_PARAM[param_i] = p_kcalmol_rnd/EVTOKCAL
        
        DIH_PARAM[param_i]  = round(DIH_PARAM[param_i] ,4)
        if( debug ): print "   rounding ",p_eV," to ",DIH_PARAM[param_i]
        
        # Zero out parameters < 0.05 kcal/mol
        #if( DIH_PARAM[param_i]  < 0.002 ): DIH_PARAM[param_i]  = 0.0 
        
    
    
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
        
    #    
    #n_dihtypes = max( ANG_IND ) + 1
    #for param_ind in range( n_dihtypes ):
    #    sq_coef = 0.0
    #    p_cnt = -1
    #    
    #    
    #    
    #    if( debug ): print "    set index ",param_ind
    #    
    #    #
    #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
    #        p_cnt += 1
    #
    #        #if( debug ): print "    Coeff  ",p_indx,  DIH_PARAM[p_indx]
    #        
    #        #F_coef.append(  DIH_PARAM[p_indx] )
    #        
    #        # Control magnitude 
    #        sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
    #        
    return resid
    
def residuals_v4(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs):
    
    debug = 0
    wt_coef = 0.1
    # Target en is not real potential energy surface!!! idiot god!
    wt_min = 0.00
    wt_transL = 100.0   # Low energy trans < max_temp
    wt_transH = 0.01   # High energy trans > max_temp
     
    # For transitions under 400 K increase weight
    #   as they are possible during simulations 
    KTOEV = 0.0000861738
    max_temp = 400  # K
    max_temp_eV  = max_temp*KTOEV
 
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(DIH_PARAM) ):
        p_eV = DIH_PARAM[param_i] #*EVTOKCAL
        #p_kcalmol_rnd = round(p_kcalmol,4)
        #DIH_PARAM[param_i] = p_kcalmol_rnd/EVTOKCAL
        
        DIH_PARAM[param_i]  = round(DIH_PARAM[param_i] ,4)
        if( debug ): print "   rounding ",p_eV," to ",DIH_PARAM[param_i]
        
        # Zero out parameters < 0.05 kcal/mol
        #if( DIH_PARAM[param_i]  < 0.002 ): DIH_PARAM[param_i]  = 0.0 
        
    
    
    resid = []
    en_fit = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        en_fit.append(tor_en)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
    #
    #
    ## Make sure each minimum is actually a minimum 
    #for indx in range( len(min_indx) ):
    #    cent_indx =  min_indx[indx]
    #    
    #    min_error = 0 
    #    
    #    
    #    if( cent_indx >  0 ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx-1] ):  min_error += wt_min
    #    if( cent_indx <    len(targets)  ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx+1] ):  min_error += wt_min
    #            
    #    resid.append(min_error)
    #    

        
    for indx in range( len(trans_list) ):
        
        max_i =  trans_indxs[indx][0]
        min_i =  trans_indxs[indx][1]
        
        trans_val = en_fit[max_i] - en_fit[min_i]
        delta_trans = trans_list[indx] - trans_val
        
        if(  trans_val < max_temp_eV ):
            trans_error = wt_transL*delta_trans*delta_trans
        else:
            trans_error = wt_transH*delta_trans*delta_trans
            
        resid.append(trans_error)
                
    #    
    #n_dihtypes = max( ANG_IND ) + 1
    #for param_ind in range( n_dihtypes ):
    #    sq_coef = 0.0
    #    p_cnt = -1
    #    
    #    
    #    
    #    if( debug ): print "    set index ",param_ind
    #    
    #    #
    #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
    #        p_cnt += 1
    #
    #        #if( debug ): print "    Coeff  ",p_indx,  DIH_PARAM[p_indx]
    #        
    #        #F_coef.append(  DIH_PARAM[p_indx] )
    #        
    #        # Control magnitude 
    #        sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
    #        
    return resid

def residuals_initial(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs):
    import math, sys 
    
    debug = 0
    wt_coef = 0.00
    # if two of the same parameters for different sets of coeffcients get too far appart 
    wt_deltaC = 0.0  # Off 
    maxdcoef  = 0.5  # eV
    maxdcoef_sq = maxdcoef*maxdcoef
    # wieght for energy controbutions from each dihedral type at each maximum
    wt_en_comp = 0.00  # Off 
    # Target en is not real potential energy surface!!! idiot god!
    wt_min = 0.00
    wt_transL = 100.0   # Low energy trans < max_temp
    wt_transH = 0.01   # High energy trans > max_temp
     
    # For transitions under 400 K increase weight
    #   as they are possible during simulations 
    KTOEV = 0.0000861738
    max_temp = 400  # K
    max_temp_eV  = max_temp*KTOEV
 
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(DIH_PARAM) ):
        p_eV = DIH_PARAM[param_i] #*EVTOKCAL
        #p_kcalmol_rnd = round(p_kcalmol,4)
        #DIH_PARAM[param_i] = p_kcalmol_rnd/EVTOKCAL
        
        DIH_PARAM[param_i]  = round(DIH_PARAM[param_i] ,4)
        if( debug ): print "   rounding ",p_eV," to ",DIH_PARAM[param_i]
        
        # Zero out parameters < 0.05 kcal/mol
        #if( DIH_PARAM[param_i]  < 0.002 ): DIH_PARAM[param_i]  = 0.0 
        
    
    
    resid = []
    en_fit = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        # tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        
        ang_val = ff_angles[cent_indx][0]
        theta = math.radians( ang_val )
        
        F_coef = DIH_PARAM
        if(   tor_type == "four" ):
            V_tor =  V_four(F_coef,theta)
        elif( tor_type == "RB" ):
            V_tor =  V_RB(F_coef,theta)
            
        elif( tor_type == "harmonic" ):
            V_tor =  V_k(F_coef,theta)
            
        elif(   tor_type == "four_F2" ):
            V_tor =  V_fourF2(F_coef,theta)
        else:
            print " unknown torsional function in sum_tor",tor_type
            sys.exit(" option set incorrectly ")
        
        tor_en = V_tor        
        
        en_fit.append(tor_en)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
        
        

    for indx in range( len(trans_list) ):
        
        max_i =  trans_indxs[indx][0]
        min_i =  trans_indxs[indx][1]
        
        trans_val = en_fit[max_i] - en_fit[min_i]
        delta_trans = trans_list[indx] - trans_val
        
        if(  trans_val < max_temp_eV ):
            trans_error = wt_transL*delta_trans*delta_trans
        else:
            trans_error = wt_transH*delta_trans*delta_trans
            
        resid.append(trans_error)
    
    return resid


def residuals_initialv2(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs):
    import math, sys 
    
    debug = 0
    
    wt_coef = 0.00
    # if two of the same parameters for different sets of coeffcients get too far appart 
    wt_deltaC = 0.0  # Off 
    maxdcoef  = 0.5  # eV
    maxdcoef_sq = maxdcoef*maxdcoef
    # wieght for energy controbutions from each dihedral type at each maximum
    wt_en_comp = 0.00  # Off 
    # Target en is not real potential energy surface!!! idiot god!
    wt_min = 0.00
    wt_transL = 100.0   # Low energy trans < max_temp
    wt_transH = 0.01   # High energy trans > max_temp
     
    # For transitions under 400 K increase weight
    #   as they are possible during simulations 
    KTOEV = 0.0000861738
    max_temp = 400  # K
    max_temp_eV  = max_temp*KTOEV
 
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(DIH_PARAM) ):
        p_eV = DIH_PARAM[param_i] #*EVTOKCAL
        #p_kcalmol_rnd = round(p_kcalmol,4)
        #DIH_PARAM[param_i] = p_kcalmol_rnd/EVTOKCAL
        
        DIH_PARAM[param_i]  = round(DIH_PARAM[param_i] ,4)
        if( debug ): print "   rounding ",p_eV," to ",DIH_PARAM[param_i]
        
        # Zero out parameters < 0.05 kcal/mol
        #if( DIH_PARAM[param_i]  < 0.002 ): DIH_PARAM[param_i]  = 0.0 
        
    
    
    resid = []
    en_fit = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        # tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        
        ang_val = ff_angles[cent_indx]
        theta = math.radians( ang_val )
        
        F_coef = DIH_PARAM
        if(   tor_type == "four" ):
            V_tor =  V_four(F_coef,theta)
        elif( tor_type == "RB" ):
            V_tor =  V_RB(F_coef,theta)
            
        elif( tor_type == "harmonic" ):
            V_tor =  V_k(F_coef,theta)
            
        elif(   tor_type == "four_F2" ):
            V_tor =  V_fourF2(F_coef,theta)
        else:
            print " unknown torsional function in sum_tor",tor_type
            sys.exit(" option set incorrectly ")
        
        tor_en = V_tor        
        
        en_fit.append(tor_en)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
        
        

    for indx in range( len(trans_list) ):
        
        max_i =  trans_indxs[indx][0]
        min_i =  trans_indxs[indx][1]
        
        trans_val = en_fit[max_i] - en_fit[min_i]
        delta_trans = trans_list[indx] - trans_val
        
        if(  trans_val < max_temp_eV ):
            trans_error = wt_transL*delta_trans*delta_trans
        else:
            trans_error = wt_transH*delta_trans*delta_trans
            
        resid.append(trans_error)
    
    return resid

    
def residuals_v5(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs):
    import math, sys 
    
    debug = 0
    wt_coef = 0.00
    # if two of the same parameters for different sets of coeffcients get too far appart 
    wt_deltaC = 0.0   # Off 
    maxdcoef  = 0.5  # eV
    maxdcoef_sq = maxdcoef*maxdcoef
    # wieght for energy controbutions from each dihedral type at each maximum
    wt_en_comp = 0.00  # Off 
    # Target en is not real potential energy surface!!! idiot god!
    wt_min = 0.00
    wt_transL = 100.0   # Low energy trans < max_temp
    wt_transH = 0.01   # High energy trans > max_temp
     
    # For transitions under 400 K increase weight
    #   as they are possible during simulations 
    KTOEV = 0.0000861738
    max_temp = 400  # K
    max_temp_eV  = max_temp*KTOEV
    
 
    resid = []
    # Constrain parameters under max transion energy 
    max_trans = max( trans_list )
    max_trans_sq = max_trans*max_trans
    wt_maxtrans = 0.0
    
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(DIH_PARAM) ):
        p_eV = DIH_PARAM[param_i] #*EVTOKCAL
        #p_kcalmol_rnd = round(p_kcalmol,4)
        #DIH_PARAM[param_i] = p_kcalmol_rnd/EVTOKCAL
        
        DIH_PARAM[param_i]  = round(DIH_PARAM[param_i] ,4)
        if( debug ): print "   rounding ",p_eV," to ",DIH_PARAM[param_i]
        
        # Zero out parameters < 0.05 kcal/mol
        #if( DIH_PARAM[param_i]  < 0.002 ): DIH_PARAM[param_i]  = 0.0
        
        
        # Constrain parameters under max transion energy
        coef_sq = DIH_PARAM[param_i]*DIH_PARAM[param_i] 
        # if( coef_sq > max_trans_sq ): resid.append( wt_maxtrans )
    
    if( tor_type == "RB" ):
        # Set C0 = -( C1 + C2 + C3 + C4 )
        n_dihtypes = max( ANG_IND ) + 1
        for param_ind in range( n_dihtypes ):
            F_coef = []            
            for p_indx in range( param_ind*n_param + 1,param_ind*n_param + n_param):
                F_coef.append(  DIH_PARAM[p_indx] )
            DIH_PARAM[param_ind*n_param ] = -1.0*sum( F_coef )
            #if( DIH_PARAM[param_ind*n_param ] != -1.0*sum( F_coef ) ): resid.append( wt_maxtrans )
                
    
    en_fit = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        # tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        
        
        tor_anglei = []
        tor_en = 0.0
        
        for angle_indx in range( len(ang_val)):
            theta = math.radians( ang_val[angle_indx] )
            param_ind = ANG_IND[angle_indx]
            F_coef = []
            
            
            for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):            
                F_coef.append(  DIH_PARAM[p_indx] )
            
            if( debug ):
                print " calc pot ",theta,param_ind,param_ind*n_param,n_param,F_coef
                
            if(   tor_type == "four" ):
                V_tor =  V_four(F_coef,theta)
            elif( tor_type == "RB" ):
                
                V_tor =  V_RB(F_coef,theta)
                
            elif( tor_type == "harmonic" ):
                V_tor =  V_k(F_coef,theta)
                
            elif(   tor_type == "four_F2" ):
                V_tor =  V_fourF2(F_coef,theta)
            else:
                print " unknown torsional function in sum_tor",tor_type
                sys.exit(" option set incorrectly ")
            
            tor_en += V_tor
            tor_anglei.append( V_tor )
        
        
        en_fit.append(tor_en)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
        
        
    #
    #
    ## Make sure each minimum is actually a minimum 
    #for indx in range( len(min_indx) ):
    #    cent_indx =  min_indx[indx]
    #    
    #    min_error = 0 
    #    
    #    
    #    if( cent_indx >  0 ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx-1] ):  min_error += wt_min
    #    if( cent_indx <    len(targets)  ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx+1] ):  min_error += wt_min
    #            
    #    resid.append(min_error)
    #    

    n_angles = len(ang_val)
    
    
    for indx in range( len(max_indx) ):
        cent_indx =  max_indx[indx]
        ang_val = ff_angles[cent_indx]
        
        t_en = targets[cent_indx]
        t_en_norm = float(n_angles)
        
        for angle_indx in range( len(ang_val)):
            theta = math.radians( ang_val[angle_indx] )
            param_ind = ANG_IND[angle_indx]
            F_coef = []
            
            
            for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):            
                F_coef.append(  DIH_PARAM[p_indx] )
            
            if(   tor_type == "four" ):
                V_tor =  V_four(F_coef,theta)
            elif( tor_type == "RB" ):
                V_tor =  V_RB(F_coef,theta)
                
            elif( tor_type == "harmonic" ):
                V_tor =  V_k(F_coef,theta)
                
            elif(   tor_type == "four_F2" ):
                V_tor =  V_fourF2(F_coef,theta)
            else:
                print " unknown torsional function in sum_tor",tor_type
                sys.exit(" option set incorrectly ")
            
            
            t_en_comp_error = t_en_norm - V_tor
            
            en_comp_sq = wt_en_comp*t_en_comp_error*t_en_comp_error
            
            
            resid.append(en_comp_sq)

        
    for indx in range( len(trans_list) ):
        
        max_i =  trans_indxs[indx][0]
        min_i =  trans_indxs[indx][1]
        
        trans_val = en_fit[max_i] - en_fit[min_i]
        delta_trans = trans_list[indx] - trans_val
        
        if(  trans_val < max_temp_eV ):
            trans_error = wt_transL*delta_trans*delta_trans
        else:
            trans_error = wt_transH*delta_trans*delta_trans
            
        resid.append(trans_error)
                
        
    n_dihtypes = max( ANG_IND ) + 1
    for param_ind in range( n_dihtypes ):
        F_coef = []
        sq_coef = 0.0
        p_cnt = -1
        
        
        
        if( debug ): print "    set index ",param_ind
        
        #
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            p_cnt += 1
    
            #if( debug ): print "    Coeff  ",p_indx,  DIH_PARAM[p_indx]
            
            F_coef.append(  DIH_PARAM[p_indx] )
            
            # Control magnitude 
            sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
            
        if( param_ind > 0 ):
            #print "    len(F_coef)  ",len(F_coef)
            #print "    len(F_coef_prev)  ",len(F_coef_prev)
            for f_indx in range( len(F_coef) ):
                delta_coef = F_coef_prev[f_indx] - F_coef[f_indx]
                
                if( debug ): print "   delta_coef   ",F_coef[f_indx],F_coef_prev[f_indx], delta_coef
                dcoef_sq = delta_coef*delta_coef
                
                # if delta is > 5 kcal/mol .2 eV place error 
                if( dcoef_sq > maxdcoef_sq ):
                    sq_coef +=  wt_deltaC
                
        # Save parameters to compare to next set 
        F_coef_prev = []
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            F_coef_prev.append(  DIH_PARAM[p_indx] )
                
                
        resid.append(sq_coef)
            
    return resid


def calc_rms(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs):
    import numpy 
    
    print " Calcualting RMSD "
    
    debug = 0
 
    # For transitions under 400 K increase weight
    #   as they are possible during simulations 
    KTOEV = 0.0000861738
    max_temp = 400  # K
    max_temp_eV  = max_temp*KTOEV
    
    deltasq_en = []
    en_fit = []
    t_cnt = 0
    for cent_indx in range(len(targets)):
        t_cnt += 1
        t_en = targets[cent_indx]
        ang_val = ff_angles[cent_indx]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        en_fit.append(tor_en)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta =  delta_en*delta_en
        
        deltasq_en.append(sq_delta)
        
        
        # print ff_angles[cent_indx][0],tor_en*1000," meV ",sq_delta

    print sum(deltasq_en),len(targets),t_cnt
    
    rmsd_en = numpy.sqrt( sum(deltasq_en)/len(targets) )
    
    print "  RMSD in energy ",rmsd_en
    #
    # Target en is not real potential energy surface!!! idiot god!
    ## Make sure each minimum is actually a minimum 
    #for indx in range( len(min_indx) ):
    #    cent_indx =  min_indx[indx]
    #    
    #    min_error = 0 
    #    
    #    
    #    if( cent_indx >  0 ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx-1] ):  min_error = 1 
    #    if( cent_indx <    len(targets)  ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx+1] ):  min_error = 1 
    #            
    #    if( min_error ):
    #        print " minimum ",indx," is not true minimum !!! "
    #        print " Tragets meV ",targets[cent_indx-1]*1000 , targets[cent_indx]*1000 ,  targets[cent_indx+1]*1000
    #        print " FItted  meV ",en_fit[cent_indx-1]*1000 , en_fit[cent_indx]*1000 ,  en_fit[cent_indx+1]*1000

    deltasq_trans = 0.0 
    for indx in range( len(trans_list) ):
        
        max_i =  trans_indxs[indx][0]
        min_i =  trans_indxs[indx][1]
        
        trans_val = en_fit[max_i] - en_fit[min_i]
        delta_trans = trans_list[indx] - trans_val
        deltasq_trans += delta_trans*delta_trans
        
        print  " Transition ",indx,"  Traget ",trans_list[indx] *1000," meV ",trans_list[indx]/KTOEV," K fitted ", trans_val*1000," meV ", trans_val/KTOEV," K "
        
    rmsd_trans = numpy.sqrt( deltasq_trans/len(trans_list) )
    print "  RMSD in transition energy ",rmsd_trans
    
    #    
    #n_dihtypes = max( ANG_IND ) + 1
    #for param_ind in range( n_dihtypes ):
    #    sq_coef = 0.0
    #    p_cnt = -1
    #    
    #    
    #    
    #    if( debug ): print "    set index ",param_ind
    #    
    #    #
    #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
    #        p_cnt += 1
    #
    #        #if( debug ): print "    Coeff  ",p_indx,  DIH_PARAM[p_indx]
    #        
    #        #F_coef.append(  DIH_PARAM[p_indx] )
    #        
    #        # Control magnitude 
    #        sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
    #        
    return rmsd_en,rmsd_trans
    
def calc_rms_v2(  param_list, fit_param, DIH_PARAM,  ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs):
    import numpy 
    import math, sys 
    
    print " Calcualting RMSD "
    
    debug = 1
    if( debug ):
        print " calc_rms_v2  is in debug mode  "
        
        
    # For transitions under 400 K increase weight
    #   as they are possible during simulations 
    KTOEV = 0.0000861738
    max_temp = 400  # K
    max_temp_eV  = max_temp*KTOEV
    
    deltasq_en = []
    en_fit = []
    t_cnt = 0
    for cent_indx in range(len(targets)):
        t_cnt += 1
        t_en = targets[cent_indx]
        ang_val = ff_angles[cent_indx]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        #tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        
        tor_en = 0.0

        p_cnt = -1        
        for angle_indx in range( len(ang_val) ): 
            theta = math.radians( ang_val[angle_indx] )
            param_ind = ANG_IND[angle_indx]
            F_coef = DIH_PARAM[param_ind]

            if( debug ):
                print "  angle_indx ",angle_indx,param_ind,F_coef
                
            if(fit_param[param_ind] == 1 ):
                p_cnt += 1
                F_coef = []
                if( debug ):
                    print "  For index ",param_ind,"  will use fitted parameters in param_list "
                for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                    #print " append biaryl_param ",p_indx
                    F_coef.append(  param_list[p_indx] )
            
            # V_tor =  V_four(F_coef,theta)
            tor_en += V_four(F_coef,theta)

        en_fit.append(tor_en)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta =  delta_en*delta_en
        
        deltasq_en.append(sq_delta)        
        
        # print ff_angles[cent_indx][0],tor_en*1000," meV ",sq_delta

    print sum(deltasq_en),len(targets),t_cnt
    
    rmsd_en = numpy.sqrt( sum(deltasq_en)/len(targets) )
    
    print "  RMSD in energy ",rmsd_en
    #
    # Target en is not real potential energy surface!!! idiot god!
    ## Make sure each minimum is actually a minimum 
    #for indx in range( len(min_indx) ):
    #    cent_indx =  min_indx[indx]
    #    
    #    min_error = 0 
    #    
    #    
    #    if( cent_indx >  0 ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx-1] ):  min_error = 1 
    #    if( cent_indx <    len(targets)  ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx+1] ):  min_error = 1 
    #            
    #    if( min_error ):
    #        print " minimum ",indx," is not true minimum !!! "
    #        print " Tragets meV ",targets[cent_indx-1]*1000 , targets[cent_indx]*1000 ,  targets[cent_indx+1]*1000
    #        print " FItted  meV ",en_fit[cent_indx-1]*1000 , en_fit[cent_indx]*1000 ,  en_fit[cent_indx+1]*1000

    deltasq_trans = 0.0 
    for indx in range( len(trans_list) ):
        
        max_i =  trans_indxs[indx][0]
        min_i =  trans_indxs[indx][1]
        
        trans_val = en_fit[max_i] - en_fit[min_i]
        delta_trans = trans_list[indx] - trans_val
        deltasq_trans += delta_trans*delta_trans
        
        print  " Transition ",indx,"  Traget ",trans_list[indx] *1000," meV ",trans_list[indx]/KTOEV," K fitted ", trans_val*1000," meV ", trans_val/KTOEV," K "
        
    rmsd_trans = numpy.sqrt( deltasq_trans/len(trans_list) )
    print "  RMSD in transition energy ",rmsd_trans
    
    #    
    #n_dihtypes = max( ANG_IND ) + 1
    #for param_ind in range( n_dihtypes ):
    #    sq_coef = 0.0
    #    p_cnt = -1
    #    
    #    
    #    
    #    if( debug ): print "    set index ",param_ind
    #    
    #    #
    #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
    #        p_cnt += 1
    #
    #        #if( debug ): print "    Coeff  ",p_indx,  DIH_PARAM[p_indx]
    #        
    #        #F_coef.append(  DIH_PARAM[p_indx] )
    #        
    #        # Control magnitude 
    #        sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
    #        
    return rmsd_en,rmsd_trans
    
    
def residuals_v6( param_list, fit_param, DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs):
    import math, sys 
    
    debug = 0
    wt_coef = 0.00
    # if two of the same parameters for different sets of coeffcients get too far appart 
    wt_deltaC = 0.0   # Off 
    maxdcoef  = 0.5  # eV
    maxdcoef_sq = maxdcoef*maxdcoef
    # wieght for energy controbutions from each dihedral type at each maximum
    wt_en_comp = 0.00  # Off 
    # Target en is not real potential energy surface!!! idiot god!
    wt_min = 0.00
    wt_transL = 10000.0   # Low energy trans < max_temp
    wt_transH = 10000.0   # High energy trans > max_temp
     
    # For transitions under 400 K increase weight
    #   as they are possible during simulations 
    KTOEV = 0.0000861738
    max_temp = 400  # K
    max_temp_eV  = max_temp*KTOEV
    
 
    resid = []
    # Constrain parameters under max transion energy 
    max_trans = max( trans_list )
    max_trans_sq = max_trans*max_trans
    wt_maxtrans = 0.0
    
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(param_list) ):
        param_list[param_i]  = round(param_list[param_i] ,4)
        
    
    en_fit = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        # tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        
        
        tor_anglei = []
        tor_en = 0.0

        p_cnt = -1        
        for angle_indx in range( len(ang_val)):
            theta = math.radians( ang_val[angle_indx] )
            param_ind = ANG_IND[angle_indx]
            F_coef = DIH_PARAM[param_ind]

            if(fit_param[param_ind] == 1 ):
                p_cnt += 1
                F_coef = []
                if( debug ):
                    print "  For index ",param_ind,"  will use fitted parameters in param_list "
                for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                    #print " append biaryl_param ",p_indx
                    F_coef.append(  param_list[p_indx] )
            
            V_tor =  V_four(F_coef,theta)
            tor_en += V_tor


            tor_anglei.append( V_tor )
        
        
        en_fit.append(tor_en)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
        
        
    #
    #
    ## Make sure each minimum is actually a minimum 
    #for indx in range( len(min_indx) ):
    #    cent_indx =  min_indx[indx]
    #    
    #    min_error = 0 
    #    
    #    
    #    if( cent_indx >  0 ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx-1] ):  min_error += wt_min
    #    if( cent_indx <    len(targets)  ):
    #        if( en_fit[cent_indx] >  en_fit[cent_indx+1] ):  min_error += wt_min
    #            
    #    resid.append(min_error)
    #    

    n_angles = len(ang_val)
    
    for indx in range( len(max_indx) ):
        cent_indx =  max_indx[indx]
        ang_val = ff_angles[cent_indx]
        
        t_en = targets[cent_indx]
        t_en_norm = float(n_angles)
        
        tor_anglei = []
        tor_en = 0.0

        p_cnt = -1        
        for angle_indx in range( len(ang_val)):
            theta = math.radians( ang_val[angle_indx] )
            param_ind = ANG_IND[angle_indx]
            F_coef = DIH_PARAM[param_ind]

            if(fit_param[param_ind] == 1 ):
                p_cnt += 1
                F_coef = []
                if( debug ):
                    print "  For index ",param_ind,"  will use fitted parameters in param_list "
                for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                    #print " append biaryl_param ",p_indx
                    F_coef.append(  param_list[p_indx] )
            
            V_tor =  V_four(F_coef,theta)
            
            tor_en += V_tor 
            
            t_en_comp_error = t_en_norm - V_tor
            
            en_comp_sq = wt_en_comp*t_en_comp_error*t_en_comp_error
            
            resid.append(en_comp_sq)

        
    for indx in range( len(trans_list) ):
        
        max_i =  trans_indxs[indx][0]
        min_i =  trans_indxs[indx][1]
        
        trans_val = en_fit[max_i] - en_fit[min_i]
        delta_trans = trans_list[indx] - trans_val
        
        if(  trans_val < max_temp_eV ):
            trans_error = wt_transL*delta_trans*delta_trans
        else:
            trans_error = wt_transH*delta_trans*delta_trans
            
        resid.append(trans_error)
                
        
    return resid
    
    
def residuals_kcos(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param):
    import math, sys
    
    wt_kpos = 10000.0 
    
    debug = 0
 
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
    
    debug = 0
    if( debug ):
        sys.exit("residuals")
    
    debug = 0
    n_dihtypes = max( ANG_IND )
    if( debug ): print " n_dihtypes ",n_dihtypes

    for param_ind in range( n_dihtypes ):
        F_coef = []
        sq_coef = 0.0
        p_cnt = -1
        #
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            p_cnt += 1
            
            F_coef.append(  DIH_PARAM[p_indx] )
            
            sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
            
            if(  DIH_PARAM[p_indx] < 0.0  ): sq_coef +=  wt_kpos
            #
            #for param_j in range( n_dihtypes ):
            #    #for p_jndx in range( param_j*n_param,param_j*n_param + n_param):
            #    p_jndx = param_j*n_param + p_cnt
            #    delta_c = DIH_PARAM[p_indx] - DIH_PARAM[p_jndx]
            #    
            #    if( debug): print " delta coef ",DIH_PARAM[p_indx] , DIH_PARAM[p_jndx] , delta_c
            #    
            #    sq_coef += wt_coef*( delta_c*delta_c )
            #
        
        #print "  Harmonic coefficents "
        #print "                        k_o ",F_coef[0]
        #print "                        mult ",F_coef[1]
        #print "                        phase ",F_coef[2] #*math.pi

        # Keep mult and phase integers
        wt_dint = 1.0 
        dint =F_coef[1] - int(F_coef[1] )
        sq_coef += wt_dint*dint*dint
        dint =F_coef[2] - int(F_coef[2] )
        sq_coef += wt_dint*dint*dint
        
        # Keep k +
        # wt_kpos =10000.0 
        # if(F_coef[0] < 0  ): sq_coef +=  wt_kpos
            
        if( debug ):
            print param_ind
            print wt_coef,F_coef[0],F_coef[1],F_coef[2] #,F_coef[3]
            print sq_coef

        resid.append(sq_coef)
        
    #sys.exit(" coef check ")
    
    return resid

def residuals_fourF2(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param):
    import math, sys
    global EVTOKCAL
    
    wt_kpos = 10000.0
    wt_deltaC = 100.0 
    wt_coef = 1.0
    
    debug = 0
    
    
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(DIH_PARAM) ):
        p_eV = DIH_PARAM[param_i] #*EVTOKCAL
        #p_kcalmol_rnd = round(p_kcalmol,4)
        #DIH_PARAM[param_i] = p_kcalmol_rnd/EVTOKCAL
        
        DIH_PARAM[param_i]  = round(DIH_PARAM[param_i] ,4)
        if( debug ): print "   rounding ",p_eV," to ",DIH_PARAM[param_i]
        
        # Zero out parameters < 0.05 kcal/mol
        if( DIH_PARAM[param_i]  < 0.002 ): DIH_PARAM[param_i]  = 0.0 
        
    
 
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
    
    debug = 0
    if( debug ):
        sys.exit("residuals")
    
    debug = 0
    n_dihtypes = max( ANG_IND ) + 1
    if( debug ): print " n_dihtypes ",n_dihtypes

    for param_ind in range( n_dihtypes ):
        F_coef = []
        sq_coef = 0.0
        p_cnt = -1
        
        
        
        if( debug ): print "    set index ",param_ind
        
        #
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            p_cnt += 1

            if( debug ): print "    Coeff  ",p_indx,  DIH_PARAM[p_indx]
            
            F_coef.append(  DIH_PARAM[p_indx] )
            
            # Control magnitude 
            sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
            
            # Keep positive 
            if(  DIH_PARAM[p_indx] < 0.0  ): sq_coef +=  wt_kpos
            #
        
        if( param_ind > 0 ):
            #print "    len(F_coef)  ",len(F_coef)
            #print "    len(F_coef_prev)  ",len(F_coef_prev)
            for f_indx in range( len(F_coef) ):
                delta_coef = F_coef_prev[f_indx] - F_coef[f_indx]
                
                if( debug ): print "   delta_coef   ",F_coef[f_indx],F_coef_prev[f_indx], delta_coef
                sq_coef +=  wt_deltaC*delta_coef*delta_coef
                
        # Save parameters to compare to next set 
        F_coef_prev = []
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            F_coef_prev.append(  DIH_PARAM[p_indx] )
                
            
        if( debug ):
            print param_ind
            print wt_coef,F_coef[0] #,F_coef[1],F_coef[2] #,F_coef[3]
            print sq_coef

        resid.append(sq_coef)
        
    #sys.exit(" coef check ")
    
    return resid

def residuals_four(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef,tor_type,n_param):
    import math, sys
    global EVTOKCAL
    
    wt_kpos = 10000.0
    wt_deltaC = 0.0 
    wt_coef = 0.0
    
    
    debug = 0
    
    
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(DIH_PARAM) ):
        p_eV = DIH_PARAM[param_i] #*EVTOKCAL
        #p_kcalmol_rnd = round(p_kcalmol,4)
        #DIH_PARAM[param_i] = p_kcalmol_rnd/EVTOKCAL
        
        DIH_PARAM[param_i]  = round(DIH_PARAM[param_i] ,4)
        if( debug ): print "   rounding ",p_eV," to ",DIH_PARAM[param_i]
        
        # Zero out parameters < 0.05 kcal/mol
        if( DIH_PARAM[param_i]  < 0.002 ): DIH_PARAM[param_i]  = 0.0 
        
    
 
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        #ang_val = ff_angles[cent_indx]
        # hack !!!
        ang_val = ff_angles[cent_indx] #[0]
        wt = wt_angle[cent_indx]
        # Sum torsional energies for each component 
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND,tor_type,n_param)
        # Calculate error
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        #print wt,t_en,tor_en ,sq_delta
        resid.append(sq_delta)
    #
    #debug = 0
    #if( debug ):
    #    sys.exit("residuals")
    #
    #debug = 0
    #n_dihtypes = max( ANG_IND ) + 1
    #if( debug ): print " n_dihtypes ",n_dihtypes
    #
    #for param_ind in range( n_dihtypes ):
    #    F_coef = []
    #    sq_coef = 0.0
    #    p_cnt = -1
    #    
    #    
    #    
    #    if( debug ): print "    set index ",param_ind
    #    
    #    #
    #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
    #        p_cnt += 1
    #
    #        if( debug ): print "    Coeff  ",p_indx,  DIH_PARAM[p_indx]
    #        
    #        F_coef.append(  DIH_PARAM[p_indx] )
    #        
    #        # Control magnitude 
    #        sq_coef += wt_coef *( DIH_PARAM[p_indx] *DIH_PARAM[p_indx]  )
    #        
    #        # Keep positive 
    #        if(  DIH_PARAM[p_indx] < 0.0  ): sq_coef +=  wt_kpos
    #        #
    #    
    #    if( param_ind > 0 ):
    #        #print "    len(F_coef)  ",len(F_coef)
    #        #print "    len(F_coef_prev)  ",len(F_coef_prev)
    #        for f_indx in range( len(F_coef) ):
    #            delta_coef = F_coef_prev[f_indx] - F_coef[f_indx]
    #            
    #            if( debug ): print "   delta_coef   ",F_coef[f_indx],F_coef_prev[f_indx], delta_coef
    #            sq_coef +=  wt_deltaC*delta_coef*delta_coef
    #            
    #    # Save parameters to compare to next set 
    #    F_coef_prev = []
    #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
    #        F_coef_prev.append(  DIH_PARAM[p_indx] )
    #            
    #        
    #    if( debug ):
    #        print param_ind
    #        print wt_coef,F_coef[0] #,F_coef[1],F_coef[2] #,F_coef[3]
    #        print sq_coef
    #
    #    resid.append(sq_coef)
    #    
    #sys.exit(" coef check ")
    
    return resid

def tor_fit(DIH_PARAM,options):
    import sys
    from scipy import optimize
    
    resid = tor_resid(DIH_PARAM)
    error = 0.0 
    for er in resid:
        error = error + er*er
    print " initial error ",error


    DIH_PARAM,success = optimize.leastsq(tor_resid, DIH_PARAM, epsfcn=0.0001)

    resid = tor_resid(DIH_PARAM)
    error = 0.0 
    for er in resid:
        error = error + er*er
    print " final error ",error
    
    return DIH_PARAM

def read_dihtypes(struct_dir,job_name,options):
    import sys 
    import file_io
    
    debug = 1
    # Set initial paratemers

    tor_type = "four"
    n_param = 4
    
    param_o = []
    for n in range(n_param):
        param_o.append(0.0)
    #
    ## Set initial paratemers
    #
    if( tor_type == "four" ):
        
        param_o[1] = options.tor_paramo
        
    elif( tor_type == "RB" ):
        
        param_o[0] = options.tor_paramo
        param_o[2] = -1.0*options.tor_paramo
        
    elif( tor_type == "harmonic" ):
        
        param_o[0] = options.tor_paramo
        param_o[1] =  1.0
        param_o[2] =  1.0
        
    else:
        
        print " unknow torsional function ",tor_type
        sys.exit(" option set incorrectly ")
        
        
        
        

    ff_file =  options.ff_tor
    
    if ( file_io.file_exists(ff_file) ): 
        F = open(ff_file,'r')
        Lines = F.readlines()
        F.close()
    else:
        sys.exit(' data file does not exist')

    # Read in all components of fitted dihedral 
    dih_comp_name = struct_dir +'/' + job_name + "_fit.list"
    f_fit = open(dih_comp_name,'r')
    Lines_fit = f_fit.readlines()
    f_fit.close()


    print Lines_fit

    ang_types = []    
    for line_fit in Lines_fit:
        col_fit = line_fit.split()
        if ( col_fit[0] == "dihedral" ):
            ang_types.append([ col_fit[5],col_fit[6],col_fit[7],col_fit[8] ])
        
    #print ang_types
    #sys.exit(" debug  331")

    # Determine types
    FF_DIHTYPES = []
    DIH_PARAM = []
    ANG_IND = []
    for ind in range( len( ang_types ) ):
        
        AT_i =  ang_types[ind][0]
        AT_j =  ang_types[ind][1]
        AT_k =  ang_types[ind][2]
        AT_l =  ang_types[ind][3]
        new_type = 1
        if( debug):
            print " checking ",ind
            print "     types ",AT_i,AT_j,AT_k,AT_l
            
        for ff_i in range (len(FF_DIHTYPES)):
            FF_i = FF_DIHTYPES[ff_i][0] #.split()
            FF_j = FF_DIHTYPES[ff_i][1] #.split()
            FF_k = FF_DIHTYPES[ff_i][2] #.split()
            FF_l = FF_DIHTYPES[ff_i][3] #.split()
            if( debug): print "      pervious types  ",FF_i,FF_j,FF_k,FF_l

            if ( AT_i == FF_i  and  AT_j == FF_j  and  AT_k == FF_k  and  AT_l ==  FF_l   ):
                new_type = 0
                break
            if ( AT_i == FF_l  and  AT_j == FF_j  and  AT_k == FF_k  and  AT_l ==  FF_i   ):
                new_type = 0
                break 

        if ( new_type ):
            FF_DIHTYPES.append( [ AT_i,AT_j,AT_k,AT_l ])
            ff_i = len(FF_DIHTYPES) - 1
            NUMB_DIHTYPES = ff_i
            #DIH_PARAM.append( param_o )
            #for val in param_o:
            #    DIH_PARAM.append( val )
            ANG_IND.append( ff_i )
            if( debug): print " new type found ",ff_i
        else:
            ANG_IND.append( ff_i )
            if( debug): print " type exists  ",ff_i

    debug = 0
    if( debug ):
        print " dihedrals found ",len(ANG_IND)
        for i in range( len(ANG_IND) ):
            print FF_DIHTYPES[ANG_IND[i]],ANG_IND[i]
        
        sys.exit(" debug ")
            

    return ang_types,FF_DIHTYPES,  ANG_IND
    
def set_weights(qm_angle,qm_en_s,min_indx,max_indx):
    import sys
    import numpy as np
    
    debug = 0 
        
    wt_o = 10.0 
    wt_max = 1000.0 
    wt_min = 1000.0 
    wt_min_pm = 30.0 
    wt_min_pm2 = 20.0 

    wt_coef = 0.0
    wt_angle = []

    n_anlges =  len(qm_angle)
    
    wt_angle = [] #np.ones(n_anlges+1)
    for cent_indx in range(n_anlges+1):
        wt_angle.append(wt_o)
    
    #h = qm_angle[1] - qm_angle[0]
    #success,  tor_en,min_indx,max_indx,trans_list, k_list = prop.calc_d2(  qm_en_s , h)


    for indx in range( len(min_indx) ):
        
        cent_indx =  min_indx[indx]
        wt_angle[cent_indx] = wt_min
        if( cent_indx < len(qm_angle)  ):  wt_angle[cent_indx + 1 ] = wt_min_pm
        if( cent_indx < len(qm_angle) - 1 ):  wt_angle[cent_indx + 2 ] = wt_min_pm2
        if( cent_indx > 0 ):  wt_angle[cent_indx - 1 ] = wt_min_pm
        if( cent_indx > 1 ):  wt_angle[cent_indx - 2 ] = wt_min_pm2
        
        if(debug): print "  Setting min at ",qm_angle[cent_indx]," to ",wt_min #," with dE = ",k_list[indx]
    
    for indx in range( len(max_indx) ):
        
        cent_indx =  max_indx[indx]
        wt_angle[cent_indx] = wt_max
        
    #
    #en_minus = qm_en_s[-1] 
    #debug = 0
    #for cent_indx in range(n_anlges):
    #    qm_tor_en =  qm_en_s[ cent_indx ]
    #    # Change wieghts at inversion points
    #    if( cent_indx > 0 ): en_minus = qm_en_s[ cent_indx - 1 ] 
    #    if( cent_indx < n_anlges - 1):
    #        en_plus = qm_en_s[ cent_indx + 1 ]
    #    else:
    #        en_plus = qm_en_s[0]
    #
    #    d_minus = qm_tor_en -  en_minus
    #    d_plus  = en_plus - qm_tor_en
    #    dm_sign = 0
    #    if( d_minus < 0.0 ): dm_sign = -1 
    #    if( d_minus > 0.0 ): dm_sign =  1
    #    dp_sign = 0
    #    if( d_plus < 0.0 ): dp_sign = -1 
    #    if( d_plus > 0.0 ): dp_sign =  1
    #        
    #    if( dm_sign >  dp_sign ): wt_angle[cent_indx] = wt_max
    #    if( dm_sign <  dp_sign ):
    #        wt_angle[cent_indx] = wt_min
    #        if( cent_indx < len(qm_angle)  ):  wt_angle[cent_indx + 1 ] = wt_min_pm
    #        if( cent_indx < len(qm_angle) - 1 ):  wt_angle[cent_indx + 2 ] = wt_min_pm2
    #        if( cent_indx > 0 ):  wt_angle[cent_indx - 1 ] = wt_min_pm
    #        if( cent_indx > 1 ):  wt_angle[cent_indx - 2 ] = wt_min_pm2
    #
    #
    ## Set large values for ends "barriers"
    #wt_angle[0] = 1000.0
    #wt_angle[len(wt_angle)-1] = 1000.0

    debug = 0
    if( debug ): 
        for cent_indx in range(len(qm_angle)):
            print  qm_angle[ cent_indx ], qm_en_s[ cent_indx ], wt_angle[cent_indx]
        sys.exit(" debug weights ")
        
    wt_coef = 1.0
    
    return ( wt_angle, wt_coef)

    
    
def main():
    import os, sys
    import jsonapy, prop
    import collections
    import math, numpy 
    from scipy import optimize
    
    debug = 0 

    EVTOKCAL = 23.0605
    EVTOkJ = 96.4853
    
    options, args = get_options()

    # Store working dir  
    work_dir = os.getcwd()

    json_files = options.json.split(',')
    print json_files
    if( len(json_files) > 0 ):
	# Read index files from args
	for json_file in json_files:
	    
	    # Verbose output
	    if( options.verbose ):
		print "The molecules specified in json file ",options.json," will be read in "
    
	    json_data,json_success = jsonapy.read_jsondata(json_file)
	    if(  json_success ):
		
		mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
		ff_dih_id_list ,ff_cent_min_list ,ff_cent_max_list ,ff_cent_step_list,ff_a_k_list, ff_a_i_list, ff_a_j_list, ff_a_l_list,ff_type_list,fftor_found = jsonapy.read_ff_tor(json_data)
		qm_dih_id_list ,qm_cent_min_list ,qm_cent_max_list ,qm_cent_step_list,qm_a_k_list, qm_a_i_list, qm_a_j_list, qm_a_l_list,qmtor_found = jsonapy.read_qm_tor(json_data)
                
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found and fftor_found and  qmtor_found ):
                   
		    if( options.verbose ):
			print " Meta data found will use specified method and basis unless others are specified in the options "
		    #
		    # Construct file names 
		    #
		    #short_name = "acc%d_%s_n%d" % (accuracy, tag, number )
		    job_name = "acc%d_%s_n%d" % (accuracy, tag, n_units )
		    struct_dir = "%s/%s/" % (mol_dir, tag )
		    
		    for dih_indx in range( len(ff_dih_id_list) ):
                        
                        dih_id = ff_dih_id_list[dih_indx] 
                        cent_min = ff_cent_min_list[dih_indx]
                        cent_max = ff_cent_max_list[dih_indx]
                        cent_step = ff_cent_step_list[dih_indx]
                        a_k = ff_a_k_list[dih_indx]
                        a_i = ff_a_i_list[dih_indx]
                        a_j = ff_a_j_list[dih_indx]
                        a_l = ff_a_l_list[dih_indx]
                        ff_type_id = ff_type_list[dih_indx]
                        
                        qm_sufix = "_qm2"  # Need to read from qm_sufix_list[dih_indx]
                        ff_software = 'lammps' # Need to read from ff_software_list[dih_indx]
                    
                        dih_qm = struct_dir +'/' +job_name + '-' + dih_id + qm_sufix + '.dat'
                        dih_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + ff_software + ff_type_id +".dat"

                        options.qm_tor = dih_qm
                        options.ff_tor = dih_ff
                        
                        qm_energy = read_energy(dih_qm)
                        ff_energy = read_energy(dih_ff)
                            
                        check_dih = 1                        
                        if( ff_dih_id_list[dih_indx] != qm_dih_id_list[dih_indx] ): check_dih = 0
                        if( ff_cent_min_list[dih_indx] != qm_cent_min_list[dih_indx] ): check_dih = 0
                        if( ff_cent_max_list[dih_indx] != qm_cent_max_list[dih_indx] ): check_dih = 0
                        if( ff_cent_step_list[dih_indx] != qm_cent_step_list[dih_indx] ): check_dih = 0
                        if( ff_a_k_list[dih_indx] != qm_a_k_list[dih_indx] ): check_dih = 0
                        if( ff_a_i_list[dih_indx] != qm_a_i_list[dih_indx] ): check_dih = 0
                        if( ff_a_j_list[dih_indx] != qm_a_j_list[dih_indx] ): check_dih = 0
                        if( ff_a_l_list[dih_indx] != qm_a_l_list[dih_indx] ): check_dih = 0
                        if( len(qm_energy) != len(ff_energy) ): check_dih = 0
                        
                        if( check_dih ):
                            
                            #  max and min
                            qm_en_max = -1000000000.0
                            qm_en_min =  1000000000.0
                            ff_en_max = -1000000000.0
                            ff_en_min =  1000000000.0
                            t_en_max  = -1000000000.0
                            t_en_min  =  1000000000.0    
                        
                            # Loop over energies and find maximum and minimum 
                            for cent_indx in range(len(qm_energy)):
                                qm_tor_en = float(qm_energy[cent_indx][2])
                                
                                if(debug): print qm_energy[cent_indx][1],qm_tor_en
                                
                                if( qm_tor_en > qm_en_max): qm_en_max = qm_tor_en
                                if( qm_tor_en < qm_en_min): qm_en_min = qm_tor_en
                        
                                ff_en = float(ff_energy[cent_indx][3])
                        
                                if( ff_en > ff_en_max): ff_en_max = ff_en
                                if( ff_en < ff_en_min): ff_en_min = ff_en
                                
                                if(debug): print ff_energy[cent_indx][1],ff_en
                        
                            debug = 1
                            
                            if(debug):print " qm max ; min ",qm_en_max,qm_en_min
                            if(debug):print " ff max ; min ",ff_en_max,ff_en_min
                        
                            # Read in energies and shift by the minimum 
                            target_en = []
                            qm_en_s   = []
                            qm_angle  = []
                            ff_en_s   = []
                            ff_angles = []
                            for cent_indx in range(len(qm_energy)):
                                qm_x = float(qm_energy[cent_indx][1])        
                                qm_en = float(qm_energy[cent_indx][2])
                                qm_s = qm_en - qm_en_min
                                ff_en = float(ff_energy[cent_indx][3])
                                ff_s = ff_en - ff_en_min
                                d_en = ( qm_s )  - ( ff_s )
                                # FF angles 
                                ang_val = []
                                ang_val.append( float(ff_energy[cent_indx][7]) )
                                ang_val.append( float(ff_energy[cent_indx][8]) ) 
                                ang_val.append( float(ff_energy[cent_indx][9]) ) 
                                ang_val.append( float(ff_energy[cent_indx][10]) )
                        
                                target_en.append(d_en)
                                qm_angle.append(qm_x)
                                qm_en_s.append(qm_s)
                                ff_en_s.append(ff_s)
                                ff_angles.append( ang_val )
                                
                                if( d_en > t_en_max): t_en_max = d_en
                                if( d_en < t_en_min): t_en_min = d_en
                                    
                            if(debug): print " target  max ; min ",t_en_max,t_en_min
                        
                            # Create list of target values 
                            target_en_s = []
                            for t_indx in range( len(target_en)):
                                # this causes functional issues 
                                #target_en_s.append( target_en[t_indx] - t_en_min )
                                target_en_s.append( target_en[t_indx] - target_en[0] )
                                
                            # Find second derivates for max/min

                            h = qm_angle[1] - qm_angle[0]
                            success,  tor_en,min_indx,max_indx,trans_list,trans_indxs, k_list = prop.calc_d2(  qm_en_s , h)
                            
                        
                            # Minimum 
                            for indx in range( len(min_indx) ):
                                cent_indx =  min_indx[indx]
                                if( cent_indx > 0 ):
                                    print " Minimum qm_en_s meV ",qm_en_s[cent_indx-1]*1000 , qm_en_s[cent_indx]*1000 ,  qm_en_s[cent_indx+1]*1000
                                    print " Minimum target_en_s meV ",target_en_s[cent_indx-1]*1000 , target_en_s[cent_indx]*1000 ,  target_en_s[cent_indx+1]*1000

                            # Set wieghts
                            wt_angle, wt_coef = set_weights(qm_angle,qm_en_s,min_indx,max_indx)
                            #
                            #    
                            #for indx in range( len(trans_list) ):
                            #    
                            #    max_i =  trans_indxs[indx][0]
                            #    min_i =  trans_indxs[indx][1]
                            #    
                            #    trans_val = qm_en_s[max_i] - qm_en_s[min_i]
                            #    
                            #    print "  trans ",indx," calcd ",trans_val," stored ",trans_list[indx]
                            #    
                            #    
                            #sys.exit(" rpint transys to check ! ")
                            #
                            if( options.verbose ):
                                print " Read in ",len(qm_en_s)," energies ",min(qm_en_s),max(qm_en_s)
                                
                            
                            # Read in dihedral information 
                            ang_types , FF_DIHTYPES, ANG_IND =  read_dihtypes(struct_dir,job_name,options)
                            n_dihtypes = max( ANG_IND ) + 1
                            
                            if( options.verbose ):
                                print "  with ",n_dihtypes," angle types "
                                for a_indx in range( len(ang_types) ):
                                    ang_i = ANG_IND[a_indx]
                                    print ang_types[a_indx]," type ",ang_i,FF_DIHTYPES[ang_i]
                                                               
                            fit_four_4 = 1
                            if( fit_four_4 ):
                                
                                os.chdir(struct_dir)
                                
                                # Initial single fourier series fit 
                                tor_type = "four"
                                n_param = 4
                                
                                param_o = [0.0,0.0,0.0,0.0]
                                
                                use_biaryl = 1
                                if( use_biaryl ):
                                    print " Using biaryl parameters "
                                
                                    DIH_PARAM = []
                                    fit_param = []
                                      
                                    for param_ind in range( n_dihtypes ):
                                        #F_coef = []
                                        #for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        #    biaryl_param.append(  0.0 )
                                        DIH_PARAM.append( param_o )
                                        
                                        fit_param.append(  1 )      # default to be fitted
                                        
                                        print "index ",param_ind," has  types ",ang_types[param_ind]
                                        
                                        print " check biarly parameters for ",ang_types[param_ind][0] , ang_types[param_ind][1] ,ang_types[param_ind][2] ,ang_types[param_ind][1]
                                        
                                        if( ang_types[param_ind][0] == 'CS' and ang_types[param_ind][1] == 'C!' and  ang_types[param_ind][2] == 'C!' and ang_types[param_ind][3] == 'CS'  ):
                                            
                                            DIH_PARAM[param_ind] = [0.0,1.49/EVTOKCAL/4.0,0.0,0.0]  # CA - C! - CP - CS 
                                            
                                            fit_param[param_ind] = 0
                                            
                                        elif( ang_types[param_ind][0] == 'S' and ang_types[param_ind][1] == 'C!' and  ang_types[param_ind][2] == 'C!' and ang_types[param_ind][3] == 'CS'  ):
                                            
                                            DIH_PARAM[param_ind] = [0.0,1.33/EVTOKCAL/4.0,0.0,0.0]  # CA - C! - CP - S 
                                            
                                            fit_param[param_ind] = 0
                                            
                                        
                                        elif( ang_types[param_ind][0] == 'S' and ang_types[param_ind][1] == 'C!' and  ang_types[param_ind][2] == 'C!' and ang_types[param_ind][3] == 'S'  ):
                                            
                                            # Fit S-C!-C!-S 
                                            DIH_PARAM[param_ind] = [0.0,1.33/EVTOKCAL/4.0,0.0,0.0]  # CA - C! - CP - S 
                                            fit_param[param_ind] = 1
                                            
                                        elif( ang_types[param_ind][0] == 'CA' and ang_types[param_ind][1] == 'C!' and  ang_types[param_ind][2] == 'C!' and ang_types[param_ind][3] == 'CA'  ):
                                            
                                            # Fit CA-C!-C!-CA 
                                            DIH_PARAM[param_ind] = [0.0,1.97/EVTOKCAL/4.0,0.0,0.0]  # CA - C! - CP - S 
                                            fit_param[param_ind] = 1
                                            
                                    fit_param_list = []
                                    
                                    for param_ind in range( n_dihtypes ):
                                        if( fit_param[param_ind] == 1 ):
                                            print " Adding ",param_ind," to fit list x  "
                                            for param_i in DIH_PARAM[param_ind]:
                                                fit_param_list.append(param_i)
        
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                p_cnt = -1
                                for param_ind in range( n_dihtypes ):
                                    F_coef = DIH_PARAM[param_ind]
                                    if(fit_param[param_ind] == 1 ):
                                        p_cnt += 1
                                        F_coef = []
                                        print "For index ",param_ind,"  will use fitted parameters "
                                        for p_indx in range( p_cnt*n_param,p_cnt*n_param + n_param):
                                            print " append fit_param_list ",p_indx
                                            F_coef.append(  fit_param_list[p_indx] )
                                    
                                    print "  Index ",param_ind," has  types ",ang_types[param_ind]
                                    print "  Initial RB coefficents from average individual fits ",F_coef
                                    # print "                                                      ",F_set
                                                                            
                                
                                
                                rmsd_en,rmsd_trans = calc_rms_v2( fit_param_list, fit_param, DIH_PARAM, ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len param_list ",len(fit_param_list)
                                #
                                resid = residuals_v6( fit_param_list, fit_param, DIH_PARAM, ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                
                                #sys.exit( "DIH_PARAM testing 1 ") 
                                
                                #                    
                                #for param_ind in range( n_dihtypes ):
                                #    F_coef = []
                                #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                #        F_coef.append(  p_rb[p_indx] )
                                #    
                                #    print "  F coefficents ",F_coef
                                #    
                                param_list_pre = []
                                for p_indx in range( len(fit_param_list)):
                                    param_list_pre.append( fit_param_list[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    fit_param_list, successp = optimize.leastsq( residuals_v6, fit_param_list, args=( fit_param, DIH_PARAM, ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(fit_param_list)):
                                        d_pf = fit_param_list[p_indx] - param_list_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    print "   fit_iter ",fit_iter," with delat_param ",d_param
                                    if( d_param <  1.0002 ):
                                        delta_param = False
                                    
                                    p_fourieparam_list_prer_pre = []
                                    for p_indx in range( len(fit_param_list)):
                                        param_list_pre.append( fit_param_list[p_indx] )
                                        
                                resid = residuals_v6( fit_param_list, fit_param, DIH_PARAM, ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                #
                                # Update DIH_PARAM to fitted
                                #
                                p_cnt = -1
                                for param_ind in range( n_dihtypes ):
                                    F_coef = DIH_PARAM[param_ind]
                                    if(fit_param[param_ind] == 1 ):
                                        p_cnt += 1
                                        F_coef = []
                                        print "For index ",param_ind,"  will use fitted parameters "
                                        for p_indx in range( p_cnt*n_param,p_cnt*n_param + n_param):
                                            print " append fit_param_list ",p_indx
                                            F_coef.append(  fit_param_list[p_indx] )
                                        DIH_PARAM[param_ind] = F_coef
                                        
                                    print "  Index ",param_ind," has  types ",ang_types[param_ind]
                                    print "  Initial RB coefficents from average individual fits ",F_coef
                                    # print "                                                      ",F_set
                                                                            
                                                                

                                rmsd_en,rmsd_trans = calc_rms_v2( fit_param_list, fit_param, DIH_PARAM, ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                                                                                    
                                for param_ind in range( n_dihtypes ):
                                    print "  F coefficents ",DIH_PARAM[param_ind]
                                                           
                                # sys.exit(" F2 fit test ")

                                #
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor_v2(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                                    
                                
                                # Print parameters for itp file
                                #  the factor of 2.0 will be removed once the parameters are normalized
                                
                                print "   Gromacs format x2 for compatability with biaryl dih parameters "
                                for a_indx in range( len(ang_types) ):
                                    param_ind = ANG_IND[a_indx]
                                    dih_id = ""
                                    for dih_indx in FF_DIHTYPES[param_ind]:
                                        dih_id += str( dih_indx ) + "  " 
                                
                                    if( tor_type == "RB" ):
                                            
                                        F_coef = ''
                                        for p_indx in range( len( DIH_PARAM[param_ind]) ):
                                            F_coef +=  str(round( DIH_PARAM[param_ind][p_indx]*EVTOkJ*2.0 ,4) ) + "  "
                                            
                                        print dih_id,"  3  ",F_coef, "   0.000000  # fit 1 "
                                    
                                    elif(   tor_type == "four" ):
                                            
                                        F_coef = []
                                        for p_indx in range( len( DIH_PARAM[param_ind]) ):
                                            # F_coef +=  str(round( DIH_PARAM[param_ind][p_indx]*EVTOkJ*2.0 ,4) ) + "  "
                                            F_coef.append( round( DIH_PARAM[param_ind][p_indx]*EVTOkJ*2.0 ,4) )
                                            
                                        C0 = F_coef[1]  + 0.5*( F_coef[0] + F_coef[2] )
                                        C1 =0.5*( -1.0*F_coef[0] + 3.0*F_coef[2] )
                                        C2 = -1.0*F_coef[1] + 4.0*F_coef[3]
                                        C3 = -2.0*F_coef[2]
                                        C4 = -4.0*F_coef[3]
                                        C5 = 0.0
                                        print dih_id,"  3  ",C0,C1,C2,C3,C4,C5, "  # fit 1 "
                            
                                print "  Lammps format  "
                                for a_indx in range( len(ang_types) ):
                                    param_ind = ANG_IND[a_indx]
                                    dih_id = ""
                                    for dih_indx in FF_DIHTYPES[param_ind]:
                                        dih_id += str( dih_indx ) + "  "
                                    
                                    if( tor_type == "RB" ):
                                            
                                        RB_coef = []
                                        
                                        for p_indx in range( len( DIH_PARAM[param_ind]) ):
                                            # F_coef +=  str(round( DIH_PARAM[param_ind][p_indx]*EVTOkJ*2.0 ,4) ) + "  "
                                            RB_coef.append( round( DIH_PARAM[param_ind][p_indx]*EVTOKCAL ,4) )
                                            
                                        F1 = -1.0*( 2.0*RB_coef[1] + 3.0*RB_coef[3]/2.0)
                                        F2 = -1.0*( RB_coef[2] + RB_coef[4])
                                        F3 = -0.5*RB_coef[3]
                                        F4 = -0.25*RB_coef[4]
                        
                                        print "  1  ",F1,F2,F3,F4, "  #  ",dih_id
                                    
                                    elif(   tor_type == "four" ):
                                        
                                        
                                        F_coef = ''
                                        for p_indx in range( len( DIH_PARAM[param_ind]) ):
                                            F_coef +=  str(round( DIH_PARAM[param_ind][p_indx]*EVTOKCAL ,4) ) + "  "
                                            
                                        print "  1  ",F_coef, "  #  ",dih_id
                                        
                                sys.exit(" Fourier fitting finished ")
                                    
                            #
                            ## Set dihdral types
                            ##tor_type = "four"
                            ##n_param = 4 
                            ##tor_type = "RB"
                            ##n_param = 5
                            #
                            #fit_cos = 0
                            #if( fit_cos ):
                            #    # Fit cos "harmonic" function 
                            #    tor_type = "harmonic" 
                            #    n_param = 3
                            #    
                            #    param_o = []
                            #    for n in range(n_param):
                            #        param_o.append(0.0)
                            #    #
                            #    # Set initial parmeters
                            #    #
                            #    param_o[0] = options.tor_paramo #*10000
                            #    param_o[1] =  2.0
                            #    param_o[2] =  1.0
                            #    #
                            #    # set parameter list consisting of each different dihedral's parameters
                            #    #
                            #    p_kcos = []
                            #    #
                            #    for iparam_indnd in range(n_dihtypes):
                            #        for val in param_o:
                            #            p_kcos.append( val )
                            #            
                            #    
                            #    print "  n_dihtypes ",n_dihtypes
                            #    print "  param_o",len(param_o)
                            #    print "  len p_kcos ",len(p_kcos)
                            #    print "  len p_kcos ",len(p_kcos)
                            #    #
                            #    resid = residuals_kcos(p_kcos,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                            #    
                            #    error = 0.0
                            #    for er in resid:
                            #        error += er*er
                            #    print " Initial error ",error
                            #    #sys.exit( "DIH_PARAM") 
                            #    
                            #                        
                            #    for param_ind in range( n_dihtypes ):
                            #        F_coef = []
                            #        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                            #            F_coef.append(  p_kcos[p_indx] )
                            #        
                            #        print "  Harmonic coefficents "
                            #        print "                        k_o ",float(F_coef[0])*EVTOKCAL," kcal/mol "
                            #        print "                        mult ",F_coef[1]
                            #        print "                        phase ",F_coef[2],float(F_coef[2])*math.pi
                            #
                            #
                            #    p_kcos,successp = optimize.leastsq(residuals_kcos,p_kcos,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param),epsfcn=0.0001)
                            #
                            #
                            #    resid = residuals_kcos(p_kcos,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                            #    error = 0.0 
                            #    for er in resid:
                            #        error += er*er
                            #    print " final error ",error
                            #    
                            #    
                            #                        
                            #    for param_ind in range( n_dihtypes ):
                            #        F_coef = []
                            #        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                            #            F_coef.append(  p_kcos[p_indx] )
                            #        
                            #        print "  Harmonic coefficents "
                            #        print "                        k_o ",float(F_coef[0])/10000*EVTOKCAL," kcal/mol "
                            #        print "                        mult ",F_coef[1]
                            #        print "                        phase ",F_coef[2],float(F_coef[2])*math.pi
                            #
                            #    os.chdir(struct_dir)
                            #        
                            fit_rb = 0
                            if( fit_rb ):
                                
                                os.chdir(struct_dir)
                                
                                tor_type = "RB"
                                n_param = 5
                                
                                param_o = []
                                for n in range(n_param):
                                    param_o.append(0.0)
                                #
                                # Set initial parmeters
                                #
                                param_o[0] =  options.tor_paramo 
                                param_o[2] =  -1.0*options.tor_paramo 
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                p_rb = []
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    #F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        p_rb.append(  0.0 )
                                    #p_fourier[param_ind*n_param +  1 ] =  p_fourierF2[param_ind] 



                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  RB coefficents ",F_coef
                                    
                                                                
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len p_fourier ",len(p_rb)
                                #
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  F coefficents ",F_coef
                                    
                                p_fourier_pre = []
                                for p_indx in range( len(p_rb)):
                                    p_fourier_pre.append( p_rb[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    p_rb,successp = optimize.leastsq(residuals_v5,p_rb,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(p_rb)):
                                        d_pf = p_rb[p_indx] - p_fourier_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    if( d_param <  0.0002 ):
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        delta_param = False
                                    
                                    p_fourier_pre = []
                                    for p_indx in range( len(p_rb)):
                                        p_fourier_pre.append( p_rb[p_indx] )
                                        
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  round(p_rb[p_indx]*EVTOKCAL,4) )
                                    
                                    print "  F coefficents ",F_coef
                                    #print "  F coefficents ",F_coef[0]*EVTOKCAL,F_coef[1]*EVTOKCAL,F_coef[2]*EVTOKCAL,F_coef[3]*EVTOKCAL
                                    
                                rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                DIH_PARAM = p_rb                                    
                                    
                                #sys.exit(" F2 fit test ")

                                #
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                  
                            fit_rb_F2 = 0
                            if( fit_rb_F2 ):
                                
                                os.chdir(struct_dir)
                                
                                tor_type = "four_F2"
                                n_param = 1
                                
                                param_o = []
                                for n in range(n_param):
                                    param_o.append(0.0)
                                    
                                #
                                # Set initial parmeters
                                #
                                param_o[0] =  options.tor_paramo 
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                p_fourierF2 = []
                                #
                                for iparam_indnd in range(n_dihtypes):
                                    for val in param_o:
                                        p_fourierF2.append( val )
                                        
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  param_o",len(param_o)
                                print "  len p_fourierF2 ",len(p_fourierF2)
                                #
                                resid = residuals_fourF2(p_fourierF2,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourierF2[p_indx] )
                                    
                                    print "  F2 coefficents "
                                    print "                        F2 ",float(F_coef[0])," eV ",float(F_coef[0])*EVTOKCAL," kcal/mol "
                            
                                p_fourierF2,successp = optimize.leastsq(residuals_fourF2,p_fourierF2,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param),epsfcn=0.0001)
                            
    
                                resid = residuals_fourF2(p_fourierF2,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourierF2[p_indx] )
                                    
                                    print "  F2 coefficents "
                                    print "                        F2 ",float(F_coef[0])," eV ",float(F_coef[0])*EVTOKCAL," kcal/mol "
                                    
                                DIH_PARAM = p_fourierF2                                    
                                
                                # 
                                tor_type = "RB"
                                n_param = 5
                                
                                param_o = []
                                for n in range(n_param):
                                    param_o.append(0.0)
                                #
                                # Set initial parmeters
                                #
                                param_o[0] =  options.tor_paramo 
                                param_o[2] =  -1.0*options.tor_paramo 
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                p_rb = []
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    #F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        p_rb.append(  0.0 )
                                    #p_fourier[param_ind*n_param +  1 ] =  p_fourierF2[param_ind] 
                                    p_rb[param_ind*n_param  ] =  1.0*p_fourierF2[param_ind] 
                                    p_rb[param_ind*n_param + 2 ] =  -1.0*p_fourierF2[param_ind] 



                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  RB coefficents ",F_coef
                                    
                                                                
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len p_fourier ",len(p_rb)
                                #
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  F coefficents ",F_coef
                                    
                                p_fourier_pre = []
                                for p_indx in range( len(p_rb)):
                                    p_fourier_pre.append( p_rb[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    p_rb,successp = optimize.leastsq(residuals_v5,p_rb,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(p_rb)):
                                        d_pf = p_rb[p_indx] - p_fourier_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    if( d_param <  0.0002 ):
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        delta_param = False
                                    
                                    p_fourier_pre = []
                                    for p_indx in range( len(p_rb)):
                                        p_fourier_pre.append( p_rb[p_indx] )
                                        
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  round(p_rb[p_indx]*EVTOKCAL,4) )
                                    
                                    print "  F coefficents ",F_coef
                                    #print "  F coefficents ",F_coef[0]*EVTOKCAL,F_coef[1]*EVTOKCAL,F_coef[2]*EVTOKCAL,F_coef[3]*EVTOKCAL
                                    
                                rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                DIH_PARAM = p_rb                                    
                                    
                                #sys.exit(" F2 fit test ")

                                #
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                            fit_rb_2 = 0
                            if( fit_rb_2 ):
                                
                                os.chdir(struct_dir)
                                
                                # Initial single fourier series fit 
                                tor_type = "RB"
                                n_param = 5
                                
                                param_o = []
                                for n in range(n_param):
                                    param_o.append(0.0)
                                
                                guess_ind = 0
                                if( guess_ind ):
                                    # Guess initial parameters by fitting each component individually
                                        
                                    #
                                    # Set initial parmeters
                                    #
                                    # param_o[0] =  options.tor_paramo 
                                    #
                                    # set parameter list consisting of each different dihedral's parameters
                                    #
                                    p_initial = []
                                    #
                                    print "  param_o",len(param_o)
                                    #
                                    resid = residuals_initial(param_o,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                    error = 0.0
                                    for er in resid:
                                        error += er*er
                                    print " Initial error ",error
                                    #sys.exit( "DIH_PARAM") 
                                    
                                    F_coef = []
                                    print "  Initial fourier coefficents "
                                    for p_indx in range( n_param ):
                                        F_coef.append(  round(param_o[p_indx]*EVTOKCAL,4) )
                                        
                                    print "  F coefficents ",F_coef
                                
                                    #p_fourierF2,successp = optimize.leastsq(residuals_fourF2,param_o,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param),epsfcn=0.0001)
                                    param_o,successp = optimize.leastsq(residuals_initial,param_o,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                
        
                                    resid = residuals_initial(param_o,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    error = 0.0 
                                    for er in resid:
                                        error += er*er
                                    print " final error ",error
                                    
                                    
                                                                                        
                                    F_coef = []
                                    for p_indx in range( n_param ):
                                        F_coef.append(  round(param_o[p_indx],4) )
                                        
                                    print "  Initial coefficents ",F_coef
                                    
                                    # 
                                    #
                                    # set parameter list consisting of each different dihedral's parameters
                                    #
                                    p_rb = []
                                    
                                                        
                                    for param_ind in range( n_dihtypes ):
                                        #F_coef = []
                                        initial_i = -1
                                        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                            initial_i += 1 
                                            p_rb.append(  param_o[initial_i]/float(len(ang_types)) )
                                        
                                    rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                        
                                    DIH_PARAM = p_rb
                                    
                                else:
                                    # set initial parameters to zero
                                    # set parameter list consisting of each different dihedral's parameters
                                    #
                                    p_rb = []
                                   
                                                       
                                    for param_ind in range( n_dihtypes ):
                                        #F_coef = []
                                        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                            p_rb.append(  0.0 )
                                    
                                #sys.exit(" F2 fit test ")

                                # Collapse ff_angles into a single set of angle 1
                                # 
                                single_angles = []
                                for ang_i in range(len(ff_angles)):
                                    single_angles.append( [ ff_angles[ang_i][0] ] )
                                
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,param_o,ANG_IND, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                #if( options.plot_tor ):
                                #    plot_tor(options,param_o,ANG_IND, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                #    plot_tor_comp(options,param_o,ANG_IND,FF_DIHTYPES, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                # sys.exit(" Fourier fitting initialization 3 ")

                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  RB coefficents ",F_coef
                                    
                                                                
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len p_fourier ",len(p_rb)
                                #
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  F coefficents ",F_coef
                                    
                                p_fourier_pre = []
                                for p_indx in range( len(p_rb)):
                                    p_fourier_pre.append( p_rb[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    p_rb,successp = optimize.leastsq(residuals_v5,p_rb,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(p_rb)):
                                        d_pf = p_rb[p_indx] - p_fourier_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    if( d_param <  0.0002 ):
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        delta_param = False
                                    
                                    p_fourier_pre = []
                                    for p_indx in range( len(p_rb)):
                                        p_fourier_pre.append( p_rb[p_indx] )
                                        
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  round(p_rb[p_indx]*EVTOKCAL,4) )
                                    
                                    print "  F coefficents ",F_coef
                                    #print "  F coefficents ",F_coef[0]*EVTOKCAL,F_coef[1]*EVTOKCAL,F_coef[2]*EVTOKCAL,F_coef[3]*EVTOKCAL
                                    
                                rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                DIH_PARAM = p_rb                                    
                                    
                                #sys.exit(" F2 fit test ")

                                #
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                            fit_rb_3 = 0
                            if( fit_rb_3 ):
                                
                                os.chdir(struct_dir)
                                
                                # Initial single fourier series fit 
                                tor_type = "RB"
                                n_param = 4
                                
                                guess_ind = 0
                                if( guess_ind ):
                                        
                                    #
                                    p_initial = []
                                    #
                                    # Fit a set of fourier coefficients for each angle
                                    for ang_i in range( len(ff_angles[0]) ):
                                        angles_list = ff_angles[:][ang_i]
                                        # print " angle set ",ang_i
                                        dtype_angles = []
                                        for angle_indx in range(len(ff_angles)):
                                            #print ff_angles[angle_indx][ang_i]
                                            dtype_angles.append( ff_angles[angle_indx][ang_i]  )
                                    
                                            
                                        #
                                        # Set initial parmeters
                                        #
                                        param_o = []
                                        for n in range(n_param):
                                            param_o.append(0.0)
                                        #
                                        print "  param_o",len(param_o)
                                        #
                                            
                                        F_coef = []
                                        print "  Initial fourier coefficents "
                                        for p_indx in range( n_param ):
                                            F_coef.append(  round(param_o[p_indx]*EVTOKCAL,4) )
                                            
                                        print "  F coefficents ",F_coef
                                    
                                    
                                        resid = residuals_initialv2(param_o,ANG_IND,target_en_s,dtype_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                        
                                        error = 0.0
                                        for er in resid:
                                            error += er*er
                                        print " Initial error ",error
                                        #sys.exit( "DIH_PARAM") 
    
                                    
                                        #p_fourierF2,successp = optimize.leastsq(residuals_fourF2,param_o,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param),epsfcn=0.0001)
                                        param_o,successp = optimize.leastsq(residuals_initialv2,param_o,args=(ANG_IND,target_en_s,dtype_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                
        
                                        resid = residuals_initialv2(param_o,ANG_IND,target_en_s,dtype_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                        error = 0.0 
                                        for er in resid:
                                            error += er*er
                                        print " final error ",error
                                    
                                    
                                                                                        
                                        F_coef = []
                                        for p_indx in range( n_param ):
                                            F_coef.append(  round(param_o[p_indx],4) )
                                            
                                        print "  Initial fitted coefficents ",F_coef
                                        
                                        
                                        p_initial.append( param_o )
                                        
                                    for a_indx in range( len(ang_types) ):
                                        ang_i = ANG_IND[a_indx]
                                        print ang_types[a_indx]," type ",ang_i," parameters ",p_initial[a_indx]
                                        
                                        
                                    # Find dihedrals for each type and average initial parameters
                                    p_rb = []
                                    for a_type in range( n_dihtypes):
                                        param_sum = []
                                        for a_indx in range( len(ang_types) ):
                                            ang_i = ANG_IND[a_indx]
                                            if( a_type ==  ang_i ):
                                                param_sum.append( p_initial[a_indx] )
                                                
                                        print a_type , param_sum,len(param_sum )
                                        for p_i in range(n_param):
                                            sum_p_i = 0 
                                            for indx in range( len(param_sum )) :
                                                sum_p_i += param_sum[indx][p_i] 
                                            
                                            
                                                print a_type,p_i,param_sum[indx][p_i]  #numpy.mean( param_sum[p_i]  )
                                                
                                            print "   average of ",p_i ,  sum_p_i/len(param_sum ) 
                                            p_rb.append(  sum_p_i/float(len(param_sum ))/float(len(ang_types) ) )

                                else:
                                    # set initial parameters to zero
                                    # set parameter list consisting of each different dihedral's parameters
                                    #
                                    p_rb = []
                                   
                                                       
                                    for param_ind in range( n_dihtypes ):
                                        #F_coef = []
                                        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                            p_rb.append(  0.0 )
                                                                                
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  Initial RB coefficents from average individual fits ",F_coef
                                    
                                                                
                                                                                                        
                                #sys.exit(" fit_rb_3 1  ")
                                
                                                    
                                rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                DIH_PARAM = p_rb                                    
                                    
                                #sys.exit(" F2 fit test ")

                                # Collapse ff_angles into a single set of angle 1
                                # 
                                #single_angles = []
                                #for ang_i in range(len(ff_angles)):
                                #    single_angles.append( [ ff_angles[ang_i][0] ] )
                                
                                # Fitted torsional energies
                                #fitted_ptor , fit_toten , fitted_comp = print_tor(options,param_o,ANG_IND, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                #if( options.plot_tor ):
                                #    plot_tor(options,param_o,ANG_IND, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                #    plot_tor_comp(options,param_o,ANG_IND,FF_DIHTYPES, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                #rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                # Fitted torsional energies
                                #fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                #if( options.plot_tor ):
                                #    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                #    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                #sys.exit(" Fourier fitting initialization 3 ")
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len p_fourier ",len(p_rb)
                                #
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                #                    
                                #for param_ind in range( n_dihtypes ):
                                #    F_coef = []
                                #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                #        F_coef.append(  p_rb[p_indx] )
                                #    
                                #    print "  F coefficents ",F_coef
                                #    
                                p_fourier_pre = []
                                for p_indx in range( len(p_rb)):
                                    p_fourier_pre.append( p_rb[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    p_rb,successp = optimize.leastsq(residuals_v5,p_rb,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(p_rb)):
                                        d_pf = p_rb[p_indx] - p_fourier_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    if( d_param <  0.0002 ):
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        delta_param = False
                                    
                                    p_fourier_pre = []
                                    for p_indx in range( len(p_rb)):
                                        p_fourier_pre.append( p_rb[p_indx] )
                                        
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  round(p_rb[p_indx]*EVTOKCAL,4) )
                                    
                                    print "  F coefficents ",F_coef
                                    #print "  F coefficents ",F_coef[0]*EVTOKCAL,F_coef[1]*EVTOKCAL,F_coef[2]*EVTOKCAL,F_coef[3]*EVTOKCAL
                                    
                                rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                DIH_PARAM = p_rb                                    
                                    
                                #sys.exit(" F2 fit test ")

                                #
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                            fit_fourier_3 = 0
                            if( fit_fourier_3 ):
                                
                                os.chdir(struct_dir)
                                
                                # Initial single fourier series fit 
                                tor_type = "four"
                                n_param = 4
                                
                                guess_ind = 0
                                if( guess_ind ):
                                        
                                    #
                                    p_initial = []
                                    #
                                    # Fit a set of fourier coefficients for each angle
                                    for ang_i in range( len(ff_angles[0]) ):
                                        angles_list = ff_angles[:][ang_i]
                                        # print " angle set ",ang_i
                                        dtype_angles = []
                                        for angle_indx in range(len(ff_angles)):
                                            #print ff_angles[angle_indx][ang_i]
                                            dtype_angles.append( ff_angles[angle_indx][ang_i]  )
                                    
                                            
                                        #
                                        # Set initial parmeters
                                        #
                                        param_o = []
                                        for n in range(n_param):
                                            param_o.append(0.0)
                                        #
                                        print "  param_o",len(param_o)
                                        #
                                            
                                        F_coef = []
                                        print "  Initial fourier coefficents "
                                        for p_indx in range( n_param ):
                                            F_coef.append(  round(param_o[p_indx]*EVTOKCAL,4) )
                                            
                                        print "  F coefficents ",F_coef, "  kcal/mol "
                                    
                                    
                                        resid = residuals_initialv2(param_o,ANG_IND,target_en_s,dtype_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                        
                                        error = 0.0
                                        for er in resid:
                                            error += er*er
                                        print " Initial error ",error
                                        #sys.exit( "DIH_PARAM") 
    
                                    
                                        #p_fourierF2,successp = optimize.leastsq(residuals_fourF2,param_o,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param),epsfcn=0.0001)
                                        param_o,successp = optimize.leastsq(residuals_initialv2,param_o,args=(ANG_IND,target_en_s,dtype_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                
        
                                        resid = residuals_initialv2(param_o,ANG_IND,target_en_s,dtype_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                        error = 0.0 
                                        for er in resid:
                                            error += er*er
                                        print " final error ",error
                                    
                                    
                                                                                        
                                        F_coef = []
                                        for p_indx in range( n_param ):
                                            F_coef.append(  round(param_o[p_indx],4) )
                                            
                                        print "  Initial fitted coefficents ",F_coef
                                        
                                        
                                        p_initial.append( param_o )
                                        
                                    for a_indx in range( len(ang_types) ):
                                        ang_i = ANG_IND[a_indx]
                                        print ang_types[a_indx]," type ",ang_i," parameters ",p_initial[a_indx]
                                        
                                        
                                    # Find dihedrals for each type and average initial parameters
                                    p_rb = []
                                    for a_type in range( n_dihtypes):
                                        param_sum = []
                                        for a_indx in range( len(ang_types) ):
                                            ang_i = ANG_IND[a_indx]
                                            if( a_type ==  ang_i ):
                                                param_sum.append( p_initial[a_indx] )
                                                
                                        print a_type , param_sum,len(param_sum )
                                        for p_i in range(n_param):
                                            sum_p_i = 0 
                                            for indx in range( len(param_sum )) :
                                                sum_p_i += param_sum[indx][p_i] 
                                            
                                            
                                                print a_type,p_i,param_sum[indx][p_i]  #numpy.mean( param_sum[p_i]  )
                                                
                                            print "   average of ",p_i ,  sum_p_i/len(param_sum ) 
                                            p_rb.append(  sum_p_i/float(len(param_sum ))/float(len(ang_types) ) )

                                else:
                                    # set initial parameters to zero
                                    # set parameter list consisting of each different dihedral's parameters
                                    #
                                    p_rb = []
                                   
                                                       
                                    for param_ind in range( n_dihtypes ):
                                        #F_coef = []
                                        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                            p_rb.append(  0.0 )
                                                                                
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_rb[p_indx] )
                                    
                                    print "  Initial RB coefficents from average individual fits ",F_coef
                                    
                                                                
                                                                                                        
                                #sys.exit(" fit_rb_3 1  ")
                                
                                                    
                                rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                DIH_PARAM = p_rb                                    
                                    
                                #sys.exit(" F2 fit test ")

                                # Collapse ff_angles into a single set of angle 1
                                # 
                                #single_angles = []
                                #for ang_i in range(len(ff_angles)):
                                #    single_angles.append( [ ff_angles[ang_i][0] ] )
                                
                                # Fitted torsional energies
                                #fitted_ptor , fit_toten , fitted_comp = print_tor(options,param_o,ANG_IND, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                #if( options.plot_tor ):
                                #    plot_tor(options,param_o,ANG_IND, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                #    plot_tor_comp(options,param_o,ANG_IND,FF_DIHTYPES, single_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                #rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                # Fitted torsional energies
                                #fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                #if( options.plot_tor ):
                                #    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                #    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                #sys.exit(" Fourier fitting initialization 3 ")
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len p_fourier ",len(p_rb)
                                #
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                #                    
                                #for param_ind in range( n_dihtypes ):
                                #    F_coef = []
                                #    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                #        F_coef.append(  p_rb[p_indx] )
                                #    
                                #    print "  F coefficents ",F_coef
                                #    
                                p_fourier_pre = []
                                for p_indx in range( len(p_rb)):
                                    p_fourier_pre.append( p_rb[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    p_rb,successp = optimize.leastsq(residuals_v5,p_rb,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(p_rb)):
                                        d_pf = p_rb[p_indx] - p_fourier_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    if( d_param <  0.0002 ):
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        delta_param = False
                                    
                                    p_fourier_pre = []
                                    for p_indx in range( len(p_rb)):
                                        p_fourier_pre.append( p_rb[p_indx] )
                                        
                                resid = residuals_v5(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  round(p_rb[p_indx]*EVTOKCAL,4) )
                                    
                                    print "  F coefficents ",F_coef
                                    #print "  F coefficents ",F_coef[0]*EVTOKCAL,F_coef[1]*EVTOKCAL,F_coef[2]*EVTOKCAL,F_coef[3]*EVTOKCAL
                                    
                                rmsd_en,rmsd_trans = calc_rms(p_rb,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                DIH_PARAM = p_rb                                    
                                    
                                #sys.exit(" F2 fit test ")

                                #
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param)
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                            fit_fourier = 0
                            if( fit_fourier ):
                                
                                os.chdir(struct_dir)
                                
                                tor_type = "four"
                                n_param = 4
                                
                                param_o = []
                                for n in range(n_param):
                                    param_o.append(0.0)
                                #
                                # Set initial parmeters
                                #
                                param_o[1] =  options.tor_paramo 
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                p_fourier = []
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        p_fourier.append(  0.0 )
                                    #p_fourier[param_ind*n_param +  1 ] =  p_fourierF2[param_ind] 



                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourier[p_indx] )
                                    
                                    print "  F coefficents ",F_coef
                                    
                                                                
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len p_fourier ",len(p_fourier)
                                #
                                #resid = residuals_four(p_fourier,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                                resid = residuals_v5(p_fourier,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourier[p_indx] )
                                    
                                    print "  F coefficents ",F_coef
                                    
                                p_fourier_pre = []
                                for p_indx in range( len(p_fourier)):
                                    p_fourier_pre.append( p_fourier[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    p_fourier,successp = optimize.leastsq(residuals_v5,p_fourier,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(p_fourier)):
                                        d_pf = p_fourier[p_indx] - p_fourier_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    if( d_param <  0.0002 ):
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        delta_param = False
                                    
                                    p_fourier_pre = []
                                    for p_indx in range( len(p_fourier)):
                                        p_fourier_pre.append( p_fourier[p_indx] )
                                        
                                resid = residuals_v5(p_fourier,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  round(p_fourier[p_indx]*EVTOKCAL,4) )
                                    
                                    print "  F coefficents ",F_coef
                                    #print "  F coefficents ",F_coef[0]*EVTOKCAL,F_coef[1]*EVTOKCAL,F_coef[2]*EVTOKCAL,F_coef[3]*EVTOKCAL
                                    
                                    
                                DIH_PARAM = p_fourier                                    
                                
                                rmsd_en,rmsd_trans = calc_rms(DIH_PARAM,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                #sys.exit(" F2 fit test ")

                                #
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param )
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                            fit_fourier_F2 = 0
                            if( fit_fourier_F2 ):
                                # Intially fit F2 only to provide a good initial guess
                                
                                tor_type = "four_F2"
                                n_param = 1
                                
                                param_o = []
                                for n in range(n_param):
                                    param_o.append(0.0)
                                    
                                #
                                # Set initial parmeters
                                #
                                param_o[0] =  options.tor_paramo 
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                p_fourierF2 = []
                                #
                                for iparam_indnd in range(n_dihtypes):
                                    for val in param_o:
                                        p_fourierF2.append( val )
                                        
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  param_o",len(param_o)
                                print "  len p_fourierF2 ",len(p_fourierF2)
                                #
                                resid = residuals_fourF2(p_fourierF2,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourierF2[p_indx] )
                                    
                                    print "  F2 coefficents "
                                    print "                        F2 ",float(F_coef[0])," eV ",float(F_coef[0])*EVTOKCAL," kcal/mol "
                            
                                p_fourierF2,successp = optimize.leastsq(residuals_fourF2,p_fourierF2,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param),epsfcn=0.0001)
                            
    
                                resid = residuals_fourF2(p_fourierF2,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourierF2[p_indx] )
                                    
                                    print "  F2 coefficents "
                                    print "                        F2 ",float(F_coef[0])," eV ",float(F_coef[0])*EVTOKCAL," kcal/mol "
                                    
                                DIH_PARAM = p_fourierF2                                    
                                    
                                #sys.exit(" F2 fit test ")

                                os.chdir(struct_dir)
                                
                                
                                    
                                # Fitted torsional energies
                                #fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param )
                                
                                # Plot tossional potential
                                #if( options.plot_tor ):
                                #    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                #    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                                
                                tor_type = "four"
                                n_param = 4
                                
                                param_o = []
                                for n in range(n_param):
                                    param_o.append(0.0)
                                #
                                # Set initial parmeters
                                #
                                param_o[1] =  options.tor_paramo 
                                #
                                # set parameter list consisting of each different dihedral's parameters
                                #
                                p_fourier = []
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        p_fourier.append(  0.0 )
                                    p_fourier[param_ind*n_param +  1 ] =  p_fourierF2[param_ind] 



                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourier[p_indx] )
                                    
                                    print "  F coefficents ",F_coef
                                    
                                                                
                                
                                print "  n_dihtypes ",n_dihtypes
                                print "  len p_fourier ",len(p_fourier)
                                #
                                resid = residuals_v5(p_fourier,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                #resid = residuals_four(p_fourier,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param)
                                
                                error = 0.0
                                for er in resid:
                                    error += er*er
                                print " Initial error ",error
                                #sys.exit( "DIH_PARAM") 
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  p_fourier[p_indx] )
                                    
                                    print "  F coefficents ",F_coef
                                    
                                p_fourier_pre = []
                                for p_indx in range( len(p_fourier)):
                                    p_fourier_pre.append( p_fourier[p_indx] )
                                    
                                delta_param = True
                                fit_iter = 0 
                                while delta_param:
                                    fit_iter +=  1 
                                    #p_fourier,successp = optimize.leastsq(residuals_four,p_fourier,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param),epsfcn=0.0001)                                
                                    p_fourier,successp = optimize.leastsq(residuals_v5,p_fourier,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs),epsfcn=0.0001)                                
                                    
                                    d_param = 0.0 
                                    for p_indx in range( len(p_fourier)):
                                        d_pf = p_fourier[p_indx] - p_fourier_pre[p_indx]
                                        d_param += numpy.sqrt(d_pf*d_pf )
                                    
                                    if( d_param <  0.0002 ):
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        delta_param = False
                                    
                                    p_fourier_pre = []
                                    for p_indx in range( len(p_fourier)):
                                        p_fourier_pre.append( p_fourier[p_indx] )
                                        
                                resid = residuals_v5(p_fourier,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                
                                
                                error = 0.0 
                                for er in resid:
                                    error += er*er
                                print " final error ",error
                                
                                                    
                                for param_ind in range( n_dihtypes ):
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append(  round(p_fourier[p_indx]*EVTOKCAL,4) )
                                    
                                    print "  F coefficents ",F_coef
                                    #print "  F coefficents ",F_coef[0]*EVTOKCAL,F_coef[1]*EVTOKCAL,F_coef[2]*EVTOKCAL,F_coef[3]*EVTOKCAL
                                    
                                    
                                DIH_PARAM = p_fourier                                    
                                rmsd_en,rmsd_trans = calc_rms(DIH_PARAM,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef,tor_type,n_param,min_indx,max_indx,trans_list,trans_indxs)
                                    
                                #sys.exit(" F2 fit test ")

                                #os.chdir(struct_dir)
                                
                                #sys.exit(" Fourier fitting 2")
                                
                                    
                                # Fitted torsional energies
                                fitted_ptor , fit_toten , fitted_comp = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,tor_type,n_param )
                                
                                # Plot tossional potential
                                if( options.plot_tor ):
                                    plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                                    plot_tor_comp(options,DIH_PARAM,ANG_IND,FF_DIHTYPES, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten , fitted_comp )
                                
                                
                            
                            
                            os.chdir(work_dir)
                                
                            
                            # Print parameters for itp file
                            #  the factor of 2.0 will be removed once the parameters are normalized
                            
                            print "   Gromacs format x2 for compatability with biaryl dih parameters "
                            for a_indx in range( len(ang_types) ):
                                param_ind = ANG_IND[a_indx]
                                dih_id = ""
                                for dih_indx in FF_DIHTYPES[param_ind]:
                                    dih_id += str( dih_indx ) + "  " 
                            
                                if( tor_type == "RB" ):
                                        
                                    F_coef = ''
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef +=  str(round( DIH_PARAM[p_indx]*EVTOkJ*2.0 ,4) ) + "  "
                                        
                                    print dih_id,"  3  ",F_coef, "   0.000000  # fit 1 "
                                
                                elif(   tor_type == "four" ):
                                        
                                    F_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef.append( round( DIH_PARAM[p_indx]*EVTOkJ*2.0 ,4) )
                                        
                                    C0 = F_coef[1]  + 0.5*( F_coef[0] + F_coef[2] )
                                    C1 =0.5*( -1.0*F_coef[0] + 3.0*F_coef[2] )
                                    C2 = -1.0*F_coef[1] + 4.0*F_coef[3]
                                    C3 = -2.0*F_coef[2]
                                    C4 = -4.0*F_coef[3]
                                    C5 = 0.0
                                    print dih_id,"  3  ",C0,C1,C2,C3,C4,C5, "  # fit 1 "
                        
                            print "  Lammps format  "
                            for a_indx in range( len(ang_types) ):
                                param_ind = ANG_IND[a_indx]
                                dih_id = ""
                                for dih_indx in FF_DIHTYPES[param_ind]:
                                    dih_id += str( dih_indx ) + "  "
                                
                                if( tor_type == "RB" ):
                                        
                                    RB_coef = []
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        RB_coef.append( DIH_PARAM[p_indx]*EVTOKCAL )
                                    
                                                
                                    F1 = -1.0*( 2.0*RB_coef[1] + 3.0*RB_coef[3]/2.0)
                                    F2 = -1.0*( RB_coef[2] + RB_coef[4])
                                    F3 = -0.5*RB_coef[3]
                                    F4 = -0.25*RB_coef[4]
                    
                                    print "  1  ",F1,F2,F3,F4, "  #  ",dih_id
                                
                                elif(   tor_type == "four" ):
                                        
                                    F_coef = ''
                                    for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                        F_coef +=  str(round( DIH_PARAM[p_indx]*EVTOKCAL ,4) ) + "  "
                                    
                                    print "  1  ",F_coef, "  #  ",dih_id
                                    
                            
                                    
                                
                        
                                
                            
                        else:
                            print " compopatable data files ",dih_qm," and ",dih_ff
                            
                        
    
if __name__=="__main__":
    main() 
