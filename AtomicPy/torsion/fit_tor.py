#! /usr/bin/env python
# run ab initio torsional potential 

# Dr. Travis Kemper
# NREL
# 12/09/2013
# travis.kemper@nrel.gov



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
    
    parser.add_option("--tor_paramo", type="float", dest="tor_paramo", default=1.0, help="Value to initialize torsional potential coefficients ")
    parser.add_option("--out_ftor", dest="out_ftor", type="string", default="tor_fit.dat", help="Output fitted torsional potential data ")
    parser.add_option("--plot_tor", dest="plot_tor", type="string", default="tor_fit.pdf", help="Output fitted torsional potential graph ")



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

    
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('dihedral angle (deg)')
    
    plt.plot( qm_angle,qm_en_s,'kx', label="mp2" )
    plt.plot( qm_angle,ff_en_s,"g+", label="$V_{tor}$=0" )
    plt.plot( qm_angle,target_en_s,"r^", label="mp2 - ($V_{tor}$=0)" )
    plt.plot( qm_angle,fitted_ptor,"b-", label="$V_{tor}$ fit" )
    plt.plot( qm_angle, fit_toten,"b--", label="$ff_{en}$ fit")
    plt.legend(loc=(0.67,0.72),prop={'size':10})
    
    #plt.show()
    
    plt.savefig(options.plot_tor,format='pdf')


def print_tor(options,DIH_PARAM,ANG_IND,ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s ):
    
    debug = 0
    
    fitted_ptor = []
    fit_toten = []
    
    # Calculate new torsional potential based on previous t(0)
    #
    for cent_indx in range(len(target_en_s)):
        ang_val = ff_angles[cent_indx]
        ff_en = ff_en_s[cent_indx]      
        # Sum torsional energies for each component
        fit_en = float( sum_tor(DIH_PARAM,ang_val,ANG_IND) )
        fitted_ptor.append( fit_en )
        fit_toten.append( fit_en + ff_en )

    fit_en_max = max( fit_toten ) 
    fit_en_min = min( fit_toten ) 
    
    # Write data file 
    f_fit = open(options.out_ftor,"w")
    f_fit.write(" angle ; mp2 ; V_tor ; mp2-Vtor ; V_tor fit ; ff_en fit")
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

        
        f_fit.write( "\n %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f  " % (qm_x,qm_en,ff_en,t_en, fit_en,fit_ten_s))        

    f_fit.close()

    return (fitted_ptor , fit_toten_s)

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

def sum_tor(DIH_PARAM,ang_val,ANG_IND):
    import math
    
    debug = 0
    
    n_param = 4 #len(DIH_PARAM)
    
    tor_lijk = 0.0
    
    for angle_indx in range( len(ang_val)):
        theta = math.radians( ang_val[angle_indx] )
        param_ind = ANG_IND[angle_indx]
        F_coef = []
        
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            F_coef.append(  DIH_PARAM[p_indx] )
        
        if( debug ):
            print " calc pot ",theta,param_ind,param_ind*n_param,n_param,F_coef
        
        V_tor =  V_four(F_coef,theta)
        
        tor_lijk += V_tor
        if( debug ):
            print " V_tor ",V_tor
    
    debug = 0
    if( debug ):
        print "tor_lijk ",tor_lijk
        sys.exit("sum_tor")
    
    return tor_lijk

def residuals(DIH_PARAM,ANG_IND,targets,ff_angles,wt_angle,wt_coef):
    
    debug = 0
 
    n_param = 4 #len(DIH_PARAM)
    
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
        tor_en = sum_tor(DIH_PARAM,ang_val,ANG_IND)
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

    for param_ind in range( NUMB_DIHTYPES  ):
        F_coef = []
        for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
            F_coef.append(  DIH_PARAM[p_indx] )
    
        F1 = F_coef[0]
        F2 = F_coef[1]
        F3 = F_coef[2]
        F4 = F_coef[3]

        sq_coef = wt_coef *( F1*F1 + F2*F2 + F3*F3 + F4*F4 )
        if( debug ):
            print param_ind
            print wt_coef,F_coef[0],F_coef[1],F_coef[2],F_coef[3]
            print sq_coef

        resid.append(sq_coef)

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

    param_o = []
    param_o.append(options.tor_paramo)
    param_o.append(options.tor_paramo)
    param_o.append(options.tor_paramo)
    param_o.append(options.tor_paramo)
    
    # Set initial paratemers
    
    F1=0.00
    F2=1.970035
    F3=0.00
    F4=0.00
    
    param_o = []
    param_o.append(F1)
    param_o.append(F2)
    param_o.append(F3)
    param_o.append(F4)
    

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
            for val in param_o:
                DIH_PARAM.append( val )
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
            

    return ang_types,FF_DIHTYPES, DIH_PARAM, ANG_IND
    
def set_weights(qm_angle,qm_en_s):
    import numpy as np
        
    wt_o = 10.0 
    wt_max = 20.0 
    wt_min = 50.0 
    wt_min_pm = 30.0 
    wt_min_pm2 = 20.0 

    wt_coef = 0.0
    wt_angle = []

    n_anlges =  len(qm_angle)
    
    wt_angle = [] #np.ones(n_anlges+1)
    for cent_indx in range(n_anlges+1):
        wt_angle.append(wt_o)
    
    
    en_minus = qm_en_s[-1] 
    debug = 0
    for cent_indx in range(n_anlges):
        qm_tor_en =  qm_en_s[ cent_indx ]
        # Change wieghts at inversion points
        if( cent_indx > 0 ): en_minus = qm_en_s[ cent_indx - 1 ] 
        if( cent_indx < n_anlges - 1):
            en_plus = qm_en_s[ cent_indx + 1 ]
        else:
            en_plus = qm_en_s[0]

        d_minus = qm_tor_en -  en_minus
        d_plus  = en_plus - qm_tor_en
        dm_sign = 0
        if( d_minus < 0.0 ): dm_sign = -1 
        if( d_minus > 0.0 ): dm_sign =  1
        dp_sign = 0
        if( d_plus < 0.0 ): dp_sign = -1 
        if( d_plus > 0.0 ): dp_sign =  1
            
        if( dm_sign >  dp_sign ): wt_angle[cent_indx] = wt_max
        if( dm_sign <  dp_sign ):
            wt_angle[cent_indx] = wt_min
            if( cent_indx < len(qm_angle)  ):  wt_angle[cent_indx + 1 ] = wt_min_pm
            if( cent_indx < len(qm_angle) - 1 ):  wt_angle[cent_indx + 2 ] = wt_min_pm2
            if( cent_indx > 0 ):  wt_angle[cent_indx - 1 ] = wt_min_pm
            if( cent_indx > 1 ):  wt_angle[cent_indx - 2 ] = wt_min_pm2

    debug = 0
    if( debug ): 
        for cent_indx in range(len(qm_angle)):
            print  qm_angle[ cent_indx ], qm_en_s[ cent_indx ], wt_angle[cent_indx]
    
    wt_coef = 1.0
    
    return ( wt_angle, wt_coef)

    
    
def main():
    import os, sys
    import jsonapy
    import collections
    import math    
    from scipy import optimize
    
    debug = 0 

    
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
                            qm_en_s = []
                            qm_angle = []
                            ff_en_s =  []
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
                                target_en_s.append( target_en[t_indx] - t_en_min )
                                

                            
                            if( options.verbose ):
                                print " Read in ",len(qm_en_s)," energies ",min(qm_en_s),max(qm_en_s)
                            
                            # Read in dihedral information 
                            ang_types , FF_DIHTYPES, DIH_PARAM, ANG_IND =  read_dihtypes(struct_dir,job_name,options)
                            if( options.verbose ):
                                print "  with ",max(ANG_IND)," angle types "
                                for a_indx in range( len(ang_types) ):
                                    ang_i = ANG_IND[a_indx]
                                    print ang_types[a_indx]," type ",ang_i,DIH_PARAM[ang_i],FF_DIHTYPES[ang_i]
                                
                            
                            # Set wieghts
                            wt_angle, wt_coef = set_weights(qm_angle,qm_en_s)
                        
                            #
                            #for i in range(len(target_en_s)):
                            #    print target_en_s[i],wt_angle[i]
                            #sys.exit( "debug")
                        
                            p =  DIH_PARAM
                            resid = residuals(p,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef)
                            error = 0.0
                            
                            for er in resid:
                                error += er*er
                            print " Initial error ",error
                            #sys.exit( "DIH_PARAM") 
                        
                            p_guess =  DIH_PARAM
                            DIH_PARAM,successp = optimize.leastsq(residuals,p_guess,args=(ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef),epsfcn=0.0001)
                        
                            #p=Param(*p)
                            print " fitted ",DIH_PARAM
                        
                            NUMB_DIHTYPES = max( ANG_IND )
                            n_param = 4
                            for param_ind in range( NUMB_DIHTYPES + 1  ):
                                F_coef = []
                                for p_indx in range( param_ind*n_param,param_ind*n_param + n_param):
                                    F_coef.append(  DIH_PARAM[p_indx] )
                                
                                print FF_DIHTYPES[param_ind],F_coef
                            
                        
                            resid = residuals(DIH_PARAM,ANG_IND,target_en_s,ff_angles,wt_angle,wt_coef)
                            error = 0.0 
                            for er in resid:
                                error += er*er
                            print " final error ",error
                            
                            os.chdir(struct_dir)
                            
                            # Fitted torsional energies
                            fitted_ptor , fit_toten  = print_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s )
                            
                            # Plot tossional potential
                            if( options.plot_tor ):
                                plot_tor(options,DIH_PARAM,ANG_IND, ff_angles,ff_en_s,qm_angle,qm_en_s,target_en_s,fitted_ptor , fit_toten )
                            
                            os.chdir(work_dir)
                                
                                
                            
                            
                        
                        
    
if __name__=="__main__":
    main() 
