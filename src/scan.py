"""
Potential energy scan object
"""

import numpy as np
import  sys
import collections
from scipy import optimize
import math
#import matplotlib.pyplot as plt

import pbcs


def FourierSeries(p,theta):
    """
    Evaluate Fourier Series 
    """
    # theta dihedral angle
    
    debug = 0
    theta_rad = math.radians( theta )
    fourier_sum = 0.0
    
    for n in range(len(p)):
        fourier_sum +=  p[n]*( math.cos( float(n)*theta_rad  ) )
        if( debug ):
            print n,p[n],theta_rad, math.cos( float(n)*theta_rad ), p[n]*( math.cos( float(n)*theta_rad  ) )
    
    return fourier_sum


def d_FourierSeries(p,theta):
    """
    Evaluate first and second derivatives of Fourier Series 
    """
    # theta dihedral angle
    
    debug = 0
    theta_rad = math.radians( theta )
    d_fourier_sum = 0.0
    dd_fourier_sum = 0.0

    for n in range(1,len(p)):
        d_fourier_sum +=  -1*p[n]*( math.sin( float(n)*theta_rad  ) )
        dd_fourier_sum +=  -1*p[n]*( math.cos( float(n)*theta_rad  ) )
        
    
    return round(d_fourier_sum,8),round(dd_fourier_sum,8)

def mmt_FourierSeries(p,range_coord):
    """
    Find maximum minimum and transistions for Fourier sieries for a given range
    """
    
    debug = False   
    
    fourier_sum = 0.0

    Fourier_max = []
    Fourier_min = []
    Fourier_trans = []

    # Test first point
    calc_minmaxtrans = 0 
    calc_maxmintrans = 0
    
    #print " ragne ", range_coord[0], range_coord[2], range_coord[1]
    range_coord[1] = 1.0
    step = 0.01
    r_mult = int(1/step)
    
    check_min = 1
    check_max = 1

    for coord_x in range( 0*r_mult,181*r_mult ):
        coord_i = float(coord_x)*step
        coord_m = coord_i - step
        coord_m_m = coord_m - step
        coord_p = coord_i + step
        coord_p_p = coord_p + step

        

        d_i,dd_i =  d_FourierSeries(p, coord_i )
        d_m,dd_m = d_FourierSeries(p, coord_m )
        d_m_m,dd_m_m = d_FourierSeries(p, coord_m_m )
        d_p,dd_p = d_FourierSeries(p, coord_p )
        d_p_p,dd_p_p = d_FourierSeries(p, coord_p_p )

        val_i = FourierSeries(p,coord_i)
        val_m = FourierSeries(p,coord_m)
        val_m_m = FourierSeries(p,coord_m_m)
        val_p = FourierSeries(p,coord_p)
        val_p_p = FourierSeries(p,coord_p_p)

        d_m = pbcs.d_2p( val_m,val_i,step) 
        d_m_m = pbcs.d_2p( val_m_m,val_m,step) 
        d_p = pbcs.d_2p( val_i,val_p,step) 
        d_p_p = pbcs.d_2p( val_p,val_p_p,step) 


        if( debug):
            print "coord_i",coord_i
            print "val_i",val_i
            print "coord_m",coord_m
            print "coord_m_m",coord_m_m
            print "coord_p",coord_p
            print "coord_p_p",coord_p_p
            print "d_m",d_m
            print "d_m_m",d_m_m
            print "d_p",d_p
            print "d_p_p",d_p_p

        check_min += 1 
        check_max += 1 
        
        if( d_m_m < 0 and d_m < 0  and d_p > 0  and d_p_p > 0 and check_min > 1  ):
            
            Fourier_min.append( coord_i  )
            min_en_i = val_i
            calc_minmaxtrans = 1

            if( calc_maxmintrans ):
                trans_ev = max_en_i - min_en_i
                Fourier_trans.append(  trans_ev )

                calc_maxmintrans = 0 

                if( debug):
                    print " calc_maxmintrans found ",trans_ev

            success = 1
            check_min = 0

            if( debug):
                print " min found ",coord_i


        if( d_m_m > 0 and d_m > 0  and d_p < 0  and d_p_p < 0 and check_max > 1  ):
            Fourier_max.append( coord_i  )
            max_en_i = val_i
            calc_maxmintrans = 1

            if( debug):
                print " max found ",coord_i

            if( calc_minmaxtrans ):
                trans_ev = max_en_i - min_en_i
                Fourier_trans.append(  trans_ev )
                calc_minmaxtrans = 0 

                if( debug):
                    print " calc_minmaxtrans found ",trans_ev

            success = 1
            check_max = 0

        

                                   
    return Fourier_max,Fourier_min,Fourier_trans


def residuals_FS(Fourier_coef, coord_list, targets, wt_list,wt_coef ):
    """
    fourier series residuals 
    """
    
    debug = 0
    
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(Fourier_coef) ):
        rounded_param = round(Fourier_coef[param_i] ,8 )
        Fourier_coef[param_i]   = rounded_param 
    
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        coord_val = coord_list[cent_indx] #[0]
        wt = wt_list[cent_indx]

        tor_en = FourierSeries(Fourier_coef,coord_val)
        
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        resid.append(sq_delta)
        
    return resid

class PES():
    """
    Potential energy surface
    """

    def __init__(self,method,coord_list,energy_list,method_list,basis_list,tag_list, verbose=False):
        """
        Constructor: for pes class

        Args:
            method (str) method used to calculate pes 
            verbose (bool): flag for printing status/debug info
        """
        self.verbose=verbose        
        self.method = method
        self.coord_list = np.array(coord_list)        # Coordinate point of energy calculation 
        self.energy_list = np.array(energy_list)    # Energy at coordinate 
        self.method_list = method_list
        self.basis_list = basis_list
        self.tag_list = tag_list

        self.min_energy_list = self.calc_min_energy_list()     # shifted to min
        self.first_energy_list = self.calc_first_energy_list()      # shifted to zero coord

        
        
        # List of minimums, maximums and transitions 
        self.min_indx = []
        self.max_indx = []
        self.trans_list = []
        self.trans_indxs = []
        self.k_list = []
        # Weight list 
        self.wt_coef = 1.0
        self.weights = []
        # Error
        self.error_list = []
        # Coordinate information 
        self.step_coord = self.coord_list[1] - self.coord_list[0]
        self.max_coord = max(self.coord_list)
        self.min_coord = min(self.coord_list)
        # Extra information 
        self.tag = []


    def set_tag(self,tag):
        """
        Set a tag for the PES 
        """
        self.tag = tag
        
    def get_tag(self):
        """
        Get a tag for the PES 
        """
        return self.tag 
        

    def calc_min_energy_list(self):
        """
        Calculate the pes shifted above zero
        """
        return self.energy_list - self.min_energy()

    def calc_first_energy_list(self):
        """
        Calculate the pes shifted to first point 
        """
        return self.energy_list - self.energy_list[0]
        
    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        del self.method
        del self.energy_list
        del self.coord_list
        del self.min_energy_list
        del self.first_energy_list 
        del self.method_list 
        del self.basis_list 
        del self.tag_list

        del self.min_indx 
        del self.max_indx 
        del self.trans_list 
        del self.trans_indxs 
        del self.k_list 
        # Weight list 
        del self.wt_coef 
        del self.weights 
        # Error
        del self.error_list

        del self.step_coord 
        del self.max_coord 
        del self.min_coord 
        # Extra information 
        del self.tag
        

        
    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.energy_list)

    def min_energy(self):
        """
        Find minimum energy
        """
        return min(self.energy_list)
    

    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        str_line = " # coordinate ; energy \n"
        for indx in range(len(self.energy_list)):
            str_line += " %f %f  %f %f \n"%(self.coord_list[indx],self.energy_list[indx],self.min_energy_list[indx],self.first_energy_list[indx] )
        return str_line

    def get_coord(self):
        """
        Return coord list
        """
        return self.coord_list

    def get_energy(self):
        """
        Return energy list
        """
        return self.energy_list


    def get_min_energy(self):
        """
        Return min energy list
        """
        return self.min_energy_list

   
    def overwrt_energy_w_min(self):
        """
        Over write energy list with min_energy_list
        """
        self.energy_list = self.min_energy_list


    def write(self,out_pes):
        """
        Write file of pes
        """
        
        str_line = " # coordinate ; energy ; method ; basis ; tag  \n"
        for indx in range(len(self.energy_list)):
            str_line += "%d %f %f  %s %s %s \n"%(indx,self.coord_list[indx],self.energy_list[indx],self.method_list[indx],self.basis_list[indx],self.tag_list[indx] )
        F = open(out_pes,"w")
        F.write(str_line)
        F.close()


    def calc_d2(self):
        """
        Find max and mins based on the second derivative
        """
        debug = False  
        success = 0 
        #
        min_val =  1e16 
        qm_max = -1e16
        ang_min = "nan"
        ang_max = "nan"
        #

        if( debug ):
            for val_i in self.energy_list:
                print val_i 
        
        min_val = min( self.energy_list )
        h = self.coord_list[1] - self.coord_list[0]

        if( debug):
            print "  Mim ",min_val, " Delta coord ",h
            print " length ",len(self.energy_list)  

        # Find inversion points in energy
        self.min_indx = []
        self.max_indx = []
        self.min_val = []
        self.max_val = []
        self.trans_list = []
        self.trans_indxs = []
        self.k_list = []


        # Test first point
        calc_minmaxtrans = 0 
        calc_maxmintrans = 0
        # loop over subsequent points 
        for indx_i in range( len(self.energy_list) ):
            indx_m = indx_i - 1  # _m minus 
            indx_m_m = indx_m - 1 
            indx_p = indx_i + 1  # _p plus 
            indx_p_p = indx_p + 1 

            # apply boundry conditions 
            if( indx_m < 0 ): indx_m  =  -1*indx_m 
            if( indx_m_m < 0 ): indx_m_m  = -1*indx_m_m 

            if(debug): print " position index -2 -1 0 1 2 : ",indx_m_m,indx_m,indx_i,indx_p,indx_p_p

            if( indx_p > len(self.energy_list)  -1  ): indx_p  = indx_p - ( indx_p -  len(self.energy_list) + 1 )*2  # as 0.0 == 180.0 
            if( indx_p_p > len(self.energy_list)  -1  ): indx_p_p  = indx_p_p -  ( indx_p_p -  len(self.energy_list) + 1 )*2  #  as 0.0 == 180.0 

            d_en_m   = pbcs.d_2p(self.energy_list[indx_m] ,self.energy_list[indx_i],h) 
            d_en_m_m = pbcs.d_2p(self.energy_list[indx_m_m] ,self.energy_list[indx_m],h) 
            d_en_p   = pbcs.d_2p(self.energy_list[indx_i] ,self.energy_list[indx_p],h)
            d_en_p_p = pbcs.d_2p(self.energy_list[indx_p] ,self.energy_list[indx_p_p],h)

            if(debug):
                print "  m i p ",self.energy_list[indx_m_m] - min_val,self.energy_list[indx_m] - min_val,self.energy_list[indx_i] - min_val,self.energy_list[indx_p] - min_val,self.energy_list[indx_p_p] - min_val
                print "     dm dp ",d_en_m_m,d_en_m,d_en_p,d_en_p_p
                #print "     dm dp ",d_en_m,d_en_p


            if( d_en_m_m < 0 and d_en_m < 0  and d_en_p > 0  and d_en_p_p > 0 ):
                self.min_indx.append( indx_i )
                min_en_i = self.energy_list[indx_i]
                self.min_val.append( self.coord_list[indx_i] )
                calc_minmaxtrans = 1
                min_transindx = indx_i

                if(debug): print " min found ",self.energy_list[indx_i] - min_val

                if( calc_maxmintrans ):
                    trans_ev = max_en_i - min_en_i
                    self.trans_list.append(  trans_ev )

                    calc_maxmintrans = 0 
                    self.trans_indxs.append( [max_transindx,indx_i] )

                    if(debug): print "   Found  max min trans ",max_en_i - min_val," -> ",min_en_i - min_val," trans = ",trans_ev

                d2_en_i = pbcs.d2_3p(self.energy_list[indx_m] ,self.energy_list[indx_i],self.energy_list[indx_p],h)
                self.k_list.append(  d2_en_i ) 
                if(debug): print " min found ",self.energy_list[indx_i], d2_en_i
                print " min found ",self.energy_list[indx_i], d2_en_i

                success = 1


            if( d_en_m_m > 0 and d_en_m > 0  and d_en_p < 0  and d_en_p_p < 0 ):
            #if( d_en_m > 0  and d_en_p < 0   ):
                self.max_indx.append( indx_i )
                max_en_i = self.energy_list[indx_i]
                self.max_val.append( self.coord_list[indx_i] )
                calc_maxmintrans = 1
                max_transindx = indx_i

                if(debug): print " max found ",self.energy_list[indx_i] - min_val

                if( calc_minmaxtrans ):
                    trans_ev = max_en_i - min_en_i
                    self.trans_list.append(  trans_ev )
                    calc_minmaxtrans = 0 
                    self.trans_indxs.append( [indx_i,min_transindx] )


                    if(debug): print "   Found min max trans ",min_en_i - min_val," -> ",max_en_i - min_val," trans = ",trans_ev

                success = 1


        verbose = True 
        for inv_indx in range( len(self.min_indx) ):
            indx = self.min_indx[inv_indx]
            if(verbose): print "  Min ",inv_indx," found at ",indx*h," w energy ",self.energy_list[indx]


        for inv_indx in range( len(self.max_indx) ):
            indx = self.max_indx[inv_indx]
            if(verbose): print "  Max ",inv_indx," found at ",indx*h," w energy ",self.energy_list[indx]

        for inv_indx in range( len(self.trans_list) ):
            bar = self.trans_list[inv_indx]
            print "  Transion ",inv_indx," =  ",bar," eV"

        if(debug):    sys.exit(" inversion testing ")



        # Minimum
        if( debug ):
            for indx in range( len(self.min_indx) ):
                cent_indx =  self.min_indx[indx]
                if( cent_indx > 0 ):
                    print " Minimum target_en_s meV ",energy_list[cent_indx-1]*1000 , energy_list[cent_indx]*1000 ,  energy_list[cent_indx+1]*1000

        return self.max_val,self.min_val,self.trans_list,success

        ## return (success,  numeric_func,min_indx,max_indx,trans_list,trans_indxs,k_list )


    def set_weights(self):
        """
        Set weights for fitting 
        """
        
        debug = False 

        wt_o = 10.0 
        wt_max = 1000.0 
        wt_min = 1000.0 
        wt_min_pm = 30.0
        wt_min_pm2 = 20.0 

        wt_coef = 0.0

        # Initialize wights 
        self.weights = []
        for cent_indx in range(len(self.energy_list)+1):
            self.weights.append(wt_o)

        for indx in range( len(self.min_indx) ):

            cent_indx =  self.min_indx[indx]
            self.weights[cent_indx] = wt_min
            if( cent_indx < len(self.coord_list)  ):  self.weights[cent_indx + 1 ] = wt_min_pm
            if( cent_indx < len(self.coord_list) - 1 ):  self.weights[cent_indx + 2 ] = wt_min_pm2
            if( cent_indx > 0 ):  self.weights[cent_indx - 1 ] = wt_min_pm
            if( cent_indx > 1 ):  self.weights[cent_indx - 2 ] = wt_min_pm2

            if(debug): print "  Setting min at ",self.coord_list[cent_indx]," to ",wt_min #," with dE = ",k_list[indx]

        for indx in range( len(self.max_indx) ):

            cent_indx =  self.max_indx[indx]
            self.weights[cent_indx] = wt_max

        debug = False 
        if( debug ): 
            for cent_indx in range(len(self.energy_list)):
                print  self.coord_list[ cent_indx ], self.energy_list[ cent_indx ], self.weights[cent_indx]
            sys.exit(" debug weights ")

        self.wt_coef = 1.0

        # return  ( self.weights, wt_coef)

    def fit_fourier(self):
        """
        Fit fourier series 
        """
        verbose = True 

        fourier_order = 6
        Fourier_coef_fit = []
        param_list_pre = []
        for n in range(fourier_order+1):
            Fourier_coef_fit.append( 0.0 )
            param_list_pre.append(  0.0 )
        
        print " fitting ",len(Fourier_coef_fit)," parameters to data ", len(self.energy_list)
        delta_param = True
        fit_iter = 0
        while delta_param:
            fit_iter +=  1

            Fourier_coef_fit,success = optimize.leastsq(residuals_FS,Fourier_coef_fit,args=(self.coord_list, self.energy_list, self.weights, self.wt_coef ),epsfcn=0.0001)

            d_param = 0.0 
            for p_indx in range( len(Fourier_coef_fit)):
                d_pf = Fourier_coef_fit[p_indx] - param_list_pre[p_indx]
                d_param += np.sqrt(d_pf*d_pf )

            print "   fit_iter ",fit_iter," with delat_param ",d_param
            if( d_param <  0.0001 ):
                delta_param = False

            param_list_pre = []
            for p_indx in range( len(Fourier_coef_fit)):
                param_list_pre.append( Fourier_coef_fit[p_indx] )

        if( verbose ):
            log_line = " Fourier fit finished with %d iterations "%(fit_iter)
            log_line += " Coef : %s "%(str(Fourier_coef_fit))
            print log_line
                                                
        return Fourier_coef_fit,success

    

    def rmse_FS(self,Fourier_coef):
        """
        Rout mean squared error fourier series 
        """
        
        debug = 0

        resid = []
        fit_en = []
        self.error_list = []
        for cent_indx in range(len(self.energy_list)):
            t_en = self.energy_list[cent_indx]
            coord_val = self.coord_list[cent_indx] #[0]
            tor_en = FourierSeries(Fourier_coef,coord_val)
            fit_en.append( tor_en )
            delta_en =  t_en - tor_en
            self.error_list.append(delta_en)
            sq_delta =  delta_en*delta_en

            resid.append(sq_delta)


        rmsa = sum( resid)/float( len(self.energy_list) )
        rmse = math.sqrt( rmsa )

        return rmse,fit_en,resid

    def plot_fit(self,fit_en,resid):
        """
        Rout mean squared error fourier series 
        """
        
        debug = 0

        # plt.ylabel('Energy (kcal/mol)')
        plt.ylabel('Energy (eV)')
        plt.xlabel('coordinate')

        plt.plot(self.coord_list,self.energy_list,'kx', label="target" )
        plt.plot(self.coord_list,fit_en,"g--", label="fit" )
        plt.plot( self.coord_list,self.error_list,"b^", label="Error" )
        plt.legend(loc=(0.67,0.72),prop={'size':10})

        plt.savefig("Fit.pdf",format='pdf')


def delta(pes_1,pes_2):
    """
    Take the difference in energy between two pes
    """
    verbose = True

    # Initialize PES lists 
    coord_list = []
    energy_list = []
    method_list = []
    basis_list = []
    tag_list = []

    for indx in range(len(pes_1.energy_list)):
        if( pes_1.coord_list[indx] != pes_2.coord_list[indx] ):
            error_line = " Coordinate %d of PES's %f and %f don't match "%( indx,pes_1.coord_list[indx],pes_2.coord_list[indx] )
            sys.exit(error_line)

        if( verbose ):
            print " %f %f  %f %f \n"%(pes_1.coord_list[indx],pes_1.min_energy_list[indx],pes_2.coord_list[indx],pes_2.min_energy_list[indx] )

        coord_list.append(pes_1.coord_list[indx])
        energy_list.append(pes_1.min_energy_list[indx] - pes_2.min_energy_list[indx])
        method_list.append( "%s-%s"%(pes_1.method_list[indx] ,pes_2.method_list[indx]) )
        basis_list.append("%s-%s"%(pes_1.basis_list[indx] ,pes_2.basis_list[indx]) )
        tag_list.append("%s-%s"%(pes_1.tag_list[indx] ,pes_2.tag_list[indx]) )

    return PES("delta",coord_list,energy_list,method_list,basis_list,tag_list )



def add(pes_1,pes_2):
    """
    Add two PES 
    """
    verbose = True

    # Initialize PES lists 
    coord_list = []
    energy_list = []
    method_list = []
    basis_list = []
    tag_list = []

    for indx in range(len(pes_1.energy_list)):
        if( pes_1.coord_list[indx] != pes_2.coord_list[indx] ):
            error_line = " Coordinate %d of PES's %f and %f don't match "%( indx,pes_1.coord_list[indx],pes_2.coord_list[indx] )
            sys.exit(error_line)

        if( verbose ):
            print " %f %f  %f %f \n"%(pes_1.coord_list[indx],pes_1.min_energy_list[indx],pes_2.coord_list[indx],pes_2.min_energy_list[indx] )

        coord_list.append(pes_1.coord_list[indx])
        energy_list.append(pes_1.min_energy_list[indx] + pes_2.min_energy_list[indx])
        method_list.append( "%s+%s"%(pes_1.method_list[indx] ,pes_2.method_list[indx]) )
        basis_list.append("%s+%s"%(pes_1.basis_list[indx] ,pes_2.basis_list[indx]) )
        tag_list.append("%s+%s"%(pes_1.tag_list[indx] ,pes_2.tag_list[indx]) )

    return PES("delta",coord_list,energy_list,method_list,basis_list,tag_list )


def read(in_pes):

    """
    Read file of pes
    """
    verbose = True 

    coord_list = []
    energy_list = []
    method_list = []
    basis_list = []
    tag_list = []

    F = open(in_pes , 'r' )
    en_lines = F.readlines()
    F.close()
    if( verbose ): print " %s energy file found "%(in_pes)
    for en_line in en_lines:
        col = en_line.split()
        if( col[0] != "#" and len(col) >= 4 ):
            coord_list.append(float(col[1]))
            energy_list.append( float(col[2]))
            method_list.append(col[3])
            basis_list.append(col[4])

            if( len(col) >= 6 ):
                tag = str(col[5:])
            else:
                tag = "NA"

            tag_list.append( tag)

    method = method_list[0]
    
    return PES(method,coord_list,energy_list,method_list,basis_list,tag_list )
