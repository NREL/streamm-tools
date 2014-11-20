"""
Analyize potential energy surface objects (pes)
"""

from scan import PES
import sys 

import scan, pbcs 

import matplotlib.pyplot as plt


def get_options():
    """
    Set options
    """
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    parser.add_option("-f","--fourierfit", dest="fourierfit", default=False,action="store_true", help="Fit fourier components ")

    parser.add_option("--pes", dest="pes", default="",type="string",help=" Potential energy scan data file  ")
    

    (options, args) = parser.parse_args()
        
    return options, args
   

def main():
    """
    Analysis of PES objects 

    """
    options, args = get_options()

    debug = True

    # Read in PES object
    pes = scan.read(options.pes)

    # Over write energy list with min_energy_list for fitting 
    pes.overwrt_energy_w_min()

    if( options.fourierfit ):

        # Find second derivates for max/min
        success = pes.calc_d2()

        # Set wieghts
        pes.set_weights()
        
        # Fit fourier coeff
        Fourier_coef_fit,success = pes.fit_fourier()
        
        #
        rmse_four,fit_en,resid = pes.rmse_FS(Fourier_coef_fit)

        if( options.verbose ):
            log_line = " Coef : %s \n"%(str(Fourier_coef_fit))
            log_line += " RMSE : %f \n"%(rmse_four)
            print log_line

        pes.plot_fit(fit_en,resid)
            
if __name__=="__main__":
    main()
