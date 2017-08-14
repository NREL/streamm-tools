#! /usr/bin/env python
"""
Handy routines to analize sets of values 

NoteTK This should be changed to gen.py 
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

import numpy as np
import scipy
import csv , math 


def round_sigfigs(num, sig_figs):
    """
    Round to specified number of significant figures.
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0


def calc_dec(num):
    """
    Round to specified number of significant figures.
    """
    return -int(math.floor(math.log10(abs(num)))) 

class Values(object):
    '''
    Functions for data processing of sets of values 
    '''

    def __init__(self, tag, data_raw, verbose=False):
        '''
        Constructor for Value object 
        '''
        self.tag = tag
        self.data_raw = data_raw

        self.bin_size = 1.0

        self.calc_stats()
        self.set_bins(self.bin_size)


    def __del__(self, verbose=False):
        """
        Deconstructor for a  Value object 
        """
        del self.tag
        del self.data_raw 

    def set_bins(self,bin_size):
        '''
        Set bin related values based on bin_size
        '''
        logger.debug("Setting bin related values based on bin size of {}".format(bin_size))
        
        self.bin_size = bin_size

        self.val_min = min(self.data_raw)
        self.val_max = max(self.data_raw)

        self.val_floor = int(self.val_min/self.bin_size)*self.bin_size - 1.0*self.bin_size
        self.val_ceil = int(self.val_max/self.bin_size)*self.bin_size + 1.0*self.bin_size
        self.val_range = self.val_ceil - self.val_floor
        self.n_bins = int(self.val_range/self.bin_size) #+ 1
        #self.bin_w = (self.val_ceil-self.val_floor)/float(self.n_bins)

    def calc_stats(self):
        '''
        Calculate statistical values using numpy 
        '''
        self.n = len(self.data_raw)
        self.max = max(self.data_raw)
        self.min = min(self.data_raw)
        self.ave = np.average(self.data_raw)
        self.mean = np.mean(self.data_raw)
        self.std = np.std(self.data_raw)
        self.mean_std = self.std/np.sqrt(self.n-1)
        # self.skew = scipy.stats.skew(self.data_raw)
        # self.conf_int = scipy.stats.norm.interval(0.95,loc=self.mean,scale=self.std)
        # self.conf_min, self.conf_max = self.conf_int
        
        # Error is  given by the standard deviation of the mean with the appropriate two-tailed
        #   '''Students t-distribution prefactor to give a 95 percent confidence interval'''
        self.error = 1.96*self.mean_std
        self.error_sf = round_sigfigs(self.error,1)
        error_dec = calc_dec(self.error_sf)
        self.ave_sf = round(self.ave,error_dec)
        self.std_sf = round(self.std,error_dec)
        
        # self.error = self.conf_max-self.conf_min

    def calc_hist(self,setden = True):
        '''
        Create histogram using numpy 
        '''
        self.hist_cent,self.bins = np.histogram(self.data_raw, bins=self.n_bins, range=(self.val_floor,self.val_ceil) , density=setden)
        #
        # The center value of each bin with the correct length... 
        self.bins_cent = (self.bins[:-1] + self.bins[1:]) / 2
        # Add zero to initialize step list to be used
        # ax.plot(self.bins_steps,self.hist_steps,'k*-' ,ls = "steps")
        self.hist_steps = []
        self.hist_steps.append(0.0)
        for v in self.hist_cent:
            self.hist_steps.append(v)
        self.hist_steps.append(0.0)

        self.bins_steps = []
        for v in self.bins:
            self.bins_steps.append(v)
        self.bins_steps.append(self.bins_steps[-1]+self.bin_size)
        

    def write_hist(self):
        '''
        Write histogram to data file 
        '''
        hist_file = "{}.hist".format(self.tag)
        print " Writing histogram {}".format(hist_file)
        with open(hist_file, 'wb') as struc_fout:
            hist_writer = csv.writer(struc_fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            hist_writer.writerow( ('bins_cent','hist') )
            for i in range(len(self.hist)):
                b_i = self.bins_cent[i]
                h_i = self.hist[i]
                hist_writer.writerow( [b_i, h_i])
                
        
    def write_hist_v1(self):
        '''
        Depricated version of writing histogram file. 
        '''
        hist_file = "{}.hist".format(self.tag)
        print " Calculating histogram {}".format(hist_file)

        self.hist_line = "# Ave %f  \n"%(self.ave)
        self.hist_line += "# Standard deviation  %f \n"%(self.std) #,self.sigma)
        #self.hist_line += "# Standard deviation of the mean %f  with a  90 percent confidence interval  \n"%(self.sigma_mean)
        self.hist_line += "# Center of value bin , count \n"
        self.bin_val = [] #np.array()
        self.hist_bin = [] #np.array()
        include_zeros = True 
        for i in range (self.n_bins):
            if( self.np_hist_bin_all[i] > 0.0 or include_zeros ):
                self.bin_val.append(self.np_hist_val_all[i])
                self.hist_bin.append(self.np_hist_bin_all[i])
                self.hist_line += " %f %e  \n"%(self.np_hist_val_all[i],self.np_hist_bin_all[i])
        self.bin_val.append(self.np_hist_val_all[-1] + self.bin_size )
        self.hist_bin.append(0.0)
                
        self.np_bin_val = np.array(self.bin_val)
        self.np_hist_bin = np.array(self.hist_bin)

        self.values = np.array(self.bin_val)
        self.prob_dens = np.array(self.hist_bin)

        #np_hist_val = (np_bin_val[:-1] + np_bin_val[1:]) / 2
        hist_out = open(hist_file,"w")
        hist_out.write(self.hist_line)
        hist_out.close()

        # print hist_line

        #return np_bin_val , np_hist_bin 
        self.hist_val_cent = self.np_hist_val_all + self.bin_size/2.0
        #return hist_val_cent , np_hist_bin_all 
        
    def hist_v1(self):
        
        setden = True             
        hist_file = "{}.hist".format(self.tag)
        
        print " Calculating histogram {}".format(hist_file)
                
        self.np_hist_bin_all,self.np_bin_val_all = np.histogram(self.data_raw, bins=self.n_bins+1, range=(self.val_floor-self.bin_w/2,self.val_ceil+self.bin_w/2) , density=setden)

        
        # The center value of each bin with the correct length... 
        self.np_hist_val_all = (self.np_bin_val_all[:-1] + self.np_bin_val_all[1:]) / 2

        self.hist_line = "# Ave %f  \n"%(self.ave)
        self.hist_line += "# Standard deviation  %f \n"%(self.std) #,self.sigma)
        #self.hist_line += "# Standard deviation of the mean %f  with a  90 percent confidence interval  \n"%(self.sigma_mean)
        self.hist_line += "# Center of value bin , count \n"
        self.bin_val = [] #np.array()
        self.hist_bin = [] #np.array()
        include_zeros = True 
        for i in range (self.n_bins):
            if( self.np_hist_bin_all[i] > 0.0 or include_zeros ):
                self.bin_val.append(self.np_hist_val_all[i])
                self.hist_bin.append(self.np_hist_bin_all[i])
                self.hist_line += " %f %e  \n"%(self.np_hist_val_all[i],self.np_hist_bin_all[i])
        self.bin_val.append(self.np_hist_val_all[-1] + self.bin_size )
        self.hist_bin.append(0.0)
                
        self.np_bin_val = np.array(self.bin_val)
        self.np_hist_bin = np.array(self.hist_bin)

        self.values = np.array(self.bin_val)
        self.prob_dens = np.array(self.hist_bin)

        #np_hist_val = (np_bin_val[:-1] + np_bin_val[1:]) / 2
        hist_out = open(hist_file,"w")
        hist_out.write(self.hist_line)
        hist_out.close()

        # print hist_line

        #return np_bin_val , np_hist_bin 
        self.hist_val_cent = self.np_hist_val_all + self.bin_size/2.0
        #return hist_val_cent , np_hist_bin_all 


    def read_hist(self,dat_file):
        """
        Read in a column from a file

        dat_file - file name 
        dat_col - column number

        NoteTK:
        
        This should be a json/csv 
        """
        self.values = []
        self.prob_dens = []
        self.ave = 0.0 
        self.sigma = 0.0 
        self.error = 0.0 
        try:
            with open(dat_file) as f:
                
                dat_raw = [] 
                F = open(dat_file  , 'r' )
                lines = F.readlines()
                F.close()
                print " Reading ",dat_file
                #
                # Read in data header with number of parameters 
                #
                for line in lines:
                    col = line.split()
                    if( len(col) >= 2 ):

                        if ( col[0] != "#"  ):
                            #if( float(col[1]) > 0.0 ):
                            self.values.append(float(col[0]))
                            self.prob_dens.append(float(col[1]))
                        else:
                            if( col[1] == "Ave"  ):
                                self.ave = float(col[2])
                            if( col[1] == "Standard" and  col[2] == "deviation" and  col[3] != "of" ):
                                self.sigma = float(col[3])
                            if( col[1] == "Standard" and  col[2] == "deviation" and  col[3] == "of" ):
                                self.error = float(col[6])
                                                

                self.bin_size = self.values[1] - self.values[0]
                self.values.append(self.values[-1] +  self.bin_size )
                self.prob_dens.append(0.0)

                self.values = np.array(self.values)
                self.prob_dens = np.array(self.prob_dens)


                # print "self.bin_size ",self.bin_size
                

                return 
            
        except IOError:
            logger.warning(" File %s not found "%(dat_file))
            

        
