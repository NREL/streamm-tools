"""
Class data structures for Gaussian data
"""

__author__ = "Travis W. Kemper"
__version__ = "0.1"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Alpha"


import copy, sys, os, shutil, math
import time, datetime
import json
import numpy as np
from string import replace


import structure, parameters, units, periodictable, project, resource, buildingblock
from calculation import CalculationRes

import logging
logger = logging.getLogger(__name__)
          
class NWChem(CalculationRes):
    """
    Dervied class implementing input/output methods Gaussian
    """
    def __init__(self, tag , verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """
        # Base class constructor is called
        CalculationRes.__init__(self, tag)

        self.meta['software'] = 'gaussian'
        self.units['distance'] = 'bohr'
        self.units['energy'] = 'hartree'
        # String found in log file when simulation finishes
        self.properties['finish_str'] = 'Normal termination of Gaussian'
        #
    def __del__(self):
        """
        Destructor, clears object memory
        """
        # Call base class destructor
        CalculationRes.__del__(self)


    def check(self):
        '''
        Check if simulation is finished 
        '''
        self.check_file(finish_key_str = 'Normal termination of Gaussian',key = "log")
        

    def analyze_log(self,log_file):
        """
        Read NWChem simulation log files  
        """
        # verbose = False 
        debug = False

        self.calctype = 'et'

        # self.converged = False
        
        try:
            with open(log_file,"r") as F:
                log_lines = F.readlines()
                F.close()
                # self.converged = True


                self.stinfo = []
                self.N_alpha_occ = 0  
                self.N_beta_occ = 0
                self.alpha_energies = []
                self.beta_energies = []

                # Initialize read functions as off 
                read_et = False 
                read_geom = False
                geom_name = "unknown"
                
                line_cnt = 0
                for line in log_lines:
                    col = line.split()
                    line_cnt += 1
                
        
                    if( read_et and len(col) > 0  ):
                        
                        logger.debug(" et col %d %d %s ",line_cnt,len(col),str(col))
                        
                        if( len(col) >= 5 ):
                            logger.debug(" >=5 et col %d %d %s ",line_cnt,len(col),str(col))
                            
                            if(  col[0]  == "MO" and  col[1]  == "vectors" and  col[3]  == "reactants:" ):
                                et_ij.reactantMO = str(col[4])                                
                            if(  col[0]  == "Reactants/Products" and  col[1]  == "overlap" ):                            
                                et_ij.S = float( replace(str(col[4]),"D","e") )
                                
                        if( len(col) >= 6 ):
                                
                            logger.debug(" >=6 et col %d %d %s ",line_cnt,len(col),str(col))
                                
                            if(  col[0]  == "MO" and  col[1]  == "vectors" and  col[3]  == "products" ):
                                et_ij.productMO = str(col[5])                        
                            if(  col[0]  == "Electronic" and  col[1]  == "energy" and  col[3]  == "reactants" ):
                                et_ij.reactanten = float(col[5])
                            if(  col[0]  == "Electronic" and  col[1]  == "energy" and  col[3]  == "products" ):
                                et_ij.producten = float(col[5])

                            if(  col[0]  == "Electron" and  col[1]  == "Transfer"  and  col[2]  == "Coupling"  and  col[3]  == "Energy" ):                            
                                #et_ij.V = float( replace(str(col[5]),"D","e") )                                
                                et_ij.V = float( col[5] )
                            if(  col[0]  == "Task" and  col[1]  == "times"  and  col[1]  == "times" ):
                                read_et = False
                                et_ij.cputime = col[3]
                                
                                print ">analyze_log ",et_ij
                                
                                self.et_list.append(copy.deepcopy(et_ij))
                                
                    if( len(col) >= 3 and read_et == False ):
                        if(  col[0]  == "Electron" and  col[1]  == "Transfer"  and  col[2]  == "Calculation" ):
                            read_et = True
                            et_ij = electrontransfer()
                            logger.info(" Electron Transfer Calculation found on line {} ".format(line_cnt))


                    if( len(col) > 3 ):
                        if(  col[0]  == "Geometry" and col[1].replace('"', '').strip() == geom_name and col[2] == "->"  ):
                            read_geom = False # True 
                            line_cnt = 0
                    if( line_cnt >= 7 and read_geom and len(col) >= 5 ):
                        # Need to update 
                        ASYMB.append( col[1] )
                        R.append( [float(col[3]),float(col[4]),float(col[5])] )

                    if( read_geom and len(col) < 1 and line_cnt >= 7  ):
                        read_geom = False 


                logger.debug(" ---- {} Electron Transfer Calculations found ".format(len(self.et_list)))
                    

        except IOError:
            logger.warning("Error output file %s not found "%log_file)
            #self.converged = False
                
            
        return
