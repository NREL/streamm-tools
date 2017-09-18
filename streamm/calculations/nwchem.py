# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides, Ross Larsen"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

"""
Class data structures for NWChem data
"""

import copy, sys, os, shutil, math
import time, datetime
import json
import numpy as np
from string import replace


try:
    # Import pymatgen Class 
    import pymatgen_core.core.units as units 
except:
    raise ImportError("pymatgen import error for units object")
    

import streamm.structures 
#import streamm.calculations.resource as resource 
from resource import Resource 
from resource import CalculationRes
    
import logging
logger = logging.getLogger(__name__)


def conv_float(fval_i):
    '''
    Convert a value to a float 
    '''
    sval_i = replace(str(fval_i),"D","e")
    
    try:
        val_i = float(sval_i)
    except:
        val_i = 'Nan'
        
    return val_i

class electrontransfer(object):
    """
    Calculation of electron transfer between groups of particles
    """
    def __init__(self,  verbose=False):
        """
        Constructor for  class. 
        """
        self.producten = 0.0 
        self.reactanten = 0.0 
        self.productMO = ""
        self.reactantMO = ""
        self.V = 0.0
        self.S = 0.0
        self.cputime = ""

    def __str__(self):
        log_line = ""
        log_line += "\n producten {}".format(self.producten)
        log_line += "\n reactanten {}".format(self.reactanten)
        log_line += "\n productMO {}".format(self.productMO)
        log_line += "\n reactantMO {}".format(self.reactantMO)
        log_line += "\n S {} H ".format(self.S)
        log_line += "\n V {} H ".format(self.V)
        log_line += "\n cputime {}".format(self.cputime)

        return log_line
                            
class NWChem(CalculationRes):
    """
    Dervied class implementing input/output methods Gaussian
    """
    def __init__(self, tag , verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """
        # Set units for LAMMPS 
        unit_conf = units.unit_conf
        unit_conf['engry'] = 'Ha'
        unit_conf['length'] = 'ang'
        unit_conf['charge'] = 'e'
        unit_conf['time'] = 'ns'
        
        # Base class constructor is called
        CalculationRes.__init__(self, tag,unit_conf=unit_conf)

        self.meta['software'] = 'nwchem'
        # String found in log file when simulation finishes
        self.properties['finish_str'] = 'Total times  cpu:'
        #
        self.et_list = []
        
    def __del__(self):
        """
        Destructor, clears object memory
        """
        # Call base class destructor
        CalculationRes.__del__(self)
        

    def proc_log(self,log_file):
        """
        Read NWChem simulation log files  
        """
        # verbose = False 
        debug = False

        self.calctype = 'et'
        self.properties['calctype']  = 'SP'

        # self.converged = False
        self.et_list = []
        
        try:
            with open(log_file,"r") as F:
                log_lines = F.readlines()
                F.close()
                # self.converged = True
                self.properties['stinfo'] = []
                self.properties['N_alpha_occ'] = 0
                self.properties['N_beta_occ'] = 0
                self.properties['alpha_energies'] = []
                self.properties['beta_energies'] = []
                self.properties['calctype']  = 'SP'
                self.properties['nstates']  = 0
                self.properties['energy']  = 0.0
                #
                # Initialize read functions as off 
                # 
                read_et = False 
                read_geom = False
                read_molorben_alpha = False
                read_molorben_beta = False
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
                                et_ij.S = conv_float(col[4]) # float( replace(str(col[4]),"D","e") )
                                
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
                                self.et_list.append(copy.deepcopy(et_ij))
                                
                    if( len(col) >= 3 and read_et == False ):
                        if(  col[0]  == "Electron" and  col[1]  == "Transfer"  and  col[2]  == "Calculation" ):
                            read_et = True
                            et_ij = electrontransfer()
                            logger.debug(" Electron Transfer Calculation found on line {} ".format(line_cnt))


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
                        
                    if( read_molorben_alpha and len(col) >= 4 ):
                        if( col[0] == 'Vector' ):
                            col_s = col[2].split('=')
                            occ = conv_float(col_s[1]) #.replace("D","e"))
                            if( occ > 0 ):
                                self.properties['N_alpha_occ'] += 1 
                            if( len(col) == 4  ):
                                en_return =  col[3].split('=')
                                en = conv_float(en_return[1]) #.replace("D","e")) 
                            elif( len(col) == 5  ):
                                en = conv_float(col[4]) #.replace("D","e"))
                            else:
                                error_msg = " Bad MO read {} {}".format(col,len(col))
                                raise IOError(error_msg)
                                
                            self.properties['alpha_energies'].append(en)
                            
                            #print occ,en
                            

                    if( read_molorben_beta and len(col) >= 4 ):
                        if( col[0] == 'Vector' ):
                            col_s = col[2].split('=')
                            occ = conv_float(col_s[1]) #float(col_s[1].replace("D","e"))
                            if( occ > 0 ):
                                self.properties['N_beta_occ'] += 1 
                            if( len(col) == 4  ):
                                en_return =  col[3].split('=')
                                en = conv_float(en_return[1]) #.replace("D","e")) 
                            elif( len(col) == 5  ):
                                en = conv_float(col[4]) #.replace("D","e"))
                            else:
                                error_msg = " Bad MO read {} {}".format(col,len(col))
                                raise IOError(error_msg)
                                
                            self.properties['beta_energies'].append(en)
                        
                            
                            
                    if( len(col)  >= 5  ):
                        if(  col[1]  == "Final" and col[3] == 'Molecular' and col[4] == 'Orbital' and col[5] == 'Analysis'  ):
                            print "Checking MO Energies "
                            if( col[2]  == "Alpha" ):
                               read_molorben_alpha = True # True 
                               cnt_molorben_alpha = 0
                            if( col[2]  == "Beta" ):
                               read_molorben_alpha = False # True 
                               read_molorben_beta = True # True 
                               cnt_molorben_beta = 0
                               
                    if( len(col)  >= 5  ):
                        #if( 'Total SCF energy' in line ):
                        if(  col[0]  == "Total" and col[1] == 'SCF' and col[2] == 'energy' ):
                            self.properties['energy'] =  conv_float(col[4]) # float( replace(str(col[4]),"D","e") )
                        


                logger.debug(" ---- {} Electron Transfer Calculations found ".format(len(self.et_list)))
                    

        except IOError:
            logger.warning("Error output file %s not found "%log_file)
            #self.converged = False
                
            
        return


    def analysis(self,output_key='log'):
        """
        Read in results from NWChem 
        """
        # Find output_key file 
        try:
            output_file = self.files['output'][output_key]
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])                              
                bash_command = "scp  %s:%s%s  ./ "%(ssh_id,self.dir['scratch'],output_file)
                os.system(bash_command)                
            self.proc_log(output_file)            
        except KeyError:
            print "Calculation %s No output_file file  with key %s found"%(self.tag,output_key)
        
    def load_json(self):
        '''
        Load json file for reference 
        '''        
        json_file = "%s_%s.json"%(self.prefix,self.tag)
        try:
            with open(json_file) as f:            
                json_data = json.load(f)
                f.close()
                
                self.meta = json_data['meta']
                self.units = json_data['units']
                self.files = json_data['files']
                self.data = json_data['data'] 
                self.properties = json_data['properties'] 
                self.dir = json_data['dir']

                # Load resource 
                try:
                    res_tag = json_data['meta']['resource']
                except:
                    res_tag = ''
                    print "No resource found "
                if( len(res_tag) > 0 ):
                    print "Resource tag found %s "%(res_tag)
                    resource_i = Resource(str(res_tag))
                    resource_i.load_json()
                    self.resource = resource_i
                                    
                # Load references 
                try:
                    ref_tags = json_data['references']
                    for rekey,ref_tag in ref_tags:
                        ref_i = Calculation(ref_tag)
                        ref_i.load_json()
                        self.add_refcalc(ref_i)
                        print " Need to set reference calculation type "
                except:
                    print "No references found "
                    


        except IOError:
            logger.warning(" File not found %s in %s "%(json_file,os.getcwd()))
            print " File not found %s in %s "%(json_file,os.getcwd())

            
