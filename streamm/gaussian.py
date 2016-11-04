"""
Class data structures for Gaussian data
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"

import copy, sys, os, shutil, math
import time, datetime
import json
import numpy as np
from string import replace


import structure, parameters, units, periodictable, resource, buildingblock
from calculation import CalculationRes

import logging
logger = logging.getLogger(__name__)
          
class Gaussian(CalculationRes):
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

    def proc_log(self,log_file):
        '''
        Read data from gaussian log file
        '''
        def ns_get(fline):
            """
            Parses a line containing 'td=nstates=ns/' embedded in other
                text to find and return the number of excited states from a
                Gaussian TDDFT calculation
            """
            # line has format * td=nstates=ns/* where ns = number of states
            l = fline[fline.index('nstates='):] # Remove everything before 'nstates='
            ll = l[l.index('='):] # Remove everything before =ns*
            lll = ll[:ll.index('/')] # Remove / and everything after
            ns = int(lll.strip('=')) # Remove the = and assign to integer ns
            
            return ns


        def spec_get(fline):
            """
            Parses a line containing "Excited State" to return the energy gap
                (eV), the absorption wavelength (nm) and the oscillator strength
            """
            st_list = fline.split()
            s = st_list[8] # The 9th element is f=XXX
                #st_list[8] = s[2:] # Remove the f= part from this element
            fosc= float(s[2:])
                #(st_list[4],st_list[6],st_list[8])
            egap = float(st_list[4])
            lam = float(st_list[6])
            
            return (egap,lam,fosc)

        def check_get(fline):
            """
            Parses a line containing "%chk=" to get the root part of the checkpoint
                file name.
                For example, "%chk=file.chk" returns the string "file"
            """
            l = fline[fline.index('=')+1:]
            ll = l[:l.index('.chk')]
            
            return ll
            
        f = open(log_file,'r')
        log_lines = f.readlines()
        f.close()
        # Set properties to read in 
        self.properties['stinfo'] = []
        self.properties['N_alpha_occ'] = 0
        self.properties['N_beta_occ'] = 0
        self.properties['alpha_energies'] = []
        self.properties['beta_energies'] = []
        self.properties['calctype']  = 'SP'
        self.properties['nstates']  = 0
        self.properties['energy']  = 0.0
        
        # self.properties['method']  = ''
        # self.properties['basis']  = ''

        for line in log_lines:
            llow = line.lower()
            # Get the number of excited states calculated 
            if 'td=nstates=' in llow:
                self.properties['calctype'] = 'td'
                self.properties['nstates'] = ns_get(llow)
            # Get the excited state energies  
            if 'Excited State' in line:
                ex_info = (egap,lam,fosc) = spec_get(line)
                self.properties['stinfo'].append(ex_info)
            # Get the dipole 
            if ('Tot=' in line and 'X=' in line ): # Get the dipole vector for the ground state and length (in Debye)
                dip_list = line.split()
                self.strucC.properties['dipole'] = [ dip_list[1], dip_list[3], dip_list[5] ]


            # Get the occupied and unoccupied levels
            if 'occ.' in line :
                lsplit = line.split()
                if 'Alpha' in line:
                    for v in lsplit[4:]:
                        self.properties['N_alpha_occ'] += 1 
                        self.properties['alpha_energies'].append(float(v))
                if 'Beta' in line:
                    betas = True
                    for v in lsplit[4:]:
                        self.properties['N_beta_occ'] += 1 
                        self.properties['beta_energies'].append(float(v))
            if 'virt.' in line :
                lsplit = line.split()
                if 'Alpha' in line:
                    for v in lsplit[4:]:
                        self.properties['alpha_energies'].append(float(v))
                if 'Beta' in line:
                    betas = True
                    for v in lsplit[4:]:
                        self.properties['beta_energies'].append(float(v))

            if ('SCF Done' in line ):
                lsplit = line.split()
                self.properties['energy'] = float(lsplit[4])
                l = lsplit[2].split("(")
                self.properties['method'] = l[1].strip(")")

            if ('Standard basis' in line ):
                lsplit = line.split()
                self.properties['basis'] = lsplit[2]
                     
    def proc_fchk(self,fchk_file):
        """
        Read in structure information from gaussian fchk file

        Args:
            fchk_file (str) gaussian fchk file

        """
        n_dim = self.strucC.lat.n_dim 
        
        with open(fchk_file) as F:
            Lines = F.readlines()
            F.close()
            #
            read_r = False
            read_eln = False
            read_esp = False
            read_mulliken = False
            read_alphaE = False
            #
            line_cnt = 0 
            for line in Lines :
                col = line.split()
                line_cnt += 1
                if( read_r ):
                    # Read positions 
                    if( col[0] == "Force" and col[1] == "Field" ):
                        read_r = False
                        for p_indx in range(n_particles_i):
                            vec_r_i =  R_all[p_indx*n_dim:p_indx*n_dim+n_dim]
                            self.strucC.positions[p_indx] =  np.array(vec_r_i)
                    else:
                        for r_i in  map(float,col) :
                            R_all.append(r_i) 

                if( read_eln ):
                    if ( eln_p_cnt == n_particles_i ):
                        read_eln = False
                    else:                    
                        for eln_i in  map(int,col):
                            part_i = self.strucC.particles[eln_p_cnt]
                            part_i = periodictable.element_number(eln_i)
                            eln_p_cnt += 1

                if( read_esp ):
                    if ( esp_p_cnt == n_particles_i ):
                        read_esp = False
                    else:
                        for q_i in  map(float,col):
                            part_i =  self.strucC.particles[esp_p_cnt]
                            part_i.properties['esp'] =  q_i
                            esp_p_cnt += 1

                if( read_mulliken ):
                    if ( mulliken_p_cnt == n_particles_i ):
                        read_mulliken = False
                    else:
                        for q_i in  map(float,col):
                            part_i =  self.strucC.particles[mulliken_p_cnt]
                            part_i.properties['mulliken'] =  q_i
                            mulliken_p_cnt += 1

                if( read_alphaE ):
                    if ( len( self.properties['alpha_energies'] ) == self.properties['N_alphae'] ):
                        read_alphaE = False
                    else:
                        for val_i in  map(float,col):
                            self.properties['alpha_energies'].append( val_i) 

                #
                # Read in and initialize structure and 
                #
                if( len(col) > 2 ):
                    if( col[0] == "Number" and col[2] == "electrons" ):
                        self.properties['N_electrons'] = int( col[4] ) 
                if( len(col) > 2 ):
                    if( col[0] == "Number" and col[2] == "alpha" ):
                        self.properties['N_alphae'] = int( col[5] ) 
                if( len(col) > 2 ):
                    if( col[0] == "Number" and col[2] == "beta" ):
                        self.properties['N_betae'] = int( col[5] ) 
                if( len(col) > 2 ):
                    if( col[0] == "Number" and col[2] == "basis" ):
                        self.properties['N_basisfunctions'] = int( col[5] )
                if( len(col) > 2 ):
                    if( col[0] == "Total" and col[1] == "Energy" ):
                        self.properties['energy'] = float( col[3] ) 

                if( len(col) == 5 ):
                    if( col[0] == "Number" and col[1] == "of"  and col[2] == "atoms" ):
                        n_particles_i = int(col[4])
                        #if( pt_update ):
                        #    if( n_particles_i != len(strucC.ptclC) ):
                        #        print " json file contains %d atoms and fchk file contains %d "%(len(strucC.ptclC),NA)
                        #        sys.exit("inconsistent files ")
                        #else:
                        # Create particles to be updated
                        if( self.strucC.n_particles == 0 ):
                            print "Particle container is empty will add new particles"
                            for pkey_i in range(n_particles_i):
                                pos_i = np.array([0.0,0.0,0.0])
                                particle_i = structure.Atom()
                                self.strucC.add_partpos(particle_i,pos_i)

                if( line_cnt == 2 ):
                    self.properties['calctype'] = col[0].strip()
                    self.properties['method'] = col[1].strip()
                    self.properties['basis'] = col[2].strip()                                
                #
                # Turn on reading based on key words 
                #           
                if( len(col) == 6 ):
                    if( col[0] == "Current" and col[1] == "cartesian"  and col[2] == "coordinates" ):
                        read_r = True
                        R_all = []
                if( len(col) > 2  ):
                    if( col[0] == "Atomic" and col[1] == "numbers"   ):
                        read_eln = True
                        eln_p_cnt = 0

                if( len(col) > 2  ):
                    if( col[0] == "ESP" and col[1] == "Charges"   ):
                        read_esp = True
                        esp_p_cnt = 0

                if( len(col) > 2  ):
                    if( col[0] == "Mulliken" and col[1] == "Charges"   ):
                        read_mulliken = True
                        mulliken_p_cnt = 0

                if( len(col) > 2  ):
                    if( col[0] == "Alpha" and col[1] == "Orbital" and col[2] == "Energies"   ):
                        read_alphaE = True
                        self.properties['N_alphae'] = int(col[5])
                        self.properties['alpha_energies'] = []

            return
                        
    def analysis(self,output_key='log'):
        """
        Read in results from gaussian 
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
        # Read fchk if it exists
        fchk_key = 'fchk'
        try:
            output_file = self.files['output'][fchk_key]
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])                              
                bash_command = "scp  %s:%s%s  ./ "%(ssh_id,self.dir['scratch'],output_file)
                os.system(bash_command)
            self.proc_fchk(output_file)            
        except KeyError:
            print "Calculation %s No output_file file  with key %s found"%(self.tag,fchk_key)
        
        #self.write_dat(dat_file)
        #self.add_file(dat_file,'data')


