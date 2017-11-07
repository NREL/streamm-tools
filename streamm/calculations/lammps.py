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
Class data structures for LAMMPS data
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
    

from resource import Resource 
from calculation import Calculation
from calculation import MDrun


import logging
logger = logging.getLogger(__name__)

class LAMMPS(Calculation):
    """
    Derived class implementing input/output methods for LAMMPS
    """


    def __init__(self, tag ):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """
        # Set units for LAMMPS 
        unit_conf = units.unit_conf
        unit_conf['engry'] = 'kCalmol'
        unit_conf['length'] = 'ang'
        unit_conf['charge'] = 'e'
        unit_conf['time'] = 'ns'
        
        # Base class constructor is called
        Calculation.__init__(self, tag,unit_conf=unit_conf)

        self.meta['software'] = 'lammps'
        self.properties['finish_str'] = 'Loop time of'
        #
    def __del__(self):
        """
        Destructor, clears object memory
        """
        # Call base class destructor
        Calculation.__del__(self)

    def read_data(self, data_file,
        btype = "harmonic",
        atype = "harmonic",
        dtype = "multiharmonic",
        imptype = "multiharmonic",
        debug = False):

        """
        Read Lammps data file

        Args:
            * data_file  (str) data file
            
        """

        F = open(data_file , 'r' )
        lines = F.readlines()
        F.close()
        logger.debug(" Reading {} ",data_file)
        #
        # Read in data header with number of parameters 
        #
        matrix = self.strucC.lat.matrix 
        for line in lines:
            col = line.split()

            if ( len(col) >=2 ):
                # Read in number of each topolgical component  
                if( col[1] == "atoms" ):
                    self.strucC.n_particles = int( col[0] )
                elif( col[1] == "bonds" ):
                    self.strucC.n_bonds = int( col[0] )
                elif( col[1] == "angles" ):
                    self.strucC.n_angles = int( col[0] )
                elif( col[1] == "dihedrals" ):
                    self.strucC.n_dihedrals = int( col[0] )
                elif( col[1] == "impropers" ):
                    self.strucC.n_impropers = int( col[0] )
                    
            if ( len(col) >= 3 ):
                # Read in number of each parameter type 
                if( col[1] == "atom" and   col[2] == "types" ):
                    self.paramC.n_ljtypes = int( col[0] )
                elif( col[1] == "bond" and   col[2] == "types"  ):
                    self.paramC.n_bondtypes = int( col[0] )
                elif( col[1] == "angle" and   col[2] == "types"  ):
                   self.paramC.n_angletypes = int( col[0] )
                elif( col[1] == "dihedral" and   col[2] == "types"  ):
                    self.paramC.n_dihtypes = int( col[0] )
                elif( col[1] == "improper" and   col[2] == "types"  ):
                    self.paramC.n_imptypes = int( col[0] )
            
            # Read in box size    
            if ( len(col) >= 4 ):
                if( col[2]  == "xlo"  and col[3]  == "xhi"  ):
                    matrix[0][0] = float( col[1] ) - float( col[0] )
                    matrix[0][1] = 0.0  
                    matrix[0][2] = 0.0  
                if( col[2]  == "ylo"  and col[3]  == "yhi"  ):
                    matrix[1][0] = 0.0  
                    matrix[1][1] = float( col[1] ) - float( col[0] ) 
                    matrix[1][2] = 0.0  
                if( col[2]  == "zlo"  and col[3]  == "zhi"  ):
                    matrix[2][0] = 0.0  
                    matrix[2][1] = 0.0  
                    matrix[2][2] = float( col[1] ) - float( col[0] )
                    # Set lattice
                    self.strucC.lat.set_matrix(matrix)
                   
        #      
        # Initialize blank structure and parameters in case entries are not in numerical order 
        #
        # Structure 
        pos_i = np.zeros(self.strucC.lat.n_dim)
        self.strucC.positions = [ pos_i for pkey_i in range(self.strucC.n_particles) ]
        #
        # Intialize read in boolean to off                    
        #
        read_Masses = False
        read_Pair = False
        read_Bond_coeff = False
        read_Angle_coeff = False
        read_Dihedral_coeff = False
        read_Improper_coeff = False

        read_Atoms = False
        read_Bonds = False
        read_Angles = False
        read_Dihedrals = False
        read_Impropers = False
  
        #
        # Read in data parameters 
        #
        for line in lines:
            col = line.split()
            
            if( read_Masses and  len(col) >= 2 ):
                cnt_Masses += 1
                ljkey_i = int( col[0]) -1
                
                mass_i = float(col[1])
                el_i = periodictable.element_mass(mass_i)
                fftype1 = str(el_i["symbol"]) + str(cnt_Masses)
                ljtype_i = parameters.LJtype(fftype1)
                ljtype_i.lammps_index = int(col[0])
                ljtype_i.mass = mass_i
                ljtype_i.atomic_symbol = el_i["symbol"]

                #print ">read_data ljtype_i",ljkey_i,ljtype_i.mass,ljtype_i.lammps_index,cnt_Masses ,self.paramC.n_ljtypes 
                
                self.paramC.ljtypes[ljkey_i] = copy.deepcopy(ljtype_i)
                #self.paramC.add_LJtype(ljtype_i,deepcopy = True )

                # Turn of mass read
                if(cnt_Masses ==  self.paramC.n_ljtypes ):
                    read_Masses = False 

            if( read_Pair and  len(col) >= 3 ):
                cnt_Pair += 1
                
                ljkey_i = int(col[0]) - 1
                
                #print ">read_data read_Pair",ljkey_i #,ljtype_i.epsilon,ljtype_i.sigma
                ljtype_i = self.paramC.ljtypes[ljkey_i]
                ljtype_i.epsilon = float(col[1])
                ljtype_i.sigma = float(col[2])
                
                #print ">read_data read_Pair",ljkey_i,ljtype_i.epsilon,ljtype_i.sigma,cnt_Pair ,  self.paramC.n_ljtypes
                # ljtype_i.setparam(epsilon,sigma)
                # Turn pair parameter read off 
                if( cnt_Pair ==  self.paramC.n_ljtypes ):
                    read_Pair = False


            if( read_Bond_coeff and  len(col) >= 3 ):
                cnt_Bond_coeff += 1
                bkey_i = int(col[0]) - 1
                
                bondtype_i = parameters.Bondtype(type=btype)
                bondtype_i.kb = float(col[1]) 
                bondtype_i.r0 = float(col[2])
                bondtype_i.lammps_index = int( col[0])
                self.paramC.bondtypes[bkey_i] = copy.deepcopy(bondtype_i)
                
                # btyp_i.setharmonic(r0,kb)
                if( cnt_Bond_coeff >=  self.paramC.n_bondtypes ):
                    read_Bond_coeff = False


            if( read_Angle_coeff and  len(col) >= 3 ):
                cnt_Angle_coeff += 1
                akey_i = int(col[0]) - 1
                
                angletype_i = parameters.Angletype(type=atype)
                angletype_i.kb = float(col[1]) 
                angletype_i.theta0 = float(col[2])
                angletype_i.lammps_index = int( col[0])
                self.paramC.angletypes[akey_i] = copy.deepcopy(angletype_i)
                                
                if( cnt_Angle_coeff >= self.paramC.n_angletypes ):
                    read_Angle_coeff = False


            if( read_Dihedral_coeff and  len(col) >= 3 ):
                cnt_Dihedral_coeff += 1
                dkey_i = int(col[0]) - 1

                dihtype_i = parameters.Dihtype(type=dtype)
                if( dtype == "opls" ):
                    dihtype_i.setopls(float(col[1]),float(col[2]),float(col[3]),float(col[4]))
                    dihtype_i.g_indx = int(3)
                elif( dtype == "harmonic" ):
                    dihtype_i.setharmonic(float(col[2]),float(col[1]),float(col[3]))
                    dihtype_i.g_indx = int(1)
                elif( dtype == "multiharmonic" ):
                    dihtype_i.setharmonic(float(col[2]),float(col[1]),float(col[3]))
                    dihtype_i.g_indx = int(1)
                dihtype_i.lammps_index = int( col[0])
                self.paramC.dihtypes[dkey_i] = copy.deepcopy(dihtype_i)
                
                if( cnt_Dihedral_coeff >= self.paramC.n_dihtypes ):
                    logger.debug( " %d dihedral types read in "%(cnt_Dihedral_coeff))
                    read_Dihedral_coeff = False


            if( read_Improper_coeff and  len(col) >= 3 ):
                cnt_Improper_coeff += 1
                ikey_i = int(col[0]) - 1

                imptype_i = parameters.Imptype(type=imptype)
                if( imptype == "improper" ):
                    imptype_i.setimp(float(col[2]),float(col[1]))
                    imptype_i.g_indx = int(2)
                imptype_i.lammps_index = int( col[0])
                self.paramC.imptypes[ikey_i] = copy.deepcopy(imptype_i)

                if( cnt_Improper_coeff >= self.paramC.n_imptypes ):
                    read_Improper_coeff = False



            if( read_Atoms and len(col) >= 7 ):

                # print ">read_Atoms col",col
                
                cnt_Atoms += 1
                
                pkey_i = int( col[0]) - 1
                ljkey_i = int( col[2]) - 1
                
                ljtype_i = self.paramC.ljtypes[ljkey_i]
                if( ljtype_i.lammps_index != ljkey_i +1 ):
                    lerror_msg = "Read in error for atom %s due to bad Pair type %d "%(cnt_Atoms,ljtype_i.lammps_index)
                    raise ValueError(error_msg)
                # set particle and  properties
                try:
                    particle_i = self.strucC.particles[pkey_i]
                    try:
                        ljtype_i.fftype1 = particle_i.properties['fftype']
                    except:
                        print " Particle %d has no ffytpe"%(pkey_i)
                except:
                    particle_i = buildingblock.BBatom(ljtype_i.atomic_symbol)
                particle_i.properties["mol"] = int(col[1])      
                particle_i.properties["charge"] = float(col[3])
                particle_i.properties["mass"] = ljtype_i.mass
                particle_i.properties["lmpindx"] = ljtype_i.lammps_index

                self.strucC.particles[pkey_i] = copy.deepcopy(particle_i)
                # set position 
                pos_i =  [ float(col[4]),float(col[5]),float(col[6])] 
                self.strucC.positions[pkey_i] = pos_i

                if( cnt_Atoms >=  self.strucC.n_particles ):
                    read_Atoms = False

            if(read_Bonds and len(col) >= 4 ):
                cnt_Bonds += 1
                
                bkey_i = int( col[0]) - 1
                pkey1 = int(col[2]) - 1
                pkey2 = int(col[3]) - 1
                
                bond_i = structure.Bond(pkey1,pkey2)
                bond_i.lammps_index = int(col[1])

                typekey =  bond_i.lammps_index  - 1
                try:
                    bondtype_i = self.paramC.bondtypes[typekey]
                    try:
                        bondtype_i.fftype1 = self.strucC.particles[pkey1].properties['fftype']
                        bondtype_i.fftype2 = self.strucC.particles[pkey2].properties['fftype']
                    except:
                        print " Particles %d %d has no ffytpe"%(pkey1,pkey2)
                except:
                    print "Bond type  %d is not set "%(typekey)
                    
                
                self.strucC.bonds[bkey_i] = copy.deepcopy(bond_i)

                if( cnt_Bonds >=  self.strucC.n_bonds ):
                    read_Bonds = False

            if(read_Angles and len(col) >= 5 ):
                cnt_Angles += 1
                
                akey_i = int( col[0]) - 1
                pkey1 = int(col[2]) - 1
                pkey2 = int(col[3]) - 1
                pkey3 = int(col[4]) - 1
                angle_i = structure.Angle(pkey1,pkey2,pkey3)
                angle_i.lammps_index = int(col[1])

                typekey =  angle_i.lammps_index  - 1
                try:
                    angletype_i = self.paramC.angletypes[typekey]
                    try:
                        angletype_i.fftype1 = self.strucC.particles[pkey1].properties['fftype']
                        angletype_i.fftype2 = self.strucC.particles[pkey2].properties['fftype']
                        angletype_i.fftype3 = self.strucC.particles[pkey3].properties['fftype']
                    except:
                        print " Particles %d %d has no ffytpe"%(pkey1,pkey2,pkey3)
                except:
                    print "angletype type  %d is not set "%(typekey)
                                    
                self.strucC.angles[akey_i] = copy.deepcopy(angle_i)

                if( cnt_Angles >=  self.strucC.n_angles ):
                    read_Angles = False


            if(read_Dihedrals and len(col) >= 6 ):
                cnt_Dihedrals += 1

                dkey_i = int( col[0]) - 1
                pkey1 = int(col[2]) - 1
                pkey2 = int(col[3]) - 1
                pkey3 = int(col[4]) - 1
                pkey4 = int(col[5]) - 1
                dihedral_i = structure.Dihedral(pkey1,pkey2,pkey3,pkey4)
                dihedral_i.lammps_index = int(col[1])
                # Set parameter types 
                typekey =  dihedral_i.lammps_index  - 1
                try:
                    dihtype_i = self.paramC.dihtypes[typekey]
                    try:
                        dihtype_i.fftype1 = self.strucC.particles[pkey1].properties['fftype']
                        dihtype_i.fftype2 = self.strucC.particles[pkey2].properties['fftype']
                        dihtype_i.fftype3 = self.strucC.particles[pkey3].properties['fftype']
                        dihtype_i.fftype4 = self.strucC.particles[pkey4].properties['fftype']
                    except:
                        logger.warning(" Particles %d %d has no ffytpe"%(pkey1,pkey2,pkey3,pkey4))
                except:
                    logger.warning("dihtype type  %d is not set "%(typekey))

                    
                self.strucC.dihedrals[dkey_i] = copy.deepcopy(dihedral_i)

                if( cnt_Dihedrals >=  self.strucC.n_dihedrals ):
                    read_Dihedrals = False


            if(read_Impropers and len(col) >= 2 ):
                cnt_Impropers += 1


                ikey_i = int( col[0]) - 1
                pkey1 = int(col[2]) - 1
                pkey2 = int(col[3]) - 1
                pkey3 = int(col[4]) - 1
                pkey4 = int(col[5]) - 1
                improper_i = structure.Improper(pkey1,pkey2,pkey3,pkey4)
                improper_i.lammps_index = int(col[1])
                # Set parameter types 
                typekey =  improper_i.lammps_index  - 1
                try:
                    improper_i = self.paramC.imptypes[typekey]
                    try:
                        improper_i.fftype1 = self.strucC.particles[pkey1].properties['fftype']
                        improper_i.fftype2 = self.strucC.particles[pkey2].properties['fftype']
                        improper_i.fftype3 = self.strucC.particles[pkey3].properties['fftype']
                        improper_i.fftype4 = self.strucC.particles[pkey4].properties['fftype']
                    except:
                        print " Particles %d %d has no ffytpe"%(pkey1,pkey2,pkey3,pkey4)
                except:
                    print "improper type  %d is not set "%(typekey)
                    
                self.strucC.impropers[ikey_i] = copy.deepcopy(improper_i)
                
                if( cnt_Impropers >=  self.strucC.n_impropers ):
                    read_Impropers = False

            if ( len(col) >= 1  ):
                if( col[0] == "Masses" ):
                    read_Masses = True
                    cnt_Masses = 0
                    # 
                    logger.debug("Reading Masses  ")
                    # 
                if( col[0] == "Atoms" ):
                    read_Atoms = True
                    cnt_Atoms = 0
                    # 
                    logger.debug("Reading Atoms   ")
                    # 
                    # for ljtkey_i, ljtype_i  in self.paramC.ljtypes.iteritems():
                    #    print ljtkey_i, ljtype_i.lammps_index , ljtype_i.mass ,  ljtype_i.fftype1,ljtype_i.epsilon,  ljtype_i.sigma 
                    # 
                if( col[0] == "Bonds" ):
                    # 
                    logger.debug("Reading Bonds ")
                    # 
                    read_Bonds = True
                    cnt_Bonds = 0 
                if( col[0] == "Angles" ):
                    # 
                    logger.debug("Reading Angles  ")
                    # 
                    read_Angles = True
                    cnt_Angles = 0 
                if( col[0] == "Dihedrals" ):
                    # 
                    logger.debug("Reading Dihedrals  ")
                    # 
                    read_Dihedrals = True
                    cnt_Dihedrals = 0 
                if( col[0] == "Impropers" ):
                    # 
                    logger.debug("Reading  Impropers ")
                    # 
                    read_Impropers = True
                    cnt_Impropers = 0 
            if ( len(col) >= 2 ):
                if( col[0] == "Pair" and col[1] == "Coeffs" ):
                    # 
                    logger.debug("Reading Pairs  ")
                    # 
                    read_Pair = True
                    cnt_Pair = 0 
                if( col[0] == "Bond" and col[1] == "Coeffs" ):
                    # 
                    logger.debug("Reading Bond_coeff  ")
                    # 
                    read_Bond_coeff = True
                    cnt_Bond_coeff = 0 
                if( col[0] == "Angle" and col[1] == "Coeffs" ):
                    # 
                    logger.debug("Reading Angle_coeff  ")
                    # 
                    read_Angle_coeff  = True
                    cnt_Angle_coeff = 0 
                if( col[0] == "Dihedral" and col[1] == "Coeffs" ):
                    # 
                    logger.debug("Reading Dihedral_coeff   ")
                    # 
                    read_Dihedral_coeff  = True
                    cnt_Dihedral_coeff = 0 
                if( col[0] == "Improper" and col[1] == "Coeffs" ):
                    # 
                    logger.debug("Reading Improper_coeff   ")
                    # 
                    read_Improper_coeff  = True
                    cnt_Improper_coeff = 0 
        #      
        return  

    def read_data_pos(self, data_file,
        btype = "harmonic",
        atype = "harmonic",
        dtype = "multiharmonic",
        imptype = "improper",
        debug = False):

        """
        Read positions only from Lammps data file

        Args:
            * data_file  (str) data file
            
        """

        F = open(data_file , 'r' )
        lines = F.readlines()
        F.close()
        logger.debug(" Reading %s "%(data_file))
        #
        # Read in data header with number of parameters 
        #
        matrix = self.strucC.lat.matrix 
        for line in lines:
            col = line.split()

            if ( len(col) >=2 ):
                # Read in number of each topolgical component  
                if( col[1] == "atoms" and self.strucC.n_particles != int( col[0] )):
                    print "LAMMPS data file %s has %d particles while strucC object has %d particles "%(data_file, int(col[0]), self.strucC.n_particles )
                    return 

            # Read in box size    
            if ( len(col) >= 4 ):
                if( col[2]  == "xlo"  and col[3]  == "xhi"  ):
                    matrix[0][0] = float( col[1] ) - float( col[0] )
                    matrix[0][1] = 0.0  
                    matrix[0][2] = 0.0  
                if( col[2]  == "ylo"  and col[3]  == "yhi"  ):
                    matrix[1][0] = 0.0  
                    matrix[1][1] = float( col[1] ) - float( col[0] ) 
                    matrix[1][2] = 0.0  
                if( col[2]  == "zlo"  and col[3]  == "zhi"  ):
                    matrix[2][0] = 0.0  
                    matrix[2][1] = 0.0  
                    matrix[2][2] = float( col[1] ) - float( col[0] )
                    # Set lattice
                    self.strucC.lat.matrix = matrix
             
        #
        # Intialize read in boolean to off                    
        #
        read_Atoms = False 
        #
        # Read in data parameters 
        #
        for line in lines:
            col = line.split()
            

            if( read_Atoms and len(col) >= 7 ):

                # print ">read_Atoms col",col
                
                cnt_Atoms += 1
                
                pkey_i = int( col[0]) - 1
                # set position 
                pos_i =  [ float(col[4]),float(col[5]),float(col[6])]
                
                
                self.strucC.positions[pkey_i] = pos_i

                if( cnt_Atoms >=  self.strucC.n_particles ):
                    read_Atoms = False

            if ( len(col) >= 1  ):
                if( col[0] == "Atoms" ):
                    read_Atoms = True
                    cnt_Atoms = 0
                    # 
                    logger.debug("Reading Atoms   ")
        #      
        return  


    
    def write_data(self,data_file=''):
        """
        Write data file
        """
        if( len(data_file) == 0 ):
            data_file = "%s.data"%(self.tag)
        
        self.add_file('input','data_file',data_file)
        self.properties['data_file'] = data_file

        F = open( data_file, 'w' )
        F.write('  Lammps data file \n')
        F.write('\n')
        F.write( "%10d  atoms \n" % self.strucC.n_particles )
        F.write( "%10d  bonds \n" %  self.strucC.n_bonds )
        F.write( "%10d  angles \n" % self.strucC.n_angles )
        F.write( "%10d  dihedrals \n" %  self.strucC.n_dihedrals )
        F.write( "%10d  impropers \n" % self.strucC.n_impropers  )
        F.write('\n')
        F.write( "%10d  atom types \n" % self.paramC.n_particletypes  )
        F.write( "%10d  bond types \n" % self.paramC.n_bondtypes )
        F.write( "%10d  angle types \n" % self.paramC.n_angletypes )
        F.write( "%10d  dihedral types \n" % self.paramC.n_dihtypes )
        if( self.paramC.n_imptypes > 0 ):
            F.write( "%10d  improper types \n" % self.paramC.n_imptypes )
        else:
            n=1
            F.write( "%10d  improper types \n" % n )
            
        F.write('\n')
        matrix = self.strucC.lat.matrix
        F.write( "%16.8f %16.8f   xlo xhi \n" %  (matrix[0][0]/-2.0 , matrix[0][0]/2.0) )
        F.write( "%16.8f %16.8f   ylo yhi \n" %  (matrix[1][1]/-2.0 , matrix[1][1]/2.0 ) )
        F.write( "%16.8f %16.8f   zlo zhi \n" %  (matrix[2][2]/-2.0 , matrix[2][2]/2.0) )
        F.write('\n')
        F.write( ' Masses \n')
        F.write('\n')
        # Write LJtypes mass 
        for ptk, pt  in self.paramC.particletypes.iteritems():
            F.write( "%10d %16.8f   # %5s \n" % ( pt.lammps_index, pt.mass , pt.fftype1  ) )
        F.write('\n')
        F.write(' Pair Coeffs \n')
        F.write('\n')
        # Write LJtypes pair Coeffs 
        for ptk, pt  in self.paramC.particletypes.iteritems():
            F.write( "%10d %12.6f %12.6f  \n" % (pt.lammps_index, pt.epsilon,  pt.sigma ) )
        F.write('\n')
        # Write Bond Coeffs
        if( self.paramC.n_bondtypes > 0 ):
            F.write(' Bond Coeffs \n')
            F.write('\n')
            for btkey_i,bondtype_i  in self.paramC.bondtypes.iteritems():
                if( bondtype_i.type == "harmonic"):
                    F.write( "%10d %12.6f %12.6f # %5s %5s  \n" % (bondtype_i.lammps_index,bondtype_i.kb,bondtype_i.r0, bondtype_i.fftype1, bondtype_i.fftype2 ) )
            F.write('\n')

        # Write Angle Coeffs
        if( self.paramC.n_angletypes > 0 ):
            F.write(' Angle Coeffs \n')
            F.write('\n')
            for atkey_i,angletype_i  in self.paramC.angletypes.iteritems():
                if( angletype_i.type == "harmonic"):
                    F.write( "%10d %12.6f %12.6f # %5s %5s  %5s   \n" % (angletype_i.lammps_index,angletype_i.kb,angletype_i.theta0, angletype_i.fftype1,angletype_i.fftype2,angletype_i.fftype3 ) )
            F.write('\n')

        # Write Dihedral Coeffs
        if( self.paramC.n_dihtypes > 0 ):
            F.write(' Dihedral Coeffs \n')
            F.write('\n')
            for dtkey_i, dihtype_i  in self.paramC.dihtypes.iteritems():    
                if( dihtype_i.type == "multiharmonic"):
                    # K = K[1+d cons(n theta)]
            
                    d = dihtype_i.theta_s
                    K = dihtype_i.kb
                    n = dihtype_i.mult
                    p = dihtype_i.paths
                    w = 0.0 # Weight 
                    F.write( "%10d %12.6f %12.6f  %12.6f %12.6f # %d x %5s %5s  %5s  %5s   \n" % (dihtype_i.lammps_index,K,n,d,w, p, dihtype_i.fftype1,dihtype_i.fftype2,dihtype_i.fftype3,dihtype_i.fftype4  ) )
                elif( dihtype_i.type == "rb" or  dihtype_i.type == "opls"  ):
                    # Get opls parameters                    
                    F.write( "%10d  %12.6f  %12.6f  %12.6f  %12.6f # %5s %5s  %5s %5s \n" % (dihtype_i.lammps_index,dihtype_i.k1,dihtype_i.k2,dihtype_i.k3,dihtype_i.k4, dihtype_i.fftype1,dihtype_i.fftype2,dihtype_i.fftype3,dihtype_i.fftype4  ) )
                else:
                    error_msg = " Unknow dihedral type {} ".format(dihtype_i.type )
                    raise ValueError(error_msg)
                    
            F.write('\n')

        # Write Dihedral Coeffs
        if( self.paramC.n_imptypes > 0 ):
            
            F.write(' Improper Coeffs \n')
            F.write('\n')
            if( self.paramC.n_imptypes > 0 ):
                for itkey_i, imptype_i  in self.paramC.imptypes.iteritems():    
                    if( imptype_i.type == "improper"):
                        F.write( "%10d %12.6f %12.6f # %5s %5s  %5s  %5s   \n" % (imptype_i.lammps_index,imptype_i.ke,imptype_i.e0, imptype_i.fftype1,imptype_i.fftype2,imptype_i.fftype3,imptype_i.fftype4  ) )
                    else:
                        error_msg = " Unknow improper type %s "%(imptype_i.type )
                        raise ValueError(error_msg)
                    
        else:
            F.write(' Improper Coeffs \n')
            F.write('\n')
            F.write( "    1 0.0 0.0  \n")
            
        F.write('\n')
        # Write Particles
        if( self.strucC.n_particles > 0 ):
            F.write(' Atoms \n')
            F.write('\n')
            for pkey_i, particle_i  in self.strucC.particles.iteritems():
                fftype_i = particle_i.paramkey
                mol_i = particle_i.mol
                charge_i = particle_i.charge
                lmpindx_i = particle_i.param.lammps_index
                pos_i = self.strucC.positions[pkey_i]       
                F.write( "%9d %9d %8d %12.8f %12.6f %12.6f %12.6f # %5s \n" % (pkey_i+1,mol_i+1,lmpindx_i,charge_i,pos_i[0],pos_i[1],pos_i[2] ,fftype_i)  )
            F.write('\n')
        # Write Bonds
        if( self.strucC.n_bonds > 0 ):
            F.write(' Bonds \n')
            F.write('\n')
            for bkey_i, bond_i  in self.strucC.bonds.iteritems():
                #
                b_i = bond_i.pkey1 + 1 
                b_j = bond_i.pkey2 + 1
                #
                AT_i =  self.strucC.particles[ bond_i.pkey1 ].paramkey
                AT_j =  self.strucC.particles[ bond_i.pkey2 ].paramkey
                #
                F.write(  '%9d %8d %9d %9d # %5s %5s \n' % (bkey_i+1,bond_i.lammps_index,b_i,b_j, AT_i, AT_j ) )
            F.write('\n')

        # Write Angles
        if( self.strucC.n_angles > 0 ):

            F.write(' Angles \n')
            F.write('\n')
            for akey_i, angle_i in self.strucC.angles.iteritems():
                a_k = angle_i.pkey1 + 1 
                a_i = angle_i.pkey2 + 1 
                a_j = angle_i.pkey3 + 1 
                AT_k = self.strucC.particles[ angle_i.pkey1 ].paramkey
                AT_i = self.strucC.particles[ angle_i.pkey2 ].paramkey
                AT_j = self.strucC.particles[ angle_i.pkey3 ].paramkey

                F.write(  '%9d %8d %9d %9d %9d  # %s %s %s \n' % (akey_i+1,angle_i.lammps_index,a_k,a_i,a_j,AT_k,AT_i,AT_j) )
            F.write(  '\n' )

        # Write Dihedrals
        if( self.strucC.n_dihedrals > 0 ):

            F.write(' Dihedrals \n')
            F.write('\n')
            for dkey_i,dih_i in self.strucC.dihedrals.iteritems():
                
                d_k = dih_i.pkey1 + 1 
                d_i = dih_i.pkey2 + 1 
                d_j = dih_i.pkey3 + 1 
                d_l = dih_i.pkey4 + 1 

                AT_k = self.strucC.particles[ dih_i.pkey1 ].paramkey
                AT_i = self.strucC.particles[ dih_i.pkey2 ].paramkey
                AT_j = self.strucC.particles[ dih_i.pkey3 ].paramkey
                AT_l = self.strucC.particles[ dih_i.pkey4 ].paramkey

                F.write(  '%9d %8d %9d %9d %9d %9d # %s %s %s %s \n' % (dkey_i+1,dih_i.lammps_index,d_k,d_i,d_j,d_l,AT_k,AT_i,AT_j,AT_l) )

            F.write( '\n' )

        # Write Impropers
        if( self.strucC.n_impropers > 0 ):
            F.write(' Impropers \n')
            F.write('\n')
            for ikey_i,imp_i in self.strucC.impropers.iteritems():
                d_k = imp_i.pkey1 + 1 
                d_i = imp_i.pkey2 + 1 
                d_j = imp_i.pkey3 + 1 
                d_l = imp_i.pkey4 + 1 
                AT_k = self.strucC.particles[ imp_i.pkey1 ].paramkey
                AT_i = self.strucC.particles[ imp_i.pkey2 ].paramkey
                AT_j = self.strucC.particles[ imp_i.pkey3 ].paramkey
                AT_l = self.strucC.particles[ imp_i.pkey4 ].paramkey
                F.write(  '%9d %8d %9d %9d %9d %9d # %s %s %s %s \n' % (ikey_i+1,imp_i.lammps_index,d_k,d_i,d_j,d_l,AT_k,AT_i,AT_j,AT_l) )

            F.write( '\n' )            


        F.close()
        

    def proc_in(self,in_file,data2cply=True):
        """
        Read in input parameters from lammps in file

        Args:
            * in_file (str) lammps in file

        """
        # Add new properties 
        self.properties['run_cnt'] = 0
        self.properties['run_list'] = []
        dump_cnt = 0 
        ref_struc = copy.deepcopy(self.strucC)
        
        f = open(in_file,'r')
        in_lines = f.readlines()
        f.close()

        read_en = False
        timestep_i = 1.0 # Set to default LAMMPS value 
        for line in in_lines:
            llow = line.lower()
            col = llow.split()
            if( 'timestep' in line):
                
                timestep_i = float(col[1])
                logger.debug("timestep:{}".format(timestep_i))
                
            if( 'run' in line):
                run_i = MDrun()
                run_i.timestep =  timestep_i
                run_i.n_steps = float(col[1])
                
                logger.debug("run found with {} steps".format(run_i.n_steps))
                
                self.properties['run_list'].append(copy.deepcopy(run_i))
                
            if( 'minimize' in line):
                run_i = MDrun()
                run_i.timestep =  timestep_i
                self.properties['run_list'].append(copy.deepcopy(run_i))
                
            if( 'write_data' in line):
                output_file = str(col[1])
                logger.info("LAMMPS .data file found {}".format(output_file))
                self.add_file('output','data_%d'%(len(self.properties['run_list'])),output_file)
                
            if( 'dump' in line and len(col) > 5 ):
                output_file = str(col[5])
                self.add_file('data','dump_%d'%(dump_cnt),output_file)
                dump_cnt += 1
                
            if( 'restart' in line and len(col)  ):
                output_file = "%s.*"%str(col[2])
                self.add_file('data','restart_%d'%(len(self.properties['run_list'])),output_file)

                
        self.properties['run_cnt'] = len(self.properties['run_list'] )
        
        return 
        

    def proc_log(self,log_file,verbose=False,debug=False):
        """
        Read in results from lammps log file

        Args:
            * log_file (str) lammps log file

        """
        update_run = True
        run_cnt_i = 0


        logger.info("Reading log file {}".format(log_file))
        f = open(log_file,'r')
        log_lines = f.readlines()
        f.close()

        if( len(self.properties['run_list'] ) > 0 ):
            logger.info("Using existing run_list with {} runs ".format(len(self.properties['run_list'])))
            run_i = self.properties['run_list'][run_cnt_i]
        else:
            logger.info("No runs found will create new run_list ")
            self.properties['run_list']  = []
            run_i = MDrun()
            update_run = False
            
        thermo_keywords = ['Step','Temp','PotEng','TotEng','Press','Volume']
        for line in log_lines:
            col = line.split()  
            if( 'Loop time' in str(line) ):
                logger.info(" Calc  {}/{} finished with {} frames".format(run_cnt_i+1,len(self.properties['run_list']),run_i.n_frames ))
                run_cnt_i += 1
                if( update_run and len(self.properties['run_list']) > run_cnt_i):
                    run_i = self.properties['run_list'][run_cnt_i]
                    logger.info(" update_run is on setting run_i to index {} with {} steps".format(run_cnt_i,run_i.n_steps))
                else:
                    run_i = MDrun()
                    update_run = False 

            if( 'Step Temp PotEng TotEng Press Volume' not in line):

                if( len(col) >= 17 and col[0] != 'thermo_style' ):                    
                    if(  run_i.n_frames  == 2 ):
                        # Calculate dstep
                        run_i.dstep =  run_i.timeseries['step'][1] -  run_i.timeseries['step'][0]
                        logger.info("dstep %f "%(run_i.dstep ))
                    '''
                    elif(  run_i.properties['n_frames'] > 2 ):
                        # Test for new run 
                        new_run = False
                        dstep_i = int(col[0]) -  run_i.timeseries['step'][-1]
                        if( int(col[0]) < run_i.timeseries['step'][-1] ):
                            new_run = True
                            print ">LAMMPS.proc_log next step %d less than last %d "%(int(col[0]) , run_i.timeseries['step'][-1])
                        
                        elif( dstep_i != run_i.properties['dstep']  ):
                            new_run = True
                            print ">LAMMPS.proc_log dstep %f has changed %f at %f "%(dstep_i,run_i.properties['dstep', run_i.step_list[-1])
                        if( new_run ):
                            print "> LAMMPS.proc_log  new run found"
                            if( update_run ):
                                run_cnt_i += 1
                                run_i = self.properties['run_list'][run_cnt_i]
                            else:
                                self.properties['run_list'].append(copy.deepcopy(run_i))
                                run_i = MDrun()
                    '''
                    # 
                    # Add properties to timeseries
                    # 
                    for prop_i in run_i.timeseries.keys():
                        i = run_i.prop_col[prop_i]
                        run_i.timeseries[prop_i].append(float(col[i]))        
                    logger.debug("Frames:{}".format(run_i.n_frames ))

                    run_i.n_frames +=1
            else:
                logger.debug(" Adding thermo keys from line: %s "%(line))
                # Add thermo keys to properties
                llow = line.lower()
                col = llow.split()
                run_i.prop_col = dict()
                i = 0 
                for prop_i in col:
                    run_i.timeseries[prop_i] = []
                    run_i.prop_col[prop_i] = i
                    i += 1 
    
            #run_i.calc_time_list()
            #run_i.calc_density_list(self.strucC.properties['mass'])

            if( not update_run and  run_i.n_frames > 0 ):
                logger.info("Adding new run with {} frames ".format(run_i.n_frames ))
                self.properties['run_list'].append(copy.deepcopy(run_i))

    def analysis(self,output_key='log',data2cply=True):
        """
        Read in results from LAMMPS 
        """
        # Read .in if it exists
        in_key = 'in'
        try:
            in_file = self.files['input'][in_key]
            logger.info("Reading input file {}".format(in_file))
            if( self.resource.meta['type'] == "ssh" ):
                logger.debug("scp data from {}".format(self.resource.ssh['address']))
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])                              
                bash_command = "scp  %s:%s%s  ./ "%(ssh_id,self.dir['scratch'],in_file)
                os.system(bash_command)
                data2cply=False              
            self.proc_in(in_file,data2cply=data2cply)            
        except KeyError:
            logger.warning("Calculation {} No output_file file  with key {} found".format(self.tag,in_key))
        # Find output_key file 
        try:
            output_file = self.files['output'][output_key]
            logger.info("Reading output file {}".format(output_file))
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])                              
                bash_command = "scp  %s:%s%s  ./ "%(ssh_id,self.dir['scratch'],output_file)
                os.system(bash_command)                
            self.proc_log(output_file)            
        except KeyError:
            logger.warning("Calculation %s No output_file file  with key %s found"%(self.tag,output_key))
        
        #self.write_dat(dat_file)
        #self.add_file(dat_file,'data')


    def write_param( self,param_file,verbose = True, debug = False):
        """
        Write parameters to a file 
        """
        param_lines = "# Param file \n "
        param_lines += "\n "
        param_lines += " MASS (fftype,mass) \n "
        for ljtkey_i, ljtype_i  in self.paramC.ljtypes.iteritems():
            param_lines += " {}  {} \n".format( ljtype_i.fftype1 , ljtype_i.mass  )
            
        param_lines += "\n "
        param_lines += " NONBON (fftype,epsilon,sigma) \n"
        for ljtkey_i, ljtype_i  in self.paramC.ljtypes.iteritems():
            param_lines += " {}  {} {} \n".format( ljtype_i.fftype1,ljtype_i.epsilon, ljtype_i.sigma)

        if( self.paramC.n_bondtypes > 0 ):
            param_lines += "\n "
            param_lines += " BOND (fftype1,fftype2,type,...) \n"
            for btkey_i,bondtype_i  in self.paramC.bondtypes.iteritems():
                param_lines += "%s %s  %s"%(bondtype_i.fftype1, bondtype_i.fftype2, bondtype_i.type)
                if( bondtype_i.type == "harmonic"):
                    param_lines += " %f %f  \n"%(bondtype_i.kb,bondtype_i.r0)
                else:
                    logger.warning(" Unknown bond type %s "%(bondtype_i.type))

        if( self.paramC.n_angletypes > 0 ):
            param_lines += "\n "
            param_lines += " ANGLE (fftype1,fftype2,fftype3,type,...) \n"
            for atkey_i,angletype_i  in self.paramC.angletypes.iteritems():
                param_lines += "%s %s %s  %s"%(angletype_i.fftype1,angletype_i.fftype2,angletype_i.fftype3 ,angletype_i.type)
                if( angletype_i.type == "harmonic"):
                    param_lines += " %f %f  \n"%(angletype_i.kb,angletype_i.theta0)
                else:
                    logger.warning(" Unknown angle type %s "%(angletype_i.type))

        if( self.paramC.n_dihtypes > 0 ):
            param_lines += "\n "
            param_lines += " DIHEDRAL (fftype1,fftype2,fftype3,fftype4,type,...) \n"
            for dtkey_i, dihtype_i  in self.paramC.dihtypes.iteritems():    
                param_lines += "%s %s %s %s  %s"%(dihtype_i.fftype1,dihtype_i.fftype2,dihtype_i.fftype3,dihtype_i.fftype4 ,dihtype_i.type)
                if( dihtype_i.type == "multiharmonic"):
                    K = dihtype_i.kb
                    n = dihtype_i.mult
                    d = dihtype_i.theta_s
                    w = 0.0 # Weight
                    p = dihtype_i.paths
                    param_lines += " %f %d %f %d %f \n"%(K,n,d,w,p)
                elif( dihtype_i.type == "rb" or  dihtype_i.type == "opls"  ):
                    param_lines += " %f %f %f %f \n"%(dihtype_i.k1,dihtype_i.k2,dihtype_i.k3,dihtype_i.k4)
                else:
                    logger.warning(" Unknown DIHEDRAL type %s "%(dihtype_i.type))
                    
        # Write Dihedral Coeffs
        if( self.paramC.n_imptypes > 0 ):
            param_lines += "\n "
            param_lines += " IMPROPER (fftype1,fftype2,fftype3,fftype4,type,...) \n"
            for itkey_i, imptype_i  in self.paramC.imptypes.iteritems():
                param_lines += "%s %s %s %s  %s"%(imptype_i.fftype1,imptype_i.fftype2,imptype_i.fftype3,imptype_i.fftype4 ,imptype_i.type)
                if( imptype_i.type == "multiharmonic"):
                    param_lines += " %f %f \n"%(imptype_i.ke,imptype_i.e0)
                else:
                    logger.warning(" Unknown IMPROPER type %s "%(imptype_i.type))
                    print " Unknown IMPROPER type %s "%(imptype_i.type)
                    
        F = open(param_file , 'w' )
        F.write(param_lines)
        F.close()

            
    def read_param( self,param_file):
        """
        Read parameter file 


        Args:
            * param_file  (str) parameter file

        """
        debug = False
        # 
        F = open(param_file , 'r' )
        lines = F.readlines()
        F.close()
        #
        # Read in parameters 
        #
        #
        # Intialize
        #   - read in boolean to off
        #
        read_Masses  = False
        read_Pair  = False
        read_Bond_coeff  = False
        read_Angle_coeff  = False
        read_Dihedral_coeff  = False
        read_Improper_coeff  = False
        # 
        for line in lines:
            col = line.split()
            if( read_Masses ):
                if(  len(col) >= 2 ):
                    if( debug ):
                        print " Reading MASS ",col
                    fftype1 = str(col[0])
                    mass_i = float(col[1])
                    
                    el_i = periodictable.element_mass(mass_i)
                    ljtype_i = parameters.LJtype(fftype1)
                    ljtype_i.mass = mass_i
                    ljtype_i.atomic_symbol = el_i["symbol"]
                    ljtype_i.lammps_index = self.paramC.n_ljtypes + 1 
                    self.paramC.add_LJtype(ljtype_i)
                    
                    if( debug ):
                        print ljtype_i
                else:
                    if( debug ):
                        print " Finished reading MASS "
                    read_Masses = False

            if( read_Pair ):
                if( len(col) >= 3 ):
                    if( debug ):
                        print " Reading NONBON ",col
                    ljtype_i = self.paramC.ljtypes[cnt_Pair]
                    ljtype_i.epsilon = float(col[1])
                    ljtype_i.sigma = float(col[2])
                    cnt_Pair += 1
                    if( debug ):
                        print ljtype_i
                else:
                    if( debug ):
                        print " Finished reading NONBON "
                    read_Pair = False


            if( read_Bond_coeff ):
                if(  len(col) >= 3 ):
                    if( debug ):
                        print " Reading BOND ",col
                    fftype1 = str(col[0])
                    fftype2 = str(col[1])
                    btype = str(col[2])
                    
                    bondtype_i = parameters.Bondtype(fftype1,fftype2,type=btype)
                    if( btype == "harmonic"):
                        bondtype_i.kb = float(col[3])
                        bondtype_i.r0 = float(col[4])
                        bondtype_i.g_indx = int(1)
                    else:
                        logger.warning(" Unknown bond type %s "%(bondtype_i.type))                    
                    bondtype_i.lammps_index = self.paramC.n_bondtypes + 1 
                    self.paramC.add_bondtype(bondtype_i)                    
                    if( debug ):
                        print bondtype_i
                else:
                    # blank line delimited 
                    if( debug ):
                        print " Finished reading BOND "
                    read_Bond_coeff = False


            if( read_Angle_coeff ):
                if(  len(col) >= 3 ):
                    if( debug ):
                        print " Reading ANGLE ",col
                    fftype1 = str(col[0])
                    fftype2 = str(col[1])
                    fftype3 = str(col[2])
                    atype = str(col[3])
                    
                    angletype_i = parameters.Angletype(fftype1,fftype2,fftype3,type=atype)
                    if( atype == "harmonic"):
                        angletype_i.kb = float(col[4])
                        angletype_i.theta0 = float(col[5])
                        angletype_i.g_indx = int(1)
                    else:
                        logger.warning(" Unknown ANGLE type %s "%(angletype_i.type))                    
                    angletype_i.lammps_index = self.paramC.n_angletypes + 1 
                    self.paramC.add_angletype(angletype_i)
                    if( debug ):
                        print angletype_i

                else:
                    # blank line delimited 
                    read_Angle_coeff = False

            if( read_Dihedral_coeff ):
                if( len(col) >= 3 ):
                    if( debug ):
                        print " Reading DIHEDRAL ",col
                    fftype1 = str(col[0])
                    fftype2 = str(col[1])
                    fftype3 = str(col[2])
                    fftype4 = str(col[3])
                    dtype = str(col[4])
                    
                    dihtype_i = parameters.Dihtype(fftype1,fftype2,fftype3,fftype4,type=dtype)
                    if( dtype == "multiharmonic"):
                        dihtype_i.setharmonic(float(col[6]),float(col[5]),float(col[7]))
                        dihtype_i.g_indx = int(1)
                    elif( dtype == "rb" or  dtype == "opls"  ):
                        dihtype_i.setopls(float(col[5]),float(col[6]),float(col[7]),float(col[8]))
                        dihtype_i.g_indx = int(3)
                    else:
                        logger.warning(" Unknown DIHEDRAL type %s "%(dihtype_i.type))                    
                    dihtype_i.lammps_index = self.paramC.n_dihtypes + 1 
                    self.paramC.add_dihtype(dihtype_i)
                    if( debug ):
                        print dihtype_i

                else:
                    # blank line delimited 
                    read_Dihedral_coeff = False


            if( read_Improper_coeff ):
                if(  len(col) >= 3 ):
                    if( debug ):
                        print " Reading IMPROPER ",col
                    fftype1 = str(col[0])
                    fftype2 = str(col[1])
                    fftype3 = str(col[2])
                    fftype4 = str(col[3])
                    dtype = str(col[4])
                    
                    imptype_i = parameters.Imptype(fftype1,fftype2,fftype3,fftype4,type=dtype)
                    if( dtype == "multiharmonic"):
                        imptype_i.kb = float(col[5])
                        imptype_i.e0 = float(col[6])
                    else:
                        logger.warning(" Unknown IMPROPER type %s "%(imptype_i.type))                    
                    imptype_i.lammps_index = self.paramC.n_imptypes + 1 
                    self.paramC.add_imptype(imptype_i)
                    if( debug ):
                        print imptyp_i
                else:
                    # blank line delimited 
                    read_Improper_coeff = False

            if ( len(col) >= 1  ):
                if( col[0] == "MASS" ):
                    read_Masses = True
                    self.n_ljtypes = 0
                if( col[0] == "NONBON" ):
                    read_Pair = True
                    cnt_Pair = 0 
                if( col[0] == "BOND" ):
                    read_Bond_coeff = True
                    self.n_bondtypes = 0
                if( col[0] == "ANGLE" ):
                    read_Angle_coeff = True
                    self.n_angletypes  = 0

                if( col[0] == "DIHEDRAL" ):
                    read_Dihedral_coeff = True
                    self.n_dihtypes = 0

                if( col[0] == "IMPROPER" ):
                    read_Improper_coeff = True
                    self.n_imptypes = 0


        if( debug ):
            print "Finished read in of ",param_file

        
