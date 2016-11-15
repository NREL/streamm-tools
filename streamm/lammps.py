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


class mdrun():
    '''
    Object to store the output of a single MD run 
    '''
    def __init__(self, verbose=False):

        self.properties = dict()
        self.properties['timestep'] = 0.50  # Time step in fmsec
        self.properties['n_steps'] = 0 # Total steps 
        self.properties['n_frames'] = 0 # Total frames 
        self.properties['dstep'] = 1 # step rate
                
        # Create dictionary of lists for time series data 
        self.timeseries = dict()
        self.prop_col  =  dict() # list of properties in column order


          
class LAMMPS(CalculationRes):
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

        self.meta['software'] = 'lammps'
        self.units['distance'] = 'angstroms'
        self.units['energy'] = 'kcal/mol'
        self.units['time'] = 'ns'
        self.properties['finish_str'] = 'Loop time of'
        #
    def __del__(self):
        """
        Destructor, clears object memory
        """
        # Call base class destructor
        CalculationRes.__del__(self)

    def convert_gaussian(self,calc_gaussian):
        
        '''
        Convert Gaussian simulation to LAMMPS

        Input:
           calc_gaussian (Gaussian) simulation object
        Return:
           sim_lmp (LAMMPS) simulation object

        Units:
            Bohr -> angstroms
            Hartree -> kcal/mol 

        '''
        import gaussian
        
        if isinstance(calc_gaussian, gaussian.Gaussian):
            #
            # Convert Structure 
            #
            if( calc_gaussian.strucC.n_particles > 0):            
                # Add lattice to GROMACS object
                matrix_o = calc_gaussian.strucC.lat._matrix 
                matrix_i = self.strucC.lat._matrix 
                for m in range(calc_gaussian.strucC.lat.n_dim):
                    for n in range(calc_gaussian.strucC.lat.n_dim):
                        matrix_i[m][n] = units.convert_bohr_ang(matrix_o[m][n] )
                self.strucC.lat.set_matrix(matrix_i)
                #
                # Add particles and positions to LAMMPS object 
                #
                for pkey_o, particle_o  in calc_gaussian.strucC.particles.iteritems():
                    particle_i = copy.deepcopy(particle_o)
                    pos_o = calc_gaussian.strucC.positions[pkey_o]       
                    pos_i = [ units.convert_bohr_ang(v_i) for v_i in pos_o ] 
                    self.strucC.add_partpos(particle_i,pos_i)
                for bkey_o, bond_o  in calc_gaussian.strucC.bonds.iteritems():
                    bond_i = copy.deepcopy(bond_o)
                    self.strucC.add_bond(bond_i)
                for akey_o, angle_o in calc_gaussian.strucC.angles.iteritems():
                    angle_i = copy.deepcopy(angle_o)
                    self.strucC.add_angle(angle_i)
                for dkey_o,dih_o in calc_gaussian.strucC.dihedrals.iteritems():
                    dih_i = copy.deepcopy(dih_o)
                    self.strucC.add_dihedral(dih_i)
                for ikey_o,imp_o in calc_gaussian.strucC.impropers.iteritems():
                    imp_i = copy.deepcopy(imp_o)
                    self.strucC.add_improper(imp_i)

                    
            
    def convert_GROMACS(self,sim_gro):
        
        '''
        Convert GROMACS simulation to LAMMPS

        Input:
           sim_gro (GROMACS) simulation object
        Return:
           sim_lmp (LAMMPS) simulation object

        Units:
            nm -> angstroms
            kJ/mol -> kcal/mol 

        '''
        import gromacs
        if isinstance(sim_gro, gromacs.GROMACS):
            # Initialize LAMMPS simulation object 
            #
            # Convert Structure 
            #
            if( sim_gro.strucC.n_particles > 0):
                # Add lattice to LAMMPS object
                matrix_o = sim_gro.strucC.lat._matrix 
                matrix_i = self.strucC.lat._matrix 
                for m in range(sim_gro.strucC.lat.n_dim):
                    for n in range(sim_gro.strucC.lat.n_dim):
                        matrix_i[m][n] = units.convert_nm_angstroms(matrix_o[m][n] )
                self.strucC.lat.set_matrix(matrix_i)
                #
                # Add particles and positions to LAMMPS object 
                #
                for pkey_o, particle_o  in sim_gro.strucC.particles.iteritems():
                    particle_i = copy.deepcopy(particle_o)
                    pos_o = sim_gro.strucC.positions[pkey_o]       
                    pos_i = [ units.convert_nm_angstroms(v_i) for v_i in pos_o ] 
                    self.strucC.add_partpos(particle_i,pos_i)
                for bkey_o, bond_o  in sim_gro.strucC.bonds.iteritems():
                    bond_i = copy.deepcopy(bond_o)
                    self.strucC.add_bond(bond_i)
                for akey_o, angle_o in sim_gro.strucC.angles.iteritems():
                    angle_i = copy.deepcopy(angle_o)
                    self.strucC.add_angle(angle_i)
                for dkey_o,dih_o in sim_gro.strucC.dihedrals.iteritems():
                    dih_i = copy.deepcopy(dih_o)
                    self.strucC.add_dihedral(dih_i)
                for ikey_o,imp_o in sim_gro.strucC.impropers.iteritems():
                    imp_i = copy.deepcopy(imp_o)
                    self.strucC.add_improper(imp_i)

            #
            # Convert Parameters  
            #
            if( len(sim_gro.paramC.ljtypes) > 0):
                for ljtkey_o, ljtype_o  in sim_gro.paramC.ljtypes.iteritems():
                    ljtype_i = copy.deepcopy(ljtype_o)
                    ljtype_i.sigma = units.convert_nm_angstroms(ljtype_o.sigma)
                    ljtype_i.epsilon = units.convert_kJmol_kcalmol(ljtype_o.epsilon)
                    self.paramC.add_LJtype(ljtype_i)

                for btkey_o,bondtype_o  in sim_gro.paramC.bondtypes.iteritems():
                    bondtype_i =  copy.deepcopy(bondtype_o)
                    if( bondtype_o.g_indx == 1 ):
                        r0 = units.convert_nm_angstroms( bondtype_o.r0 )
                        kb = units.convert_g_bond_kb( bondtype_o.kb )
                        bondtype_i.setharmonic(r0,kb)
                    self.paramC.add_bondtype(bondtype_i)

                for atkey_o,angletype_o  in sim_gro.paramC.angletypes.iteritems():
                    angletype_i = copy.deepcopy(angletype_o)
                    if( angletype_o.g_indx == 1 ):
                        theta0 =  angletype_o.theta0
                        kb = units.convert_g_angle_kb( angletype_o.kb )
                        angletype_i.setharmonic(theta0,kb)
                    self.paramC.add_angletype(angletype_i)
                for dtkey_o, dihtype_o  in sim_gro.paramC.dihtypes.iteritems():
                    if( dihtype_o.g_indx == 1 or dihtype_o.g_indx == 9 ):
                        dihtype_i = copy.deepcopy(dihtype_o)
                        theta0 =  dihtype_o.theat_s
                        kb = units.convert_kJmol_kcalmol( dihtype_o.kb )
                        mult =  dihtype_o.mult
                        dihtype_i.setharmonic( mult, kb,theat_s)
                        self.paramC.add_dihtype(dihtype_i)

                    elif( dihtype_o.g_indx == 2 ):
                        imptyp_i = copy.deepcopy(dihtype_o)
                        e0 = imptyp_o.e0
                        ke =  units.convert_g_angle_kb( imptyp_o.ke )
                        imptyp_i.setimp(e0,ke)
                        self.paramC.add_omptype(imptyp_i)
                    elif( dihtype_o.g_indx == 3 ):
                        dihtype_i = copy.deepcopy(dihtype_o)
                        C0 = units.convert_kJmol_kcalmol( dihtype_o.C0 )
                        C1 = units.convert_kJmol_kcalmol( dihtype_o.C1 )
                        C2 = units.convert_kJmol_kcalmol( dihtype_o.C2 )
                        C3 = units.convert_kJmol_kcalmol( dihtype_o.C3 )
                        C4 = units.convert_kJmol_kcalmol( dihtype_o.C4 )
                        C5 = units.convert_kJmol_kcalmol( dihtype_o.C5 )           
                        dihtype_i.setrb(C0,C1,C2,C3,C4,C5)  # Sets oplsa as well since they are equivalent                     
                        self.paramC.add_dihtype(dihtype_i)
                    elif( dihtype_o.g_indx == 4 ):
                        imptyp_i = copy.deepcopy(dihtype_o)
                        e0 = imptyp_o.e0
                        ke =  units.convert_g_angle_kb( imptyp_o.ke )
                        imptyp_i.setimp(e0,ke)
                        imptyp_i.pn = imptyp_o.pn                    
                        self.paramC.add_omptype(imptyp_i)
                    elif( dihtype_o.g_indx == 5 ):
                        dihtype_i = copy.deepcopy(dihtype_o)
                        k1 = units.convert_kJmol_kcalmol( dihtype_o.k1 )
                        k2 = units.convert_kJmol_kcalmol( dihtype_o.k2 )
                        k3 = units.convert_kJmol_kcalmol( dihtype_o.k3 )
                        k4 = units.convert_kJmol_kcalmol( dihtype_o.k4 )
                        dihtype_i.setopls(k1,k2,k3,k4) # Sets oplsa as well since they are equivalent                     
                        self.paramC.add_dihtype(dihtype_i)                    

        else:
            raise TypeError("1st arg should be GROMACS simulation ")

    def read_data(self, data_file,
        btype = "harmonic",
        atype = "harmonic",
        dtype = "multiharmonic",
        imptype = "multiharmonic",
        debug = False):

        """
        Read Lammps data file

        Arguments:
            data_file  (str) data file
        Return:
            None
        """

        F = open(data_file , 'r' )
        lines = F.readlines()
        F.close()
        logger.debug(" Reading {} ",data_file)
        #
        # Read in data header with number of parameters 
        #
        matrix = self.strucC.lat._matrix 
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
                ljtype_i.lmpindx = int(col[0])
                ljtype_i.mass = mass_i
                ljtype_i.atomic_symbol = el_i["symbol"]

                #print ">read_data ljtype_i",ljkey_i,ljtype_i.mass,ljtype_i.lmpindx,cnt_Masses ,self.paramC.n_ljtypes 
                
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
                bondtype_i.lmpindx = int( col[0])
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
                angletype_i.lmpindx = int( col[0])
                self.paramC.angletypes[akey_i] = copy.deepcopy(angletype_i)
                                
                if( cnt_Angle_coeff >= self.paramC.n_angletypes ):
                    read_Angle_coeff = False


            if( read_Dihedral_coeff and  len(col) >= 3 ):
                cnt_Dihedral_coeff += 1
                dkey_i = int(col[0]) - 1

                #print "self.paramC.n_dihtypes ",self.paramC.n_dihtypes,cnt_Dihedral_coeff,dkey_i
                
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
                dihtype_i.lmpindx = int( col[0])
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
                imptype_i.lmpindx = int( col[0])
                self.paramC.imptypes[ikey_i] = copy.deepcopy(imptype_i)

                if( cnt_Improper_coeff >= self.paramC.n_imptypes ):
                    read_Improper_coeff = False



            if( read_Atoms and len(col) >= 7 ):

                # print ">read_Atoms col",col
                
                cnt_Atoms += 1
                
                pkey_i = int( col[0]) - 1
                ljkey_i = int( col[2]) - 1
                
                ljtype_i = self.paramC.ljtypes[ljkey_i]
                if( ljtype_i.lmpindx != ljkey_i +1 ):
                    logger.warning("Read in error for atom %s due to bad Pair type %d "%(cnt_Atoms,ljtype_i.lmpindx))
                    sys.exit(2)
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
                particle_i.properties["lmpindx"] = ljtype_i.lmpindx

                #print ">read_data cnt_Atoms ",cnt_Atoms,pkey_i,ljkey_i,ljtype_i.mass,ljtype_i.lmpindx 
                #if( cnt_Atoms > 10 ):
                #    sys.exit("9320ur092ur0298")
                
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
                bond_i.lmpindx = int(col[1])

                typekey =  bond_i.lmpindx  - 1
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
                angle_i.lmpindx = int(col[1])

                typekey =  angle_i.lmpindx  - 1
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
                dihedral_i.lmpindx = int(col[1])
                # Set parameter types 
                typekey =  dihedral_i.lmpindx  - 1
                try:
                    dihtype_i = self.paramC.dihtypes[typekey]
                    try:
                        dihtype_i.fftype1 = self.strucC.particles[pkey1].properties['fftype']
                        dihtype_i.fftype2 = self.strucC.particles[pkey2].properties['fftype']
                        dihtype_i.fftype3 = self.strucC.particles[pkey3].properties['fftype']
                        dihtype_i.fftype4 = self.strucC.particles[pkey4].properties['fftype']
                    except:
                        print " Particles %d %d has no ffytpe"%(pkey1,pkey2,pkey3,pkey4)
                except:
                    print "dihtype type  %d is not set "%(typekey)

                    
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
                improper_i.lmpindx = int(col[1])
                # Set parameter types 
                typekey =  improper_i.lmpindx  - 1
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
                    #    print ljtkey_i, ljtype_i.lmpindx , ljtype_i.mass ,  ljtype_i.fftype1,ljtype_i.epsilon,  ljtype_i.sigma 
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

        Arguments:
            data_file  (str) data file
        Return:
            None
        """

        F = open(data_file , 'r' )
        lines = F.readlines()
        F.close()
        logger.debug(" Reading {} ",data_file)
        #
        # Read in data header with number of parameters 
        #
        matrix = self.strucC.lat._matrix 
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
                    self.strucC.lat.set_matrix(matrix)
             
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
        F.write( "%10d  atom types \n" % self.paramC.n_ljtypes  )
        F.write( "%10d  bond types \n" % self.paramC.n_bondtypes )
        F.write( "%10d  angle types \n" % self.paramC.n_angletypes )
        F.write( "%10d  dihedral types \n" % self.paramC.n_dihtypes )
        if( self.paramC.n_imptypes > 0 ):
            F.write( "%10d  improper types \n" % self.paramC.n_imptypes )
        else:
            n=1
            F.write( "%10d  improper types \n" % n )
            
        F.write('\n')
        matrix = self.strucC.lat._matrix
        F.write( "%16.8f %16.8f   xlo xhi \n" %  (matrix[0][0]/-2.0 , matrix[0][0]/2.0) )
        F.write( "%16.8f %16.8f   ylo yhi \n" %  (matrix[1][1]/-2.0 , matrix[1][1]/2.0 ) )
        F.write( "%16.8f %16.8f   zlo zhi \n" %  (matrix[2][2]/-2.0 , matrix[2][2]/2.0) )
        F.write('\n')
        F.write( ' Masses \n')
        F.write('\n')
        # Write LJtypes mass 
        for ljtkey_i, ljtype_i  in self.paramC.ljtypes.iteritems():
            F.write( "%10d %16.8f   # %5s \n" % ( ljtype_i.lmpindx, ljtype_i.mass , ljtype_i.fftype1  ) )
        F.write('\n')
        F.write(' Pair Coeffs \n')
        F.write('\n')
        # Write LJtypes pair Coeffs 
        for ljtkey_i, ljtype_i  in self.paramC.ljtypes.iteritems():
            F.write( "%10d %12.6f %12.6f  \n" % (ljtype_i.lmpindx, ljtype_i.epsilon,  ljtype_i.sigma ) )
        F.write('\n')
        # Write Bond Coeffs
        if( self.paramC.n_bondtypes > 0 ):
            F.write(' Bond Coeffs \n')
            F.write('\n')
            for btkey_i,bondtype_i  in self.paramC.bondtypes.iteritems():
                if( bondtype_i.type == "harmonic"):
                    F.write( "%10d %12.6f %12.6f # %5s %5s  \n" % (bondtype_i.lmpindx,bondtype_i.kb,bondtype_i.r0, bondtype_i.fftype1, bondtype_i.fftype2 ) )
            F.write('\n')

        # Write Angle Coeffs
        if( self.paramC.n_angletypes > 0 ):
            F.write(' Angle Coeffs \n')
            F.write('\n')
            for atkey_i,angletype_i  in self.paramC.angletypes.iteritems():
                if( angletype_i.type == "harmonic"):
                    F.write( "%10d %12.6f %12.6f # %5s %5s  %5s   \n" % (angletype_i.lmpindx,angletype_i.kb,angletype_i.theta0, angletype_i.fftype1,angletype_i.fftype2,angletype_i.fftype3 ) )
            F.write('\n')

        # Write Dihedral Coeffs
        if( self.paramC.n_dihtypes > 0 ):
            F.write(' Dihedral Coeffs \n')
            F.write('\n')
            for dtkey_i, dihtype_i  in self.paramC.dihtypes.iteritems():    
                if( dihtype_i.type == "multiharmonic"):
                    # K = K[1+d cons(n theta)]
            
                    d = dihtype_i.theat_s
                    K = dihtype_i.kb
                    n = dihtype_i.mult
                    p = dihtype_i.paths
                    w = 0.0 # Weight 
                    F.write( "%10d %12.6f %12.6f  %12.6f %12.6f # %d x %5s %5s  %5s  %5s   \n" % (dihtype_i.lmpindx,K,n,d,w, p, dihtype_i.fftype1,dihtype_i.fftype2,dihtype_i.fftype3,dihtype_i.fftype4  ) )
                elif( dihtype_i.type == "rb" or  dihtype_i.type == "opls"  ):
                    # Get opls parameters                    
                    F.write( "%10d  %12.6f  %12.6f  %12.6f  %12.6f # %5s %5s  %5s %5s \n" % (dihtype_i.lmpindx,dihtype_i.k1,dihtype_i.k2,dihtype_i.k3,dihtype_i.k4, dihtype_i.fftype1,dihtype_i.fftype2,dihtype_i.fftype3,dihtype_i.fftype4  ) )
                else:
                    error_line = " Unknow dihedral type %s "%(dihtype_i.type )
                    sys.exit(error_line)
            F.write('\n')

        # Write Dihedral Coeffs
        if( self.paramC.n_imptypes > 0 ):
            
            F.write(' Improper Coeffs \n')
            F.write('\n')
            if( self.paramC.n_imptypes > 0 ):
                for itkey_i, imptype_i  in self.paramC.imptypes.iteritems():    
                    if( imptype_i.type == "improper"):
                        F.write( "%10d %12.6f %12.6f # %5s %5s  %5s  %5s   \n" % (imptype_i.lmpindx,imptype_i.ke,imptype_i.e0, imptype_i.fftype1,imptype_i.fftype2,imptype_i.fftype3,imptype_i.fftype4  ) )
                    else:
                        error_line = " Unknow improper type %s "%(imptype_i.type )
                        sys.exit(error_line)
                    
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
                fftype_i = particle_i.properties["fftype"]
                mol_i = particle_i.properties["mol"]
                charge_i = particle_i.properties["charge"]
                lmpindx_i = particle_i.properties["lmpindx"]
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
                AT_i =  self.strucC.particles[ bond_i.pkey1 ].properties["fftype"]
                AT_j =  self.strucC.particles[ bond_i.pkey2 ].properties["fftype"]
                #
                F.write(  '%9d %8d %9d %9d # %5s %5s \n' % (bkey_i+1,bond_i.lmpindx,b_i,b_j, AT_i, AT_j ) )
            F.write('\n')

        # Write Angles
        if( self.strucC.n_angles > 0 ):

            F.write(' Angles \n')
            F.write('\n')
            for akey_i, angle_i in self.strucC.angles.iteritems():
                a_k = angle_i.pkey1 + 1 
                a_i = angle_i.pkey2 + 1 
                a_j = angle_i.pkey3 + 1 
                AT_k = self.strucC.particles[ angle_i.pkey1 ].properties["fftype"]
                AT_i = self.strucC.particles[ angle_i.pkey2 ].properties["fftype"]
                AT_j = self.strucC.particles[ angle_i.pkey3 ].properties["fftype"]

                F.write(  '%9d %8d %9d %9d %9d  # %s %s %s \n' % (akey_i+1,angle_i.lmpindx,a_k,a_i,a_j,AT_k,AT_i,AT_j) )
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

                AT_k = self.strucC.particles[ dih_i.pkey1 ].properties["fftype"]
                AT_i = self.strucC.particles[ dih_i.pkey2 ].properties["fftype"]
                AT_j = self.strucC.particles[ dih_i.pkey3 ].properties["fftype"]
                AT_l = self.strucC.particles[ dih_i.pkey4 ].properties["fftype"]

                #error_line = " dih index is not in bounds "
                #error_line += " for atoms %d  %d  %d  %d  "%(d_k,d_i,d_j,d_l)
                #error_line += " for type atoms %s %s %s %s "%(AT_k,AT_i,AT_j,AT_l)
                #sys.exit(error_line)
                F.write(  '%9d %8d %9d %9d %9d %9d # %s %s %s %s \n' % (dkey_i+1,dih_i.lmpindx,d_k,d_i,d_j,d_l,AT_k,AT_i,AT_j,AT_l) )

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
                AT_k = self.strucC.particles[ imp_i.pkey1 ].properties["fftype"]
                AT_i = self.strucC.particles[ imp_i.pkey2 ].properties["fftype"]
                AT_j = self.strucC.particles[ imp_i.pkey3 ].properties["fftype"]
                AT_l = self.strucC.particles[ imp_i.pkey4 ].properties["fftype"]
                F.write(  '%9d %8d %9d %9d %9d %9d # %s %s %s %s \n' % (ikey_i+1,imp_i.lmpindx,d_k,d_i,d_j,d_l,AT_k,AT_i,AT_j,AT_l) )

            F.write( '\n' )            


        F.close()
        

    def proc_in(self,in_file,data2cply=True,verbose=False,debug=False):
        """
        Read in input parameters from lammps in file

        Args:
            in_file (str) lammps in file

        """
        # Add new properties 
        self.properties['run_cnt'] = 0
        self.properties['run_list'] = []
        
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
                if( debug ):
                    print ">analyze_in timestep" ,line,timestep_i
                
            if( 'run' in line):
                run_i = mdrun()
                run_i.timestep =  timestep_i
                run_i.n_steps = float(col[1])
                
                if( debug ):
                    print ">analyze_in run" ,line,run_i.timestep,run_i.n_steps
                
                self.properties['run_list'].append(copy.deepcopy(run_i))
                
            if( 'minimize' in line):
                run_i = mdrun()
                run_i.timestep =  timestep_i
                self.properties['run_list'].append(copy.deepcopy(run_i))
                
            if( 'write_data' in line):
                output_file = str(col[1])
                logger.info("Reading %s data file"%(output_file))
                self.add_file('output','data_%d'%(len(self.properties['run_list'])),output_file)
                if( data2cply ):
                    
                    self.read_data(output_file)
                    # Update particle properties with reference
                    if( ref_struc.n_particles == self.strucC.n_particles ):
                        # Set properties from cply file
                        for pkey_i, particle_i in self.strucC.particles.iteritems():
                            particle_i.properties =  ref_struc.particles[pkey_i].properties
                    else:
                        logger.warning(" Data file has %d particles and reference structure has %d particle "%(self.strucC.n_particles,ref_struc.n_particles))

                        
                    run_col = output_file.split('.')
                    
                    self.strucC.tag = run_col[0]
                    self.strucC.write_xyz()
                    self.add_file('output','xyz_%d'%(len(self.properties['run_list'])),"%s.xyz"%(self.strucC.tag))
                    self.strucC.write_cply()
                    self.add_file('output','cply_%d'%(len(self.properties['run_list'])),"%s.cply"%(self.strucC.tag))
                
                
            if( 'dump' in line and len(col) > 5 ):
                if( col[3] == 'dcd' ):
                    output_file = str(col[5])
                    self.add_file('output','dcd_%d'%(len(self.properties['run_list'])),output_file)

                
        # set run cnt for check()
        self.properties['run_cnt'] = len(self.properties['run_list'] )
        
        return 
        

    def proc_log(self,log_file,verbose=False,debug=False):
        """
        Read in results from lammps log file

        Args:
            log_file (str) lammps log file

        """

        f = open(log_file,'r')
        log_lines = f.readlines()
        f.close()

        update_run = False 
        if( len(self.properties['run_list'] ) > 0 ):
            update_run = True
            run_cnt_i = 0
            run_i = self.properties['run_list'][run_cnt_i]
            print ">  using mdrun with len %d "%(len(self.properties['run_list']))
        else:
            self.properties['run_list']  = []
            run_i = mdrun()
        if( debug):
            print ">analyze_log update_run",update_run
        thermo_keywords = ['Step','Temp','PotEng','TotEng','Press','Volume']
        for line in log_lines:
            # print "line:",line
            # llow = line.lower()
            col = line.split()
  
            if( 'Loop time' in str(line) ):
                print " Calc %s finished "%(run_cnt_i)
                run_cnt_i += 1
                if( update_run and len(self.properties['run_list']) <= run_cnt_i):
                    run_i = self.properties['run_list'][run_cnt_i]
                else:
                    self.properties['run_list'].append(copy.deepcopy(run_i))
                    run_i = mdrun()                    

            if( 'Step Temp PotEng TotEng Press Volume' not in line):
                
                
            

                if( len(col) >= 17 and col[0] != 'thermo_style' ):
                    print "> col ",run_i.properties['n_frames'],col
                    
                    if(  run_i.properties['n_frames']  == 2 ):
                        # Calculate dstep
                        print  run_i.timeseries['step']
                        run_i.properties['dstep']  =  run_i.timeseries['step'][1] -  run_i.timeseries['step'][0]
                        if( debug):
                            print ">LAMMPS.analyze_log run_i. dstep %f "%(run_i.properties['dstep'] )
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
                                run_i = mdrun()
                    '''
                        
                    if( debug ):
                        print "> read_enline ",len(col),col
                    
                    for prop_i in run_i.timeseries.keys():
                        i = run_i.prop_col[prop_i]
                        run_i.timeseries[prop_i].append(float(col[i]))

                    run_i.properties['n_frames'] +=1
                    print run_i.properties['n_frames'],"9823u ",run_i.timeseries['step']
                  
            else:
                print " Adding thermo keys from line: %s "%(line)
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

        if( not update_run ):
            self.properties['run_list'].append(copy.deepcopy(run_i))

    def analysis(self,output_key='log'):
        """
        Read in results from LAMMPS 
        """
        # Read .in if it exists
        in_key = 'in'
        try:
            in_file = self.files['input'][in_key]
            if( self.resource.meta['type'] == "ssh" ):
                ssh_id = "%s@%s"%(self.resource.ssh['username'],self.resource.ssh['address'])                              
                bash_command = "scp  %s:%s%s  ./ "%(ssh_id,self.dir['scratch'],in_file)
                os.system(bash_command)
                data2cply=False              
            self.proc_in(in_file,data2cply=data2cply)            
        except KeyError:
            print "Calculation %s No output_file file  with key %s found"%(self.tag,in_key)
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
                    d = dihtype_i.theat_s
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


        Arguments:
            param_file  (str) parameter file

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
                    ljtype_i.lmpindx = self.paramC.n_ljtypes + 1 
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
                    bondtype_i.lmpindx = self.paramC.n_bondtypes + 1 
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
                    angletype_i.lmpindx = self.paramC.n_angletypes + 1 
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
                    dihtype_i.lmpindx = self.paramC.n_dihtypes + 1 
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
                    
                    imptype_i = parameters.Dihtype(fftype1,fftype2,fftype3,fftype4,type=dtype)
                    if( dtype == "multiharmonic"):
                        imptype_i.kb = float(col[5])
                        imptype_i.e0 = float(col[6])
                    else:
                        logger.warning(" Unknown IMPROPER type %s "%(imptype_i.type))                    
                    imptype_i.lmpindx = self.paramC.n_imptypes + 1 
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

        
