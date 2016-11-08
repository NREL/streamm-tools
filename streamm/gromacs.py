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
          
class GROMACS(CalculationRes):
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

        self.meta['software'] = 'gromacs'
        self.units['distance'] = 'nm'
        self.units['energy'] = 'kJ/mol'
        #
    def __del__(self):
        """
        Destructor, clears object memory
        """
        # Call base class destructor
        CalculationRes.__del__(self)


    def add_strucC(self,strucC_i):
        '''
        Add a structure.Container to the simulation 
        '''

        if isinstance(strucC_i, structure.Container):
            strucC_add = copy.deepcopy(strucC_i)
            # Convert structure from default units
            #   angstorms -> nm
            matrix = strucC_add.lat._matrix
            for m in range(strucC_add.lat.n_dim):
                for n in range(strucC_add.lat.n_dim):
                    matrix[m][n] = units.convert_angstroms_nm(matrix[m][n] )
            strucC_add.lat.set_matrix(matrix)
            for pkey_i, particle_i  in strucC_add.particles.iteritems():
                pos_i = strucC_add.positions[pkey_i]       
                pos_i_nm = [units.convert_angstroms_nm(pos_i[0]) ,units.convert_angstroms_nm(pos_i[1]) ,units.convert_angstroms_nm(pos_i[2]) ]
                strucC_add.positions[pkey_i]  = pos_i_nm

            self.strucC += strucC_add   
        else:
            raise TypeError("1st arg should be a structure.Container object")



    def write_xyz(self, xyz_file=''):
        '''
        Write a structure  to an xyz file

        Args:
            xyz_file    (str) xyz file tag
            
        Reutrns:
            null
            
        '''
        if( len(xyz_file) == 0 ):
            xyz_file = "%s.xyz"%(self.tag)

        F = open(xyz_file,"w")

        # Loop over structures
        F.write(" %d \n" % self.strucC.n_particles )
        F.write(" %s \n"%" structure.Container %s  "%(self.tag))
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            pos_o = self.strucC.positions[pkey_i]
            pos_i = [ units.convert_nm_angstroms(v_i) for v_i in pos_o ]
            F.write( " %5s %16.8f %16.8f %16.8f \n"  % (particle_i.properties["symbol"] ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]) ) )
        F.close()
        
    def write_gro(self,gro_file=''):
        """
        Write gromacs structure file 
        """
        #
        if( len(gro_file) == 0 ):
            gro_file = "%s.gro"%(self.tag)
        
        matrix = self.strucC.lat._matrix 

        gro_lines =  " com2gro \n"
        gro_lines += " %-2i \n"  %( self.strucC.n_particles )
        atom_indx = 0 
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            atom_indx += 1
            if( atom_indx > 10000): atom_indx = 1
            pos_i = self.strucC.positions[pkey_i]       
            gro_lines += "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"  % ( particle_i.properties["residue"],particle_i.properties["resname"][:5],particle_i.properties["label"][:5],atom_indx,pos_i[0],pos_i[1],pos_i[2] )
        gro_lines += " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (matrix[0][0],matrix[1][1],matrix[2][2],matrix[0][1],matrix[0][2],matrix[1][0],matrix[1][2],matrix[2][0],matrix[2][1]) 

        F = open( gro_file, 'w' )
        F.write(gro_lines)
        F.close()


    def read_itp(self,ff_file):
        """
        Read a gromacs paramter file

        Args:
            parmC (parameter.Container)  to add parameters to
            ff_file (str) GROMACS itp file
        Return:
            parmC (parameter.Container) with read in parameters 
        """

        F = open(ff_file , 'r' )
        Lines = F.readlines()
        F.close()
        #
        # Get DEFAULTS
        #
        read_DEFAULTS = False
        for line in Lines:
            col = line.split()
            if ( len(col) > 1 ):
                if( read_DEFAULTS ):
                    if ( len(col) >= 2 ):
                        if ( col[0][:1] != ';' and  col[0][:1] != '[' ):
                            col_int = 0

                            self.paramC.nbfunc = int( col[col_int] )
                            self.paramC.combmixrule = int( col[col_int+1] )
                            self.paramC.genpairs = str( col[col_int+2] )
                            self.paramC.fudgeLJ = float( col[col_int+3] )
                            self.paramC.fudgeQQ = float( col[col_int+4] )

                    if ( col[0][:1] == '[' ):
                        read_DEFAULTS = False

                if ( col[0][:1] == '[' and col[1] == 'defaults' ):
                    read_DEFAULTS = True
        #
        # Get atom types
        # 
        read_line = False 
        ID='atomtypes'
        for line in Lines:
            col = line.split()
            if ( len(col) > 2 ):
                if ( col[0][:1] == '[' and col[1] == ID ):
                    read_line = True 
                if ( col[0][:1] == '[' and col[1] != ID ):
                    read_line = False 
                if read_line :
                    if ( len(col) >= 7 ):
                        if ( col[0][:1] != ';' ):

                            col_int = 0 
                            if(  len(col[0]) > 2 ):
                                col_int = 1

                            fftype1 = col[col_int]
                            ljtype_i = parameters.LJtype(fftype1)
                            ljtype_i.mass = float( col[col_int+2] )
                            ljtype_i.sigma = float(col[col_int+5])
                            ljtype_i.epsilon  = float(col[col_int+6])
                            self.paramC.add_LJtype(ljtype_i)

        #
        # Get bond types
        # 
        read_line = False
        ID='bondtypes'
        for line in Lines:
            col = line.split()
            if ( len(col) > 2 ):
                if ( col[0][:1] == '[' and col[1] == ID ):
                    read_line = True
                if ( col[0][:1] == '[' and col[1] != ID ):
                    read_line = False
                if read_line :
                    if ( len(col) >= 5 ):
                        if ( col[0][:1] != ';' ):
                            fftype1 = col[0]
                            fftype2 = col[1]
                            # Set parameters according to type 
                            g_type = int( col[2] )                    # Gromacs id 
                            if( g_type == 1 ):
                                r0 =  float(col[3]) 
                                kb =  float(col[4]) 
                                btype = "harmonic"
                            btyp_i = parameters.Bondtype(fftype1,fftype2,btype)
                            btyp_i.g_indx =  g_type

                            if( g_type == 1 ):
                                btyp_i.setharmonic(r0,kb)
                            self.paramC.add_bondtype(btyp_i)

        #
        # Get angle types
        # 
        read_line = False
        ID='angletypes'
        for line in Lines:
            col = line.split()
            if ( len(col) > 2 ):
                if ( col[0][:1] == '[' and col[1] == ID ):
                    read_line = True
                if ( col[0][:1] == '[' and col[1] != ID ):
                    read_line = False
                if read_line :
                    if ( len(col) >= 5 ):
                        if ( col[0][:1] != ';' ):
                            fftype1 = col[0]
                            fftype2 = col[1]
                            fftype3 = col[2]
                            atype = "blank"
                            # Set parameters according to type 
                            gfunc_type = int( col[3] )                    # Gromacs id 
                            if( gfunc_type == 1 ):
                                theta0 = float( col[4] )        # degrees 
                                kb = float( col[5] ) 
                                atype = "harmonic"
                            atyp_i = parameters.Angletype(fftype1,fftype2,fftype3,atype)
                            atyp_i.g_indx =  gfunc_type

                            if( gfunc_type == 1 ):
                                atyp_i.setharmonic(theta0,kb)

                            self.paramC.add_angletype(atyp_i)

        #
        # Get dihedral types
        # 
        read_line = False
        ID='dihedraltypes'
        for line in Lines:
            col = line.split()
            if ( len(col) > 2 ):
                if ( col[0][:1] == '[' and col[1] == ID ):
                    read_line = True
                if ( col[0][:1] == '[' and col[1] != ID ):
                    read_line = False

                if read_line :
                    if ( len(col) > 6 ):
                        if ( col[0][:1] != ';' and col[0][:1] != '#' ):
                            fftype1 = col[0]
                            fftype2 = col[1]
                            fftype3 = col[2]
                            fftype4 = col[3]
                            # Set parameters according to type
                            gfunc_type = int( col[4] )                    # Gromacs id

                            if( gfunc_type == 1 or gfunc_type == 9 ):
                                theat_s = float( col[5] )
                                kb = float( col[6] ) 
                                mult = float( col[7] )
                                dtype = "multiharmonic"
                                # Vd(theta) = kb[1 + cos(mult theta - theat_s)]

                            elif( gfunc_type == 2 ):
                                e0 = float( col[5] )
                                ke =  float( col[6] ) 
                                dtype = "improper"
                                # V = 1/2 ke( e-e0)^2

                            elif( gfunc_type == 3 ):
                                C0 =  float( col[5] ) 
                                C1 =  float( col[6] ) 
                                C2 =  float( col[7] ) 
                                C3 =  float( col[8] ) 
                                C4 =  float( col[9] ) 
                                C5 =  float( col[10] ) 
                                dtype = "rb"
                                # Ryckaert-Bellemans function
                                # Vrb(theta) = \sum_n=0^5 C_n [ cos(theata - 180 )^n ]

                            elif( gfunc_type == 4 ):
                                e0 = float( col[5] )
                                ke = float( col[6] ) 
                                pn = int( col[7] )
                                dtype = "periodicimproper"

                            elif( gfunc_type == 5 ):
                                k1 =  float( col[5] )
                                k2 =  float( col[6] )
                                k3 =  float( col[7] )
                                k4 =  float( col[8] )
                                dtype = "opls"
                                # opls function

                            else:
                                raise TypeError(" unknow dihedral type %d "%gfunc_type)

                            if( gfunc_type == 2 ):
                                imptyp_i = parameters.Imptype(fftype1,fftype2,fftype3,fftype4,dtype)
                                imptyp_i.g_indx =  gfunc_type
                                imptyp_i.setimp(e0,ke)
                                self.paramC.add_imptype(imptyp_i)
                            elif( gfunc_type == 4 ):
                                imptyp_i = parameters.Imptype(fftype1,fftype2,fftype3,fftype4,dtype)
                                imptyp_i.g_indx =  gfunc_type
                                imptyp_i.setimp(e0,ke)
                                imptyp_i.set_pn(pn)
                                self.paramC.add_imptype(imptyp_i)
                            elif( gfunc_type == 1 or  gfunc_type == 9  ):
                                dtyp_i = parameters.Dihtype(fftype1,fftype2,fftype3,fftype4,dtype)
                                dtyp_i.g_indx =  gfunc_type
                                dtyp_i.setharmonic( mult, kb,theat_s)
                                self.paramC.add_dihtype(dtyp_i)
                            elif( gfunc_type == 3 ):
                                dtyp_i = parameters.Dihtype(fftype1,fftype2,fftype3,fftype4,dtype)
                                dtyp_i.g_indx =  gfunc_type
                                dtyp_i.setrb(C0,C1,C2,C3,C4,C5)  # Sets oplsa as well since they are equivalent 
                                self.paramC.add_dihtype(dtyp_i)
                            elif(  gfunc_type == 5  ):
                                dtyp_i = parameters.Dihtype(fftype1,fftype2,fftype3,fftype4,dtype)
                                dtyp_i.g_indx =  gfunc_type
                                dtyp_i.setopls(k1,k2,k3,k4)   # Sets rb as well since they are equivalent 
                                self.paramC.add_dihtype(dtyp_i)
                            else:
                                raise TypeError(" unknow dihedral type %d "%gfunc_type)

