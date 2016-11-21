#! /usr/bin/env python
"""
This module defines the classes relating to classical particles in box

The default units are

distance - Angstroms 
mass - AMU

Units are changed when a structure is associated with a simulation 

"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"


import numpy as np
from numpy import pi #, dot, transpose, radians
import copy, math, random , sys, json, csv, os
try:
    import cPickle as pickle
except:
    import pickle
    
from datetime import datetime

import periodictable

import logging
logger = logging.getLogger(__name__)


class Particle():
    """
    Data structure for describing any localized object in QM/MD simulations
    A 'Particle' has a type and dict of specifiers to it's property 
    """

    def __init__(self, type="blank"):
        """
        Constructor for a general particle. 
        """
        self.type = type
        self.tag = "blank"
        # Tags dictionary. To be set by caller
        self.properties=dict()

    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.type
        del self.tag
        del self.properties

        
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s "%(self.type)

# Then you have types particles
class Atom(Particle):
    """
    A derived type of particles for atoms 
    """

    def __init__(self,symbol="X", type="atom"):
        Particle.__init__(self, type=type)

        # Get properties of element based on symbol 
        self.properties = periodictable.element_symbol(symbol)

        self.properties["mol"] = 0  
        self.properties["charge"] = 0.0     
        self.properties["fftype"] = self.properties["symbol"]
        self.properties["ffmass"] = self.properties["mass"]
        self.properties["lmpindx"] = -1 
        self.properties["group"] = 0
        self.properties["ring"] = 0
        self.properties["residue"] = 0
        self.properties["resname"] = "RES"
        self.properties["qgroup"] = 0
        self.properties["label"] =  self.properties["symbol"]
        
        self.tag = self.properties["symbol"]
        
class Bond:
    """
    Data structure for describing any 2-point associatiaon of Particle-s
    """

    def __init__(self, pkey1, pkey2):
        """
        Constructor for a general bond. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pkey1   (int)   Dictionary key of Particle object in bond
            pkey2   (int)   Dictionary key of Particle object in bond
        """
        
        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            print "2nd arg should be int type"
            raise TypeError

        self.lmpindx = 0 
        self.g_indx = 0

        self.properties = dict()
        # self.border = 1 # 1-single,2-double,3-triple

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2
        del self.lmpindx
        del self.g_indx
        del self.properties

    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s  "%(self.pkey1,self.pkey2 )



class Angle:
    """
    Data structure for describing any 3-point associatiaon of Particle-s
    """

    def __init__(self, pkey1, pkey2, pkey3):
        """
        Constructor for a general angle. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pkey1   (int)   Dictionary key of Particle object in angle
            pkey2   (int)   Dictionary key of Particle object in angle
            pkey3   (int)   Dictionary key of Particle object in angle
        """
        
        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            print "2nd arg should be int type"
            raise TypeError

        if isinstance(pkey3, int):
            self.pkey3 = pkey3
        else:
            print "3rd arg should be int type"
            raise TypeError

        self.lmpindx = 0
        self.g_indx = 0
        self.properties = dict()
        
    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2
        del self.pkey3
        del self.lmpindx
        del self.g_indx
        del self.properties


    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s - %s  "%(self.pkey1,self.pkey2,self.pkey3)

class Dihedral:
    """
    Data structure for describing any 4-point associatiaon of Particle-s
    """

    def __init__(self, pkey1, pkey2, pkey3, pkey4):
        """
        Constructor for a general dihedral. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pkey1   (int)   Dictionary key of Particle object in dihedral
            pkey2   (int)   Dictionary key of Particle object in dihedral
            pkey3   (int)   Dictionary key of Particle object in dihedral
            pkey4   (int)   Dictionary key of Particle object in dihedral
        """
        
        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            print "2nd arg should be int type"
            raise TypeError

        if isinstance(pkey3, int):
            self.pkey3 = pkey3
        else:
            print "3rd arg should be int type"
            raise TypeError

        if isinstance(pkey4, int):
            self.pkey4 = pkey4
        else:
            print "4rd arg should be int type"
            raise TypeError

        self.lmpindx = 0
        self.g_indx = 0
        
        self.properties = dict()


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2 
        del self.pkey3 
        del self.pkey4
        del self.lmpindx
        del self.g_indx 
        del self.properties

    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s - %s - %s "%(self.pkey1,self.pkey2,self.pkey3,self.pkey4 )


class Improper:
    """
    Data structure for describing any 4-point associatiaon of Particle-s
    """

    def __init__(self, pkey1, pkey2, pkey3, pkey4):
        """
        Constructor for a general improper. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pkey1   (int)   Dictionary key of Particle object in improper
            pkey2   (int)   Dictionary key of Particle object in improper
            pkey3   (int)   Dictionary key of Particle object in improper
            pkey4   (int)   Dictionary key of Particle object in improper
        """

        if isinstance(pkey1, int):
            self.pkey1 = pkey1
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(pkey2, int):
            self.pkey2 = pkey2
        else:
            print "2nd arg should be int type"
            raise TypeError

        if isinstance(pkey3, int):
            self.pkey3 = pkey3
        else:
            print "3rd arg should be int type"
            raise TypeError

        if isinstance(pkey4, int):
            self.pkey4 = pkey4
        else:
            print "4rd arg should be int type"
            raise TypeError
        #
        self.lmpindx = 0
        self.g_indx = 0
        #
        self.properties = dict()

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pkey1
        del self.pkey2 
        del self.pkey3 
        del self.pkey4
        del self.lmpindx
        del self.g_indx 
        del self.properties

    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s - %s - %s "%(self.pkey1,self.pkey2,self.pkey3,self.pkey4 )
        
class Lattice():
    
    def __init__(self):
        """        
        Create a lattice 
        """
        self.n_dim = 3 # Number of spcail dimensions is set to 3
        self._matrix = np.zeros((self.n_dim,self.n_dim)) #np.array(matrix, dtype=np.float64).reshape((3, 3))
        
        self._angles = np.zeros(self.n_dim)
        self._lengths = np.zeros(self.n_dim)
        self.pbcs = [ False for d in range(self.n_dim) ] # Initialize periodic boundries as off


    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """

        del self._angles
        del self._lengths
        del self._matrix
                
    def set_matrix(self,matrix):
        '''
        Set the matrix values and reset lengths and angles accordingly 
        '''

        def abs_cap(val, max_abs_val=1):
            """
            Returns the value with its absolute value capped at max_abs_val.
            Particularly useful in passing values to trignometric functions where
            numerical errors may result in an argument > 1 being passed in.

            Args:
                val (float): Input value.
                max_abs_val (float): The maximum absolute value for val. Defaults to 1.

            Returns:
                val if abs(val) < 1 else sign of val * max_abs_val.
            """
            return max(min(val, max_abs_val), -max_abs_val)

        m = np.array(matrix, dtype=np.float64).reshape((self.n_dim, self.n_dim))
        self._matrix = m
        self._lengths = np.sqrt(np.sum(m ** 2, axis=1))
        self._angles = [ 90.0 for d in range(self.n_dim) ]
        if( sum(self._lengths) > 0.0 ):
            angles = np.zeros(self.n_dim)
            for i in range(self.n_dim):
                j = (i + 1) % self.n_dim
                k = (i + 2) % self.n_dim
                self._angles[i] = abs_cap(np.dot(m[j], m[k]) / (self._lengths[j] * self._lengths[k]))

            self._angles = np.arccos(angles) * 180. / pi

        self.pbcs = [ True for d in range(self.n_dim) ] # Initialize periodic boundries as off


    def LCtoLV(self,box):
        """
        Convert lattice constants  to lattice  vectors
        """
        A = box[0]
        B = box[1]
        C = box[2]
        alpha = np.deg2rad( box[3] )
        beta  = np.deg2rad( box[4] )
        gamma = np.deg2rad( box[5] )

        ax = A
        bx = math.cos( gamma )* B
        by = math.sin( gamma )* B
        cx = C*math.cos(beta)
        cy = C*(math.cos(alpha)- math.cos(beta)*math.cos(gamma))/math.sin(gamma)
        cz = np.sqrt( C*C - cx*cx - cy*cy )

 
        self._lengths[0] = A
        self._lengths[1] = B
        self._lengths[2] = C

        self._angles[0] = alpha
        self._angles[1] = beta
        self._angles[2] = gamma
        
        
        # Set lattice vectors
        self._matrix =  np.array([ np.zeros(self.n_dim) for dim in range(self.n_dim) ])
        self._matrix[0][0] = ax
        self._matrix[1][0] = bx
        self._matrix[1][1] = by
        self._matrix[2][0] = cx
        self._matrix[2][1] = cy
        self._matrix[2][2] = cz

        return 

    def deltasq_pos(self,pos_i,pos_j):
        """
        Difference between two positions 
        """
        
        dr_ij  = np.array(pos_j) -  np.array(pos_i)

        return  dr_ij

    def deltasq_pos_c(self,pos_i,pos_j):
        """
        Difference between two positions in cubic lattice  
        """
        dr_ij = self.deltasq_pos(pos_i,pos_j)
        
        dr_ij[0] = dr_ij[0] - self._matrix[0][0] * round( dr_ij[0]/  self._matrix[0][0] )
        dr_ij[1] = dr_ij[1] - self._matrix[1][1] * round( dr_ij[1]/  self._matrix[1][1] )
        dr_ij[2] = dr_ij[2] - self._matrix[2][2] * round( dr_ij[2]/  self._matrix[2][2] )
        
        return dr_ij

    def norm_delta_pos_c(self,pos_i,pos_j):
        '''
        Normalized difference between two positions in cubic lattice 
        '''
        dr_ij = self.deltasq_pos_c(pos_i, pos_j)
        return (dr_ij)/np.linalg.norm(dr_ij)
        
        

    def delta_pos_c(self,pos_i,pos_j):
        """
        Difference between two positions in cubic lattice  
        """
        dr_ij = self.deltasq_pos_c(pos_i,pos_j)
        mag_dr_ij = np.sqrt(dr_ij.dot(dr_ij))

        return dr_ij,mag_dr_ij


    def delta_pos(self,pos_i,pos_j):
        """
        Difference between two positions 
        """
        
        dr_ij  = self.deltasq_pos(pos_i,pos_j)
        mag_dr_ij = np.sqrt(dr_ij.dot(dr_ij))

        return dr_ij,mag_dr_ij


    def proximitycheck(self,npos_i,npos_j,pos_cut,p=None):
        """
        Difference between two position lists in  cubic lattice  
        """
        n_i = len(npos_i)
        n_j = len(npos_j)
        key_list_i = range(n_i)
        key_list_j = range(n_j)
        pos_cut_sq = pos_cut*pos_cut
        #
        # MPI setup
        #
        if( p == None ):
            rank = 0
            size = 0
            key_list_i_p = key_list_i
        else:
            rank = p.getRank()
            size = p.getCommSize()
            key_list_i_p =  p.splitListOnProcs(key_list_i)
            
        overlap_p = 0 
        # if any pbc's are on 
        if( any ( pbcs_i for pbcs_i in self.pbcs ) ):
            logging.debug(" Using cubic pbcs")
            for m in key_list_i_p:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.deltasq_pos_c(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    if( dot_dr_ij < pos_cut_sq ):
                        overlap_p += 1 
        else:
            logging.debug(" Using no pbcs")
            for m in key_list_i_p:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_i[n]
                    dr_ij = self.deltasq_pos(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    if( dot_dr_ij < pos_cut_sq ):
                        overlap_p += 1 
        if( p == None ):
            if( overlap_p == 0 ):
                return True
            else:
                return False
        else:
            overlap_sum = p.allReduceSum(overlap_p)
            if( overlap_sum == 0 ):
                return True
            else:
                return False
                                
                            
    def delta_npos(self,npos_i,npos_j):
        """
        Difference between two position lists in  cubic lattice  
        """
        n_i = len(npos_i)
        n_j = len(npos_j)
        key_list_i = range(n_i)
        key_list_j = range(n_j)
        #

        n_ij = n_i*n_j
        #print n_i,n_j

        start_dpos = datetime.now()
        logging.debug(" Taking difference %d x %d "%(n_i,n_j))
        
        npos_ij = np.empty(n_ij*self.n_dim,dtype='float64')
        #nd_ij = np.empty(n_ij,dtype='float64')
        nd_ij=  np.zeros(shape=(n_i,n_j),dtype='float64')

        # if any pbc's are on 
        if( any ( pbcs_i for pbcs_i in self.pbcs ) ):
            logging.debug(" Using cubic pbcs")
            for m in key_list_i:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.deltasq_pos_c(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    #npos_ij[m][n] = dr_ij
                    nd_ij[m][n] = np.sqrt(dot_dr_ij)
        else:
            logging.debug(" Using no pbcs")
            for m in key_list_i:
                pos_i = npos_i[m]
                for n in key_list_j:
                    pos_j = npos_j[n]
                    dr_ij = self.deltasq_pos(pos_i,pos_j)
                    dot_dr_ij = dr_ij.dot(dr_ij)
                    #npos_ij[m][n] = dr_ij
                    nd_ij[m][n] = np.sqrt(dot_dr_ij)

        finish_dpos = datetime.now()
        delt_t = finish_dpos - start_dpos
        logging.debug(" Finished taking difference %s sec "%(delt_t.seconds))

          
        return npos_ij,nd_ij

    def random_pos(self):
        '''
        Generate random position in lattice

        random.seed(seed) need to be initialized
        
        '''
        
        n_dec = int(6)
        int_mult = 10**float(n_dec)
        pos_o = np.zeros(self.n_dim)
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                if( self._matrix[m][n] > 0.0 ):
                    pos_o[m] += float(random.randrange(0,int(int_mult*self._matrix[m][n]) ))/float(int_mult)
        return pos_o


    def expand_matrix(self,exlat_frac):
        '''
        Increase size of lattice by certain fraction

        Arguments:
        exlat_frac (float) fraction to increase lattice by
        '''
        matrix_i = self._matrix
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                matrix_i[m][n] += matrix_i[m][n]*exlat_frac
                
        self.set_matrix(matrix_i)
        
        return 


    def fractoreal(self,frac_o):
        '''
        Translate fractional coordinates to real 

        Arguments:
            frac_o (np.array) fraction coordinates
        '''
        pos_o = np.zeros(self.n_dim)
        for m in range(self.n_dim):
            for n in range(self.n_dim):
                pos_o[m] += self._matrix[n][m]*frac_o[n]
                        
        return pos_o
                                
class NBlist():
    """
    Class for neighbor list
    """
    def __init__(self,verbose=False):
        """
        Constructor
        """
        self.verbose = verbose
        self.list = []
        self.index = []
        self.cnt = -1 

    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        del self.verbose
        del self.list
        del self.index
        del self.cnt

    

    def getnbs(self,key_i):
        '''
        Return list of keys of neighbors of key_i
        '''
        try:
            nbs_i = self.list[self.index[key_i]:self.index[key_i+1]]
        except:
            logger.warning(" Neighbor list not set ")
            nbs_i = []
            
        return nbs_i

    def calc_nnab(self,key_i):
        '''
        Calculate the number of neighbors for particle key_i
        '''
        try:
            nnab_i = self.index[key_i+1] - self.index[key_i]
        except:
            logger.warning(" Neighbor list not set ")
            nnab_i = 0
            
        return nnab_i

    def build_nblist(self,particles,bonds):
        """
        Create neighbor list of bonded particles

        Arguments:
            particles (dict) of particles with int keys
            bonds (dict) of bonds with int keys of particles 

        """
        self.list = []
        self.index = []
        self.cnt = -1 
        
        # Create 2D list of lists for each particle 
        nd2D = [ [] for pkey_i  in particles.keys() ]
        # Fill each particle list with it's neighbors based on the bonds
        for bkey_i, bond_i  in bonds.iteritems():            
            nd2D[bond_i.pkey2].append( bond_i.pkey1 )
            nd2D[bond_i.pkey1].append( bond_i.pkey2 )
        # Loop over all particles and add it's neighbors to 1D list  (NBlist.list)          
        #   while tracking the index in the 1D in the index list  (NBlist.index)
        for pkey_i  in particles.keys():
            self.index.append(self.cnt + 1)
            for pkey_j in nd2D[pkey_i]:
                if( pkey_i != pkey_j):
                    self.cnt += 1
                    self.list.append(pkey_j)
        # Add extra index positions for key+1 call made by final key 
        self.index.append(self.cnt + 1)
        # Clear 2D list from memory 
        del nd2D

    def guess_nblist(self,lat,particles,positions,radii_key,radii_buffer=1.25):
        """
        Create neighbor list of particles based on distance and radius of each particle 
        
        Arguments:
            lat (Lattice) object 
            particles (dict) of particles with int keys
            positions (list) of  particle position numpy arrays
            radii_key (str) of key for particle radius in the particles dict
            radii_buffer (float) to multiply radii cut off 

        """
        self.list = []
        self.index = []
        self.cnt = -1 
        # Create 2D list of lists of inter particle distances
        npos_i = positions
        npos_j = positions
        dr_matrix, dist_matrix  = lat.delta_npos(npos_i,npos_j)
        # Loop over all particles
        for pkey_i,particle_i  in particles.iteritems():
            radii_i = particle_i.properties[radii_key]
            self.index.append(self.cnt + 1)
            for pkey_j,particle_j in particles.iteritems():
                if( pkey_i != pkey_j):
                    radii_j = particle_j.properties[radii_key]
                    dr_cut = radii_i + radii_j
                    dr_cut = dr_cut*radii_buffer
                    logger.debug("Particles  i_%d - j_%d dr %f cut %f "%(pkey_i,pkey_j,dist_matrix[pkey_i,pkey_j],dr_cut))
                    if( dist_matrix[pkey_i,pkey_j] <= dr_cut ):
                        self.cnt += 1
                        self.list.append(pkey_j)
                    
        # Add extra index positions for key+1 call made by final key 
        self.index.append(self.cnt + 1)
        # Clear list from memory 
        del dr_matrix
        del dist_matrix


    def radii_nblist(self,lat,positions,radii,radii_buffer=1.25,write_dr=True,del_drmatrix=False):
        """
        Create neighbor list of particles based on distance and radius of each particle 
        
        Arguments:
            lat (Lattice) object 
            positions (list) of  particle position numpy arrays
            radii (list) of  radius of particle
            radii_buffer (float) to multiply radii cut off 

        """
        self.index = []
        self.list = []
        self.cnt = -1

        # Record rd
        if( write_dr ):
            dr_file = 'dr.csv'
            fout = open(dr_file,'wb')
            pair_writer = csv.writer(fout,delimiter=',')
            header = ['key_i','key_j','dr']
            #if( rank == 0 ):
            pair_writer.writerow(header)            
        
        # Create 2D list of lists of inter particle distances
        npos_i = positions
        npos_j = positions
        self.dr_matrix, self.dist_matrix  = lat.delta_npos(npos_i,npos_j)
        # Loop over all particles
        for key_i  in range(len(npos_i)):
            self.index.append(self.cnt + 1)
            radi_i = radii[key_i]
            for key_j  in range(len(npos_j)):
                radi_j = radii[key_j]
                if( key_i != key_j):
                    dr_cut = radi_i + radi_j
                    dr_cut = dr_cut*radii_buffer
                    dr = self.dist_matrix[key_i,key_j] 
                    if( dr <= dr_cut ):
                        self.cnt += 1
                        self.list.append(key_j)
                        if( write_dr  ):
                            row_i = [key_i,key_j,dr]
                            pair_writer.writerow(row_i)


        # Record rd
        if( write_dr ):
            fout.close()
                        
        # Add extra index positions for key+1 call made by final key 
        self.index.append(self.cnt + 1)
        if( del_drmatrix ):
            # Clear list from memory
            del self.dr_matrix
            del self.dist_matrix
        else:
            logger.debug(" Saving dr and dist matrix for groupset ")
        

        
class Replication():
    '''
    Object to recode the replication of structure 
    '''
    def __init__(self,name_i,name_j,name_ij,method,n):

        self.name_i = name_i
        self.name_j = name_j
        self.name_ij = name_ij
        self.method = method
        self.n = n
    
    def __del__(self):
        """
        'Magic' method for deleting contents of container
        """
        del self.name_ij
        del self.name_i
        del self.name_j
        del self.method
        del self.n
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s +  %s x %d ( %s ) -> %s "%(self.name_i,self.name_j,self.n,self.method,self.name_ij)

class GroupSet():
    """
    Set of groups within a structureContainer
    """
    def __init__(self,tag):
        """
        Constructor
        """
        # Set pointer for easy reference
        self.tag = tag 
        self.groups = dict()
        self.group_nblist = NBlist()        # Creates nblist object for other groups
        self.properties = dict()   # Dictionary of group properties

    def __del__(self):
        """
        Destructor
        """
        del self.tag 
        del self.groups
        del self.group_nblist
        del self.properties
        
    def calc_cent_mass(self):
        '''
        Calculate center of mass of groups and store them in list
        '''
        self.properties['cent_mass'] = []
        for gkey,group_i in self.groups.iteritems():
            group_i.centerofmass()
            self.properties['cent_mass'].append(group_i.properties['cent_mass'])

    def calc_dl(self):
        '''
        Calculate the maximum end to end distance of each group 
        '''
        self.properties['dl_sq'] = []
        for gkey,group_i in self.groups.iteritems():
            group_i.calc_dl()
            self.properties['dl_sq'].append(group_i.properties['dl_sq'])

    def calc_radius(self):
        '''
        Calculate radius of groups and store them in list
        '''
        self.properties['radius'] = []
        self.properties['r_gy_sq'] = []
        self.properties['Q_mn'] = []
        self.properties['Rgy_eignval'] = []
        self.properties['A_sphere'] = []
        self.properties['A_sphere_num'] = []
        self.properties['A_sphere_dem'] = []
        for gkey,group_i in self.groups.iteritems():
            group_i.calc_radius()
            group_i.calc_asphericity()
            self.properties['radius'].append(group_i.properties['radius'])
            self.properties['r_gy_sq'].append(group_i.properties['r_gy_sq'])
            self.properties['Q_mn'].append(group_i.properties['Q_mn'])
            self.properties['Rgy_eignval'].append(group_i.properties['Rgy_eignval'])
            self.properties['A_sphere'].append(group_i.properties['A_sphere'])
            self.properties['A_sphere_num'].append(group_i.properties['A_sphere_num'])
            self.properties['A_sphere_dem'].append(group_i.properties['A_sphere_dem'])

    def calc_dl(self):
        '''
        Calculate radius of groups and store them in list
        '''
        self.properties['dl_sq'] = []
        for gkey,group_i in self.groups.iteritems():
            group_i.calc_dl()
            self.properties['dl_sq'].append(group_i.properties['dl_sq'])
            
    def write_cm_xyz(self,group_file=""):
        '''
        Write center of mass of each group into an xyz file
        '''
        if( len(group_file) == 0 ):
            group_file = "%s_cm.xyz"%(self.tag)
        group_out = open(group_file,"w")
        group_line = " %d \n"%(len(self.groups.keys()))
        group_line += " \n"
        group_out.write(group_line)        
        for gkey,group_i in self.groups.iteritems():
            cent_mass_i = group_i.properties['cent_mass']
            group_line = " Ar   {} {} {}  \n".format(cent_mass_i[0],cent_mass_i[1],cent_mass_i[2])
            group_out.write(group_line)
        group_out.close()

    def write_xyzs(self):
        '''
        Write group coordinates into an xyz file
        '''
        for gkey,group_i in self.groups.iteritems():
            group_i.write_xyz()

    def group_pbcs(self):
        """
        Apply PBC's to create whole groups 

        Assumes group length is shorter than box length

        """
        for gkey,group_i in self.groups.iteritems():
            # Make sure group has particles 
            if( len( group_i.pkeys ) > 0 ):
                strucC = group_i.strucC
                # Get position of first particle in molecule
                pid_o  = group_i.pkeys[0]
                r_o = strucC.positions[pid_o]

                part_shifted = [False]*strucC.n_particles 

                r_mol_mass = np.zeros(strucC.lat.n_dim)
                shift = np.zeros(strucC.lat.n_dim)
                total_mass = 0.0 

                # shift all atoms to be conected 

                for pid_i in sorted(group_i.pkeys):
                    particle_i = strucC.particles[pid_i]
                    a_mass_i = particle_i.properties['mass']
                    r_i = strucC.positions[pid_i]
                    r_io = strucC.lat.deltasq_pos(r_i,r_o)
                    # sum center of mass
                    total_mass += a_mass_i
                    
                    shifted = False 
                    for dim in range(strucC.lat.n_dim):
                        shift_dim = round( r_io[dim]/  strucC.lat._matrix[dim][dim] )
                        r_i[dim] = r_i[dim]  + strucC.lat._matrix[dim][dim] * shift_dim
                        if( shift_dim != 0 ):
                            shifted = True 
                        r_mol_mass[dim] = r_mol_mass[dim]  + a_mass_i*r_i[dim] 

                    group_i.strucC.positions[pid_i] = r_i
                    r_o = r_i
                    pid_o = pid_i

                # Shift molecular center of mass into box 
                for dim in range(strucC.lat.n_dim):
                    cent_mass_i = r_mol_mass[dim] /total_mass
                    shift[dim] = strucC.lat._matrix[dim][dim] * round( cent_mass_i /  strucC.lat._matrix[dim][dim] )


                for pid_i in sorted(group_i.pkeys):
                    r_i = strucC.positions[pid_i]
                    for dim in range(strucC.lat.n_dim):
                        r_i[dim] = r_i[dim] - shift[dim] 

    def dump_json(self):
        '''
        Write group coordinates into an xyz file
        '''
        json_data = dict()
                
        for gkey,group_i in self.groups.iteritems():
            json_data[gkey] = group_i.pkeys

        f = open("groupset_%s.json"%(self.tag), 'w')
        json.dump(json_data,f, indent=2)
        f.close()
               
class Group():
    """
    Sets of particles within a structureContainer
    """

    def __init__(self,strucC,verbose=False):
        """
        Constructor
        """
        # Set pointer for easy reference 
        self.strucC = strucC
        self.n_dim = strucC.lat.n_dim
        #self.pid_list  = []
        self.pkeys  = []
        self.gkey = int(0)                   # group number/key
        self.bonded_nblist = NBlist()       # Creates nblist object for  bonded particles
        self.nonbonded_nblist = NBlist()    # Creates nblist object for nonbonded particles
        # 
        # Group properties
        # 
        self.properties = dict()
        
        '''
        self.cent_mass = np.zeros( self.n_dim)
        self.cent_charge = np.zeros( self.n_dim)
        self.total_mass = 0.0
        self.charge = 0.0

        
        self.radius = 0.0
        self.dl = 0.0
        self.dl_sq = 0.0 #  np.zeros( self.n_dim)
        self.r_gy = 0.0
        self.r_gy_sq = 0.0
        self.S_mn = np.zeros([self.n_dim,self.n_dim])
        self.Rgy_eignval = np.zeros(self.n_dim)
        self.asphere = 0.0
        '''

    def __del__(self):
        """
        Destructor
        """

        del self.n_dim 
        del self.pkeys
        del self.gkey
        del self.bonded_nblist
        del self.nonbonded_nblist
        del self.properties
        
        '''
        del self.cent_mass
        
        del self.cent_charge
        del self.total_mass
        del self.charge
        del self.radius
        del self.dl
        del self.r_gy
        del self.r_gy_sq
        del self.S_mn
        del self.Rgy_eignval
        del self.asphere
        '''

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
        F.write(" %d \n" % len(self.pkeys) )
        F.write("group  %s \n"%(self.tag))
        for pkey_i in self.pkeys:
            particle_i = self.strucC.particles[pkey_i]
            pos_i = self.strucC.positions[pkey_i]
            F.write( " %5s %16.8f %16.8f %16.8f \n"  % (particle_i.tag ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]) ) )
        F.close()

        return

  
    def centerofmass(self ):
        """
        Calculate the center of mass

        Center of mass of the molecule 

                   \sum_i r_i*mass_i 
        r_cmas =  __________________
                    \sum_i mass_i

        where r_i is the position of particle i and mass_i is it's mass 

        """
        #
        # Intialize center of mass list 
        self.properties['cent_mass'] = np.zeros( self.n_dim)
        self.properties['total_mass'] = 0.0
        #
        for pkey_i in self.pkeys:
            particle_i = self.strucC.particles[pkey_i]
            r_i = self.strucC.positions[pkey_i]
            mass_i = particle_i.properties["mass"]
            self.properties['total_mass'] += mass_i
            for dim in range(self.n_dim):
                self.properties['cent_mass'][dim] += mass_i*r_i[dim]
        # Normalize
        for dim in range(self.n_dim):
            self.properties['cent_mass'][dim] = self.properties['cent_mass'][dim]/self.properties['total_mass']
        #
        return

    def calc_radius(self):
        """
        Calculate the maximum radius

        and


        Radius of gyration^2 (r_gy_sq)

                  \sum_i  (r_i - r_cmas )^2
        r_gy_sq = __________________________
                      \sum_i


        Gyration tensor

                    \sum_i  (r_i_m - r_cmas_m ) (r_i_n - r_cmas_n )
        Q(m,n) = __________________________
                      \sum_i

        where m amd n are component of the r vector

        https://arxiv.org/pdf/1009.4358.pdf
        
        """
        # Intialize sums 
        self.properties['radius'] = 0.0
        self.properties['r_gy_sq'] = 0.0
        self.properties['Q_mn'] = np.zeros([self.n_dim,self.n_dim])
        
        
        for pkey_i in self.pkeys:
            particle_i = self.strucC.particles[pkey_i]
            r_i = self.strucC.positions[pkey_i]
            dr_cmass = self.strucC.lat.deltasq_pos_c(r_i,self.properties['cent_mass'])
            dot_dr_ij = dr_cmass.dot(dr_cmass)
            if( dot_dr_ij > self.properties['radius']  ):
                self.properties['radius'] = dot_dr_ij
            # \sum_i  (r_i_n - r_cmas_n ) (r_i_m - r_cmas_m )
            for d_m in range(self.n_dim):
                for d_n in range(self.n_dim):
                    self.properties['Q_mn'][d_m,d_n] += dr_cmass[d_m]*dr_cmass[d_n]
        # Normalize
        self.properties['radius'] = np.sqrt(self.properties['radius'])
        for d_m in range(self.n_dim):
            for d_n in range(self.n_dim):
                self.properties['Q_mn'][d_m,d_n]  = self.properties['Q_mn'][d_m,d_n]  /float( len(self.pkeys))
                
        for d_m in range(self.n_dim):
             self.properties['r_gy_sq'] += self.properties['Q_mn'][d_m,d_m]
        return

    def calc_asphericity(self):
        """
        Calculate the eigen values of the Gyration tensor and the Asphericity

         Soft Matter, 2013, 9, 3976-3984
         
        """

        eign_raw = np.linalg.eigvals(self.properties['Q_mn'])
        # Sort eignvalues 
        eign_raw.sort()
        Rgy_eignval = [x for x in reversed(eign_raw)]
        num = (Rgy_eignval[0] - Rgy_eignval[2])**2 + (Rgy_eignval[1] - Rgy_eignval[2])**2 + (Rgy_eignval[0] - Rgy_eignval[1])**2
        dem = (Rgy_eignval[0] + Rgy_eignval[1] + Rgy_eignval[2])**2
        # 
        # Add properties to property dictionary
        # 
        self.properties['Rgy_eignval'] = Rgy_eignval
        self.properties['A_sphere_num'] = num
        self.properties['A_sphere_dem'] = dem
        if( dem != 0.0 ):
            self.properties['A_sphere'] = num/(2.0*dem)
        else:
            self.properties['A_sphere'] = 0.0 

    def calc_dl(self ):
        """
        Calculate the maximum end to end distance 

        """

        # Intialize center of mass list 
        self.properties['dl_sq'] = 0.0 #  np.zeros( self.n_dim)
        
        for pkey_i in self.pkeys:
            particle_i = self.strucC.particles[pkey_i]
            r_i = self.strucC.positions[pkey_i]
            for pkey_j in self.pkeys:
                particle_j = self.strucC.particles[pkey_j]
                r_j = self.strucC.positions[pkey_j]
                dr_ij = self.strucC.lat.deltasq_pos_c(r_i,r_j)
                dot_dr_ij = dr_ij.dot(dr_ij)
                if( dot_dr_ij > self.properties['dl_sq'] ): self.properties['dl_sq'] = dot_dr_ij
                          
        

    def hterm_group(self,debug=False):
        """
        Hydrogen terminate group  

        (j)    (l)
            \  /
            (i)    - 109.5 deg
            /  \
          (k)   (m)
          
        Hydrogens will be added to sp3 carbons in a tetrahedral configuration

            ^         (l)
            |        /  |
       A    |   C  /    | 
            |    /      |
            |  /        |
            (i) -------->  
                B

        If atom i is bonded to two other atoms j and k
        r_ij = position_i - position_j
        r_ik = position_i - position_k
            A = r_ij x r_ik
            B = r_ij + r_ik

        For the angle (l)(i)(m) to be 109.5

        |A| = sin(  109.5/2 ) |C|
        |B| = cos(  109.5/2 ) |C|

        By scaling A and B they can be added to get the bonds between i and l and m

        r_il = A + B
        r_im = -A + B

        """


        latticevec = self.strucC.lat._matrix


        hb_length = 1.09
        hb_angle =  109.5
        tetrahedral_angle = np.deg2rad(hb_angle/2.0 )
        scale_vcross = np.sin(tetrahedral_angle) * hb_length
        scale_vadd = np.cos(tetrahedral_angle) * hb_length

        pt_H = Atom('H')
        pt_H.properties = periodictable.element_number(1)
        pt_H.properties["fftype"] = "HC"
        pt_H.properties["resname"] = "TERM"
        pt_H.properties["mol"] = self.properties["mol"] 
        pt_H.properties["residue"] = self.properties["residue"] 

        original_ref_mod = []
        sub_ref_mod = []
        pid_i_mod = 0

        Htermed = self.strucC.getSubStructure(self.pkeys,tag='%s_hterm'%(self.tag))

        for pkey_i in Htermed.particles.keys():
            NNAB_i = Htermed.bonded_nblist.calc_nnab(pkey_i)  # Group neighbors
            pkey_o = self.pkeys[pkey_i]
            NNAB_o = self.strucC.bonded_nblist.calc_nnab(pkey_o)# Structure container neighbors 
            if( NNAB_o != NNAB_i):
                # If there has been neighbors removed terminate 
                dB = NNAB_o - NNAB_i
                particle_i = Htermed.particles[pkey_i]
                r_i = Htermed.positions[pkey_i]
                r_ij_array = []
                for pkey_j in Htermed.bonded_nblist.getnbs(pkey_i):
                    r_j =  Htermed.positions[pkey_j]
                    r_ij = Htermed.lat.deltasq_pos(r_i,r_j)
                    r_ij_array.append(r_ij)

                    # print r_i,r_j,r_ij
                
                if( len(r_ij_array) != NNAB_i ):
                    error_line = " len(r_ij_array) {} != NNAB_i {} in groups.hterm() ".format(len(r_ij_array),NNAB_i)
                    sys.exit(error_line)

                if( NNAB_o == 3 and dB == 1 ):
                    logger.debug("Adding hterm_conjugated sp2")
                    pos_j = hterm_Csp2(hb_length,r_i,r_ij_array)
                    pt_H.properties["fftype"] = "HA"
                    Htermed.add_partpos(pt_H,pos_j,deepcopy = True)
                    p_j = Htermed.n_particles -1

                    Bond_iH = Bond(pkey_i,p_j)
                    Htermed.add_bond(Bond_iH)

                elif( NNAB_o == 4 and dB == 1 ):
                    logger.debug("Adding hterm_sp3 ")
                    pos_j = hterm_Csp3(hb_length,r_i,r_ij_array)
                    pt_H.properties["fftype"] = "HC"
                    Htermed.add_partpos(pt_H,pos_j,deepcopy = True)
                    p_j = Htermed.n_particles -1 
                    Bond_iH = Bond(pkey_i,p_j)
                    Htermed.add_bond(Bond_iH)


                elif( NNAB_o == 4 and dB == 2 ):

                    logger.debug("Adding 2 H to hterm_sp3 ")
                    #r_new = hterm_sp3v2(hb_length,r_i,r_ij_array)

                    cros_jk = np.cross(r_ij_array[0],r_ij_array[1])
                    add_jk = -1.0*( r_ij_array[0] + r_ij_array[1] )
                    cros_scale = scale_vcross*cros_jk/np.linalg.norm(cros_jk)
                    add_scale = scale_vadd*add_jk/np.linalg.norm(add_jk)

                    hbond_1 = cros_scale + add_scale
                    hbond_2 =  -1.0*cros_scale + add_scale

                    hpos_1 = r_i + hbond_1
                    hpos_2 = r_i + hbond_2

                    pt_H.properties["fftype"] = "HC"
                    Htermed.add_partpos(pt_H,hpos_1,deepcopy = True)
                    p_j = Htermed.n_particles -1 
                    Bond_iH = Bond(pkey_i,p_j)
                    Htermed.add_bond(Bond_iH)

                    pt_H.properties["fftype"] = "HC"
                    Htermed.add_partpos(pt_H,hpos_2,deepcopy = True)
                    p_j = Htermed.n_particles -1 
                    Bond_iH = Bond(pkey_i,p_j)
                    Htermed.add_bond(Bond_iH)

                    #original_ref_mod.append(-1)
                    #pid_i_mod += 1 
                    #sub_ref_mod.append(pid_i_mod)

                    #  Add H-(i)-H angle
                    pid_H_j = Htermed.n_particles
                    a_i = Angle( pid_H_j-1 ,pkey_i, pid_H_j )            
                    Htermed.add_angle(a_i)

                    #original_ref_mod.append(-1)
                    #pid_i_mod += 1 
                    #sub_ref_mod.append(pid_i_mod)
                    
                else:
                    error_line =  " Nubmer of missing atoms %d has yet to be accounted for in groups.hterm \n"%(dB)
                    error_line +=  " {} -> {}".format(NNAB_o,NNAB_i)
                    Htermed.write_xyz("hterm_failed.xyz")
                    sys.exit(error_line)
                # Redo neighbor list
                Htermed.bonded_nblist.build_nblist(Htermed.particles,Htermed.bonds )
            
        return Htermed #original_ref_mod,sub_ref_mod 

                 
class Container():
    """
    Data structure for describing a collection of Particles that have associated
    positions within a Lattice, and consistent set keys corresponding to
    Bond, Angle, Dihedral and Improper descriptions 
    """

    def __init__(self,tag="blank",verbose=False):
        """
        Constructor for a composite structure. 
        """
        if isinstance(tag, str):
            self.tag = tag
        else:
            raise TypeError("1st arg (tag) in %s Container initialization should be string"%(__name__))
        
        self.lat = Lattice()                                  # Creates lattice object for structure
        self.bonded_nblist = NBlist()                         # Creates nblist object for  bonded particles
        self.nonbonded_nblist = NBlist()                      # Creates nblist object for nonbonded particles
        self.particles = dict()                               # Creates empty dict struc
        self.prop_particles =  dict() 
        self.positions = []                                   # Creates empty array
        self.bonds = dict()                                   # Creates empty dict struc
        self.angles = dict()                                  # Creates empty dict struc
        self.dihedrals = dict()                                # Creates empty dict struc
        self.impropers = dict()                                # Creates empty dict struc
        # Int count of the length of each dictionary
        #   mostly for internal use 
        self.n_particles = 0    
        self.n_bonds = 0    
        self.n_angles = 0    
        self.n_dihedrals = 0    
        self.n_impropers = 0    

        # Container properties
        self.properties = dict()
        
        self.properties['mass'] = 0.0 
        self.properties['volume'] = 0.0 
        self.properties['density'] = 0.0 
        self.properties['center_mass'] = np.zeros(self.lat.n_dim)
        self.properties['dipole'] = np.zeros(self.lat.n_dim)
        # Reference information 
        self.properties['name'] = ""   # Tag of structure to be set by file read in 
        self.properties['chemicalformula'] = ""
        self.properties['IUPAC'] = ""
        self.properties['common_tag'] = ""
        self.properties['deptag'] = ""
        self.properties['ctag'] = ""
        self.properties['moltype'] = ""
        self.properties['backbone'] = ""
        self.properties['composition'] = ""
        # Particles with dangling bond
        self.properties['danglkey'] = -1
        
        # Groups within structure 
        self.groupsets = dict()
        # Track replications 
        self.replications = []
        
    def __del__(self, verbose=False):
        """
        Constructor for a composite structure. 
        """
        del self.lat 
        del self.particles
        del self.positions
        del self.bonds
        del self.angles
        del self.dihedrals
        del self.impropers
        # Del counts 
        del self.n_particles
        del self.n_bonds
        del self.n_angles
        del self.n_dihedrals
        del self.n_impropers
        # Del properties 
        del self.properties
        del self.groupsets
        del self.replications
        
    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s "%(self.tag)

    def dump_pickle(self):
        '''
        Pickle object
        '''
        file_i = open("%s.pkl"%(self.tag),'w')
        pickle.dump(self,file_i)
        file_i.flush()


    def add_particle(self, particle_i, deepcopy = False ):
        """
        Add 'Particle' object to this container and update n_particles accordingly 
        """
        if( deepcopy ):
                self.particles[self.n_particles] = copy.deepcopy(particle_i) # index 0 -> (N-1)
        else:
                self.particles[self.n_particles] = particle_i # index 0 -> (N-1)
                
        self.n_particles = len(self.particles)

    def add_position(self, pos_i):
        """
        Append 'position' as numpy array to this container. 
        """
        if( len(pos_i) == self.lat.n_dim ):
            self.positions.append(np.array(pos_i))
        else:
            print "Attempting to add non-%d-dimension position to container"%(n_dim)
            raise TypeError

    def add_partpos(self, particle_i, pos_i, deepcopy = True):
        """
        Add 'Particle' object to this container and update n_particles accordingly and
         append 'position' array to this container. 
        """
        self.add_particle( particle_i, deepcopy = deepcopy)
        self.add_position( pos_i)


    def add_bond(self, bond_i, deepcopy = True ):
        """
        Add 'Bond' object to bonds dict in this container and update n_bonds accordingly
        """
        if isinstance(bond_i, Bond):
            self.n_bonds = len(self.bonds)
            if( deepcopy ):
                self.bonds[self.n_bonds] = copy.deepcopy(bond_i) # index 0 -> (N-1)
            else:
                self.bonds[self.n_bonds] = bond_i # index 0 -> (N-1)
                
            self.n_bonds = len(self.bonds)
        else:
            print "Attempting to add non-Bond type to container"
            raise TypeError


    def add_angle(self, angle_i, deepcopy = True ):
        """
        Add 'Angle' object to angles dict in this container and update n_angles accordingly
        """
        if isinstance(angle_i, Angle):
            self.n_angles = len(self.angles)
            if( deepcopy ):
                self.angles[self.n_angles] = copy.deepcopy(angle_i) # index 0 -> (N-1)
            else:
                self.angles[self.n_angles] = angle_i # index 0 -> (N-1)
                
            self.n_angles = len(self.angles)
        else:
            print "Attempting to add non-Angle type to container"
            raise TypeError


    def add_dihedral(self, dihedral_i, deepcopy = True ):
        """
        Add 'Dihedral' object to dihedrals dict in this container and update n_dihedrals accordingly
        """
        if isinstance(dihedral_i, Dihedral):
            self.n_dihedrals = len(self.dihedrals)
            if( deepcopy ):
                self.dihedrals[self.n_dihedrals] = copy.deepcopy(dihedral_i) # index 0 -> (N-1)
            else:
                self.dihedrals[self.n_dihedrals] = dihedral_i # index 0 -> (N-1)
                
            self.n_dihedrals = len(self.dihedrals)
        else:
            print "Attempting to add non-Dihedral type to container"
            raise TypeError


    def add_improper(self, improper_i, deepcopy = True ):
        """
        Add 'Improper' object to impropers dict in this container and update n_impropers accordingly
        """
        if isinstance(improper_i, Improper):
            self.n_impropers = len(self.impropers)
            if( deepcopy ):
                self.impropers[self.n_impropers] = copy.deepcopy(improper_i) # index 0 -> (N-1)
            else:
                self.impropers[self.n_impropers] = improper_i # index 0 -> (N-1)
                
            self.n_impropers = len(self.impropers)
        else:
            print "Attempting to add non-Improper type to container"
            raise TypeError

    def write_coord(self):
        """
        Write a coordnates into string 
        """
        coord = ''.join([" %5s %16.8f %16.8f %16.8f \n"%(particle_i.tag,self.positions[pkey_i][0],self.positions[pkey_i][1],self.positions[pkey_i][2] ) for pkey_i,particle_i in self.particles.iteritems()])

        return coord
    
    def write_xyz_str(self):
        '''
        Write xyz file string
        '''

        xyz_str = " %d \n" % self.n_particles
        xyz_str += " %s \n"%(self.tag)
        xyz_str += self.write_coord()
        return str(xyz_str)
    
    def write_xyz(self, xyz_file=''):
        '''
        Write a structure  to an xyz file

        Args:
            xyz_file (str) xyz file tag
            
        Reutrns:
            null
            
        '''
        if( len(xyz_file) == 0 ):
            xyz_file = "%s.xyz"%(self.tag)
            
        xyz_str = self.write_xyz_str()
        
        F = open(xyz_file,"w")
        F.write(xyz_str)
        F.close()


    def write_xyz_list(self, list_i,xyz_file=''):
        '''
        Write a list of certain particles of the structure  to an xyz file

        Args:
            xyz_file    (str) xyz file tag
            
        Reutrns:
            null
            
        '''
        if( len(xyz_file) == 0 ):
            xyz_file = "%s.xyz"%(self.tag)

        F = open(xyz_file,"w")
        
        # Loop over structures
        F.write(" %d \n" % len(list_i) )
        F.write(" %s \n"%" structure.Container  ")
        for pkey_i  in list_i:
            particle_i = self.particles[pkey_i]
            pos_i = self.positions[pkey_i]
            F.write(" %5s %16.8f %16.8f %16.8f \n"%(particle_i.tag,pos_i[0],pos_i[1],pos_i[2] ))
        F.close()

    def read_xyz(self, xyz_file=''):
        '''
        Read a structure  to an xmol file
        '''
        if( len(xyz_file) == 0 ):
            xyz_file = "%s.xyz"%(self.tag)
            
        line_cnt = 0
        try:
            
            with open(xyz_file) as f:
                for line in f:

                    line_cnt += 1
                    col = line.split()
                    if( line_cnt > 2 and len(col) >= 4 ):
                        # Read lines and add particles to structure
                        symbol = str(col[0])
                        pos_i = np.array( [float(col[1]),float(col[2]),float(col[3])] )
                        pt_i = Atom(symbol)
                        self.add_partpos(pt_i,pos_i,deepcopy = True)
        except:
            logger.warning(" File not found %s in %s "%(xyz_file,os.getcwd()))
        
        return
    
    def write_list(self,list_i,tag_i):
        '''
        Write list of particle keys to  file to use in remote analysis
        '''
        list_str = [str(pkey) for pkey in list_i]
        list_file = '%s.list'%(tag_i)
        outfile = open(list_file,'wb')
        outfile.write("\n".join(list_str))
        outfile.close()        

        return list_file

    def shift(self, pkey, vec):
        """
        Shift position of pkey by vector

        Arguments:
            pkey (int) particle key
            vec  (np.array) vector 
        """
        self.positions[pkey] += vec

    def shift_pos(self,vec):
        '''
        Shift position of all particles by vecx
        
        '''
        for pkey_i in self.particles.keys():
            self.shift( pkey_i, vec)


    def pbc_pos(self):
        '''
        Apply periodic boundry conditions to 
        '''
        for r_i in self.positions:
            for d in range(self.lat.n_dim ):
                r_i[d] = r_i[d] - self.lat._matrix[d][d] * round( r_i[d]/  self.lat._matrix[d][d] )
                    
    def lat_cubic(self,len_o):
        '''
        Set lattice to cubic with lattice constant len
        '''
        matrix = np.zeros((self.lat.n_dim, self.lat.n_dim))
        for d in range(self.lat.n_dim):
            matrix[d][d] = len_o
        self.lat.set_matrix(matrix)
                
    def calc_mass(self):
        """
        Calculate total mass of structure  
        """
        self.mass = float(0.0)
        
        for pkey_i, particle_i  in self.particles.iteritems():
            self.mass += particle_i.properties["mass"]

        return

    def calc_charge(self):
        """
        Calculate total charge of structure  
        """
        self.charge = float(0.0)
        
        for pkey_i, particle_i  in self.particles.iteritems():
            self.charge += particle_i.properties["charge"]

        return


    def calc_volume(self):
        """
        Calculate volume of structure  
        Volume = ( v_i x v_j ) . v_k 
        """
        v_i = self.lat._matrix[0] 
        v_j = self.lat._matrix[1] 
        v_k = self.lat._matrix[2]
        
        v_ij = np.cross(v_i,v_j)
        self.volume = np.dot(v_ij,v_k)
        
        return 

    def calc_density(self):
        """
        Calculate density of structure  
        """
        self.density = self.mass/self.volume

                

    def calc_composition(self):
        """
        Calculate composition
        """
        # Find max mol and residue numbers
        self.mol_max = -1
        self.residue_max = -1
        for pkey_i, particle_i  in self.particles.iteritems():
            if( particle_i.properties["mol"] > self.mol_max ): self.mol_max =  particle_i.properties["mol"]
            if( particle_i.properties["residue"] > self.residue_max ): self.residue_max =  particle_i.properties["residue"]


        self.composition = np.zeros(periodictable.n_elements,dtype=np.int)    
        for pkey_i, particle_i  in self.particles.iteritems():
            el_i = int( particle_i.properties["number"] )
            if( el_i >= 0 ):
                self.composition[el_i] += 1
    
    def calc_formula(self):
        """
        Calculate chemical formula
        """

        self.chemicalformula = ""
        self.calc_composition()

        el_n_list = [6,1]
        el_n_list += [ i for i in range(1,periodictable.n_elements) if( i != 6 and i != 1 ) ]
        for n_i in el_n_list:
            if( self.composition[n_i] > 0 ):
                el_i = periodictable.element_number(n_i)
                self.chemicalformula += "%s%d"%(el_i["symbol"],self.composition[n_i])
                            
    def calc_center_mass(self):
        """
        Find center of mass of a structure
        """

        self.center_mass = np.zeros(self.lat.n_dim)

        for pkey_i, particle_i  in self.particles.iteritems():
            mass_i = particle_i.properties["mass"]
            # print self.positions[pkey_i][0],self.positions[pkey_i][1],self.positions[pkey_i][2],mass_i
            self.center_mass += mass_i*np.array(self.positions[pkey_i])

        self.center_mass = self.center_mass/self.mass

        return

    def sum_prop(self,propkey,pkey_i,pkeyp_j):
        '''
        Sum property of particle i into particle j
        '''
        try:
            # Sum charges of particles to be removed into attachment points
            #print " Summing ",self.particles[pkeyp_j].properties['symbol'],self.particles[pkeyp_j].properties['charge']
            #print " into ",self.particles[pkey_i].properties['symbol'],self.particles[pkey_i].properties['charge']
            self.particles[pkey_i].properties['charge'] += self.particles[pkeyp_j].properties['charge']
            self.particles[pkeyp_j].properties['charge'] = 0.0
        except:
            logger.warning("Container %s does not have particle charges"%(pkeyp_j.tag))
                                 
    def maxtags(self):
        """
        Find max mol and residue numbers
        """
        self.mol_max = -1
        self.residue_max = -1
        for pkey_i, particle_i  in self.particles.iteritems():
            if( particle_i.properties["mol"] > self.mol_max ): self.mol_max =  particle_i.properties["mol"]
            if( particle_i.properties["residue"] > self.residue_max ): self.residue_max =  particle_i.properties["residue"]

    def mol_mult(self):
        """
        Find value to multiply mol index by
        so reisude index can be added to get a unique group value
        """
        
        self.mol_multiplier =  float( len( str( abs( round(self.mol_max,0) )))*10.0/len( str( abs( round(self.mol_max,0) ))) )*10.0

    def group_prop(self,prop,tag,particles_select=[]):
        """
        Create groups of mols or residues

        Args:
            prop  (str) property key
            tag (str) key for Groups object in groups dict
            particles_select  (list) list of particle keys to be included in groups
                                     defaults to all particles.keys if none are specified
        Returns:
            None 
        """
        supported_group_props = ['mol','residue']
        if( prop not in supported_group_props):
            logger.warning(" Unsupported property selection for groups %s "%( prop))
            sys.exit(2)        

        if( len(particles_select) == 0 ):
            particles_select = self.particles.keys()
        # 
        self.maxtags()
        self.mol_mult()
        #
        # Group particles 
        #
        groupset_i = GroupSet(tag)
        groupset_i.properties['keys'] = []
        groupset_i.properties['tags'] = []
        group_uniqueid = dict()
        #gkey = 0 
        for pkey_i, particle_i  in self.particles.iteritems():
            if( pkey_i in particles_select ):
                # If in selected set of particles 
                uniqueid_i = pkey_i
                if( prop == "mol" ): 
                    uniqueid_i =  int(particle_i.properties["mol"] )
                elif( prop == "residue" ): 
                    uniqueid_i =  int( particle_i.properties["mol"]*self.mol_multiplier + particle_i.properties["residue"] )
                else:
                    logger.warning(" Unsupported property selection for groups %s "%( prop))
                    sys.exit(2)
                if( uniqueid_i not in group_uniqueid.keys() ):
                    gkey = len(groupset_i.groups)
                    group_uniqueid[uniqueid_i] = gkey
                    group_i = Group(self)
                    group_i.gkey = gkey
                    group_i.tag = "%s_%s"%(tag,gkey)
                    group_i.properties["mol"]  = particle_i.properties["mol"] 
                    group_i.properties["residue"] = particle_i.properties["residue"] 
                    group_i.properties["resname"] = particle_i.properties["resname"] 
                    groupset_i.groups[gkey] = group_i
                    # 
                    groupset_i.properties['keys'].append(group_i.gkey)
                    groupset_i.properties['tags'].append(group_i.tag)
                else:
                    gkey = group_uniqueid[uniqueid_i]
                    group_i =  groupset_i.groups[gkey]
                group_i.pkeys.append(pkey_i)
        # Store set of groups in dict
        self.groupsets[tag] = groupset_i
        
        return  


            

    def n_molecules(self):
        """
        Number of molecules
        """
        max_mol = 0
        for pkey_i, particle_i  in self.particles.iteritems():
            if( max_mol < particle_i.properties["mol"] ): max_mol = particle_i.properties["mol"]
        return max_mol

    def rotate_xz(self,theta_xz,direction="counterclockwise",verbose=False):
        """

        Rotate around the y-axis in the xz plane
        
             |   cos theta_xz 0 -sin theta_xz   |
        Ry = |           0    1      0         |
             |_  sin theta_xz 0  cos theta_xz _|

        Arguments:
            theta_xz (float) angle in radians
            direction (str)  counterclockwise or clockwise around y axis 

        """
        
        def Rxzdotv(v_i,cos_xz,sin_xz,sinprefix12,sinprefix21):
            verbose_line = " Rotating vector {} ".format(v_i)
            
            v_j = np.zeros(len(v_i))
            v_j[0] = cos_xz*v_i[0] + sinprefix12*sin_xz*v_i[2] 
            v_j[1] = v_i[1] 
            v_j[2] = sinprefix21*sin_xz*v_i[0] + cos_xz*v_i[2] 
            return v_j
            
        if( verbose ):
            print " Rotating particle {} around y-axis ".format(direction)
            
        
        if( self.n_particles > 0 ):
            cos_xz = math.cos(theta_xz)
            sin_xz = math.sin(theta_xz)
            if(verbose):
                print direction
                print "cos_xz ",cos_xz
                print "sin_xz ",sin_xz
                
            if( direction == "clockwise"  ):
                sinprefix12 = 1.0
                sinprefix21 = -1.0
            elif( direction == "counterclockwise"  ):
                sinprefix12 = -1.0
                sinprefix21 = 1.0
            else:
                error_line = "!!! Error in structurContainer.rotate_xz() !!! \n"
                error_line += "!!! Unknow direction selected {} please select counterclockwise or clockwise !!!".format(direction)
                sys.exit(error_line)
            
            for pkey_i  in self.particles.keys():
                self.positions[pkey_i] = Rxzdotv(self.positions[pkey_i],cos_xz,sin_xz,sinprefix12,sinprefix21)

        return
    
    def rotate_xy(self,theta_xy,direction="counterclockwise",verbose=False,debug = False):
        """

        Rotate around the z-axis in the xy plane
        
            z        v_i 
            |       /
            |      /
            |    /
            |  /
            |/_______________ y
             \
              \
               \
                \
                 \
                  x
              _                                _
             |   cos theta_xy  sin theta_xy  0  |
        Rz = |   -sin theta_xy  cos theta_xy  0  |
             |_         0           0        1 _|

        Arguments:
            theta_xy (float) angle in radians
            direction (str)  counterclockwise or clockwise around z axis 
        """

        
        def Rxydotv(v_i,cos_xy,sin_xy,sinprefix12,sinprefix21):
            verbose_line = " Rotating vector {} ".format(v_i)
            
            v_j = np.zeros(len(v_i))
            v_j[0] = cos_xy*v_i[0] + sinprefix12*sin_xy*v_i[1] 
            v_j[1] = sinprefix21*sin_xy*v_i[0] + cos_xy*v_i[1] 
            v_j[2] = v_i[2] 
            return v_j
            
        if( verbose ):
            print " Rotating particle {} around z-axis ".format(direction)

        
        if( debug ):

            v_i = [1.0,1.0,0.0]

            cos_xy = math.cos(np.pi/4.0)
            sin_xy = math.sin(np.pi/4.0)
            sinprefix12 = 1.0
            sinprefix21 = -1.0
            print "v_i",v_i
            print "cos_xy",cos_xy
            print "sin_xy",sin_xy
            print Rxydotv(v_i,cos_xy,sin_xy,sinprefix12,sinprefix21)
            sys.exit(" debug Rxydotv in rotate_xy")

        
        if( self.n_particles > 0 ):
            cos_xy = math.cos(theta_xy)
            sin_xy = math.sin(theta_xy)
            if( verbose ):
                print direction
                print "cos_xy ",cos_xy
                print "sin_xy ",sin_xy
                
            if( direction == "clockwise"  ):
                sinprefix12 = 1.0
                sinprefix21 = -1.0
            elif( direction == "counterclockwise"  ):
                sinprefix12 = -1.0
                sinprefix21 = 1.0
            else:
                error_line = "!!! Error in structurContainer.rotate_xy() !!! \n"
                error_line += "!!! Unknow direction selected {} please select counterclockwise or clockwise !!!".format(direction)
                sys.exit(error_line)

            for pkey_i  in self.particles.keys():
                self.positions[pkey_i] = Rxydotv(self.positions[pkey_i],cos_xy,sin_xy,sinprefix12,sinprefix21)
                    
        return

    def rotate_yz(self,theta_yz,direction="counterclockwise",verbose=False):
        """

        Rotate around the x-axis in the yz plane
        
            z        v_i 
            |       /
            |      /
            |    /
            |  /
            |/_______________ y
             \
              \
               \
                \
                 \
                  x
               _                                    _
              |      1        0           0          |
        Rx =  |      0  cos theta_xy  sin theta_xy   |
              |_     0  -sin theta_xy  cos theta_xy _|
              
        Arguments:
            theta_yz (float) angle in radians
            direction (str)  counterclockwise or clockwise around x axis 
        
        """

        
        def Ryzdotv(v_i,cos_yz,sin_yz,sinprefix12,sinprefix21):
            verbose_line = " Rotating vector {} ".format(v_i)
            
            v_j = np.zeros(len(v_i))
            v_j[0] = v_i[0] 
            v_j[1] = cos_yz*v_i[1] + sinprefix12*sin_yz*v_i[2] 
            v_j[2] = sinprefix21*sin_yz*v_i[1] + cos_yz*v_i[2] 
            return v_j
            
        if( verbose ):
            print " Rotating particle {} around x-axis ".format(direction)

        debug = False 
        if( debug ):

            v_i = [0.0,-1.0,1.0]

            cos_yz = math.cos(np.pi/4.0)
            sin_yz = math.sin(np.pi/4.0)
            sinprefix12 = 1.0
            sinprefix21 = -1.0
            print "v_i",v_i
            print "cos_yz",cos_yz
            print "sin_yz",sin_yz
            print Ryzdotv(v_i,cos_yz,sin_yz,sinprefix12,sinprefix21)
            sys.exit(" debug Ryzdotv in rotate_yz")

        
        if( self.n_particles > 0 ):
            cos_yz = math.cos(theta_yz)
            sin_yz = math.sin(theta_yz)
            if(verbose):
                print direction
                print "cos_yz ",cos_yz
                print "sin_yz ",sin_yz
                
            if( direction == "clockwise"  ):
                sinprefix12 = 1.0
                sinprefix21 = -1.0
            elif( direction == "counterclockwise"  ):
                sinprefix12 = -1.0
                sinprefix21 = 1.0
            else:
                error_line = "!!! Error in structurContainer.rotate_yz() !!! \n"
                error_line += "!!! Unknow direction selected {} please select counterclockwise or clockwise !!!".format(direction)
                sys.exit(error_line)
            
            for pkey_i  in self.particles.keys():
                self.positions[pkey_i] = Ryzdotv(self.positions[pkey_i],cos_yz,sin_yz,sinprefix12,sinprefix21)
                        
        return

    def del_particle(self,pkey):
        '''
        Remove a particle from Container 
        '''

        self.keyupdate = dict()  # Need to keep track of all ptcl key changes at once
        
        # Create copies to iterate over 
        particles_i = copy.deepcopy( self.particles )
        positions_i = copy.deepcopy( self.positions )
        bonds_i = copy.deepcopy(self.bonds)
        angles_i = copy.deepcopy(self.angles)
        dihedrals_i = copy.deepcopy(self.dihedrals)
        impropers_i = copy.deepcopy(self.impropers)
        # Remove particle from copy 
        del particles_i[pkey]
        # Re initialized container list's and counts 
        self.particles = dict()                               # Creates empty dict struc
        self.positions = []                                   # Creates empty array
        self.n_particles = 0    
        self.bonds = dict()                               # Creates empty dict struc
        self.angles = dict()                                  # Creates empty dict struc
        self.dihedrals = dict()                                # Creates empty dict struc
        self.impropers = dict()                                # Creates empty dict struc
        self.n_bonds = 0                                   # Creates empty array
        self.n_angles = 0    
        self.n_dihedrals = 0    
        self.n_impropers = 0    
        # Add particles to container and track new keys 
        toPtclID = 0 
        for  pkey_i in particles_i.keys():
            self.keyupdate[pkey_i] = toPtclID               # Store ID changes
            toPtclID += 1
            particle_i = particles_i[pkey_i]
            pos_i = positions_i[pkey_i]
            self.add_partpos(particle_i,pos_i)
            
        # Add bonds to container and use new keys 
        for bkey_i,bond_i in bonds_i.iteritems():
            if( bond_i.pkey1 != pkey  and  bond_i.pkey2 != pkey ):
                bond_i.pkey1 = self.keyupdate[bond_i.pkey1]
                bond_i.pkey2 = self.keyupdate[bond_i.pkey2]
                self.add_bond(bond_i)
        # Add angles to container and use new keys 
        for key_i,angle_i in angles_i.iteritems():
            if( angle_i.pkey1 != pkey  and  angle_i.pkey2 != pkey  and  angle_i.pkey3 != pkey ):
                angle_i.pkey1 = self.keyupdate[angle_i.pkey1]
                angle_i.pkey2 = self.keyupdate[angle_i.pkey2]
                angle_i.pkey3 = self.keyupdate[angle_i.pkey3]
                self.add_angle(angle_i)
                
        # Add dihedrals to container and use new keys 
        for key_i,dih_i in dihedrals_i.iteritems():
            if( dih_i.pkey1 != pkey  and  dih_i.pkey2 != pkey  and  dih_i.pkey3 != pkey  and  dih_i.pkey4 != pkey ):
                dih_i.pkey1 = self.keyupdate[dih_i.pkey1]
                dih_i.pkey2 = self.keyupdate[dih_i.pkey2]
                dih_i.pkey3 = self.keyupdate[dih_i.pkey3]
                dih_i.pkey4 = self.keyupdate[dih_i.pkey4]
                self.add_dihedral(dih_i)
        # Add impropers to container and use new keys 
        for key_i,imp_i in impropers_i.iteritems():
            if( imp_i.pkey1 != pkey  and  imp_i.pkey2 != pkey  and  imp_i.pkey3 != pkey  and  imp_i.pkey4 != pkey ):
                imp_i.pkey1 = self.keyupdate[imp_i.pkey1]
                imp_i.pkey2 = self.keyupdate[imp_i.pkey2]
                imp_i.pkey3 = self.keyupdate[imp_i.pkey3]
                imp_i.pkey4 = self.keyupdate[imp_i.pkey4]
                self.add_improper(imp_i)

        # Remake neighbor list based on updated bonds 
        self.bonded_nblist = NBlist() 
        self.bonded_nblist.build_nblist(self.particles,self.bonds )


    def get_max(self):
        """
        get max mol, residue and charge group number 
        """
        self.max_mol = 0 
        self.max_residue = 0 
        self.max_qgroup = 0 
        self.max_ring = 0
        for j,p_j in self.particles.iteritems():
            if( p_j.properties["mol"] > self.max_mol ): self.max_mol =   p_j.properties["mol"] 
            if(  p_j.properties["residue"] > self.max_residue ): self.max_residue =   p_j.properties["residue"] 
            if(  p_j.properties["qgroup"] > self.max_qgroup ): self.max_qgroup =   p_j.properties["qgroup"]
            if(  p_j.properties["ring"] > self.max_ring ): self.max_ring =   p_j.properties["ring"]


    def add_mol(self,max_ref_mol ):
        """
        Set minimum mol number to references
        """
        for pid, ptclObj in self.ptclC :
            ptclObj.properties["mol"] +=   max_ref_mol  - self.max_mol 

    def shift_tag(self,tag,tag_min):
        """
        shift tag by an number  
        """
        for pid, ptclObj  in self.particles.iteritems():
                # if( ptclObj.properties[tag] > 0 ):
                ptclObj.properties[tag] += tag_min

    def getSubStructure(self,pkeys,tag="blank"):
        """
        Create new structure container from list of particle keys
        """
        new_strucC = Container(str(tag))
        
        key_update = dict()
        # Set lattice 
        new_strucC.lat = self.lat
        # Set particles
        for pkey_i in pkeys:
            p_i = self.particles[pkey_i]
            pos_i = self.positions[pkey_i]            
            new_strucC.add_partpos(p_i,pos_i, deepcopy = True)
            key_update[pkey_i]  = new_strucC.n_particles -1
            
        if( len(self.bonded_nblist.index) > 0 ):
            # Update bonded nieghbor list
            new_strucC.bonded_nblist = NBlist() 
            for pkey_i in pkeys:
                new_strucC.bonded_nblist.index.append(new_strucC.bonded_nblist.cnt + 1)
                for pkey_j in self.bonded_nblist.getnbs(pkey_i):
                    if( pkey_j in pkeys ):
                        new_strucC.bonded_nblist.cnt += 1 
                        new_strucC.bonded_nblist.list.append(key_update[pkey_j])

            new_strucC.bonded_nblist.index.append(new_strucC.bonded_nblist.cnt + 1)
            new_strucC.bonded_bonds()

        return new_strucC
        

    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator (eg struc1 += struc2)
        
        Compare global IDs of particles and reassign globalIDs for particle
        container using the max ID between the two lists. Tracks these changes
        for all other (bond, angle, dihedral and improper ) containers that reference particleIDs
        """
        
        # Empty container checks
        if other.n_particles == 0:             # If struc2 is empty (has no particles)
            return self                 #   simply return unchanged current container
            
        new_strcC = copy.deepcopy(other) #    full copy other and return
        if self.n_particles == 0:              # If struc1 (this struc) is empty (has no particles)
            return new_strcC

        self.keyupdate = dict()  # Need to keep track of all ptcl key changes at once
        #
        # Add particles 
        for pkey_other,particle_other in new_strcC.particles.iteritems():
            fromPtclID = pkey_other                           # Track IDs from--->to
            toPtclID   = self.n_particles                  #  --> toID (to is the maxid of this ptclC)
            self.add_partpos(particle_other,new_strcC.positions[pkey_other]) # Pushes ptcl to this struc's ptcl container
            self.keyupdate[fromPtclID]=toPtclID               # Store ID changes

        # Add other neighbor lists to self with new particle keys
        #  Bonded
        if( new_strcC.bonded_nblist.cnt > 0 ):
            # Remove last place holder in nblist index to add other neighbors to be appended properly
            self.bonded_nblist.index.pop()
            for pkey_i, particle_i in new_strcC.particles.iteritems():
                self.bonded_nblist.index.append(self.bonded_nblist.cnt + 1)
                for pkey_j in   new_strcC.bonded_nblist.getnbs(pkey_i):
                    self.bonded_nblist.cnt += 1
                    self.bonded_nblist.list.append(self.keyupdate[pkey_j])
            self.bonded_nblist.index.append(self.bonded_nblist.cnt + 1)
        
        #  Non-Bonded 
        if( new_strcC.nonbonded_nblist.cnt > 0 ):
            # Remove last place holder in nblist index
            #  to add other neighbors to be appended properly
            self.nonbonded_nblist.index.pop()
            for pkey_i, particle_i in new_strcC.particles.iteritems():
                self.nonbonded_nblist.index.append(self.nonbonded_nblist.cnt + 1)
                for pkey_j in   new_strcC.nonbonded_nblist.getnbs(pkey_i):
                    self.nonbonded_nblist.cnt += 1
                    self.nonbonded_nblist.list.append(self.keyupdate[pkey_j])
            self.nonbonded_nblist.index.append(self.nonbonded_nblist.cnt + 1)
        #
        # Add bonds
        for bkey_other,bond_other in new_strcC.bonds.iteritems():
            pkey1 = self.keyupdate[bond_other.pkey1]
            pkey2 = self.keyupdate[bond_other.pkey2]
            bond_i = Bond(pkey1,pkey2)
            self.add_bond(bond_i)
        #
        # Add angles 
        for akey_other,angle_other in new_strcC.angles.iteritems():
            pkey1 = self.keyupdate[angle_other.pkey1]
            pkey2 = self.keyupdate[angle_other.pkey2]
            pkey3 = self.keyupdate[angle_other.pkey3]
            angle_i = Angle(pkey1,pkey2,pkey3)
            self.add_angle(angle_i)
        #
        # Add dihedrals 
        for dkey_other,dihedral_other in new_strcC.dihedrals.iteritems():
            pkey1 = self.keyupdate[dihedral_other.pkey1]
            pkey2 = self.keyupdate[dihedral_other.pkey2]
            pkey3 = self.keyupdate[dihedral_other.pkey3]
            pkey4 = self.keyupdate[dihedral_other.pkey4]
            dihedral_i = Dihedral(pkey1,pkey2,pkey3,pkey4)
            self.add_dihedral(dihedral_i)
        #
        # Add impropers 
        for ikey_other,improper_other in new_strcC.impropers.iteritems():
            pkey1 = self.keyupdate[dihedral_other.pkey1]
            pkey2 = self.keyupdate[dihedral_other.pkey2]
            pkey3 = self.keyupdate[dihedral_other.pkey3]
            pkey4 = self.keyupdate[dihedral_other.pkey4]
            improper_i = Improper(pkey1,pkey2,pkey3,pkey4)
            self.add_improper(improper_i)
            
        return self
                

    def add_struc_grid(self,other,n_i,p=None,tag="blank",verbose=False,calc_overlap = True ):
        """
        Add structure other to self n times on a grid
        
        Args:
            n_i  (int) number of times to replicate 
            
        """
        #
        # Initialize 
        rank = 0
        size = 0
        pos_cut = 10.0 
        ang_acc = 1000  # number of digets in random angle
        exlat_cnt = 0                   # Number of lattice expantions 
        exlat_indx = 0 # Initial lattice vector index to expand 
        exlat_frac = 0.10 # Fraction to expand the lattice by if max_sys is exceeded 
        #
        # Initilize
        strucC_new = copy.deepcopy(self) # Create copy to reset to if addition of other fails 
        strucC_new.lat = self.lat 
        # Set lattice to be recoppied over the strucC_new during initial add
        other.lat = strucC_new.lat
        # 
        # Create a list of atomic indices for each processor 
        # particle_keys = self.particles.keys()  
        # particle_keys_p  = p.splitListOnProcs(particle_keys)
        #
        # Guess size of volume need for each added structure 
        other.calc_mass()
        strucC_new.calc_volume()
        n_vol = strucC_new.volume/float(n_i)  # Number of other strucC per volume
        l_n = n_vol**(1.0/3.0)
        n_lat = [] #np.zeros(strucC_new.lat.n_dim)
        for d in range(strucC_new.lat.n_dim):
            n_lat.append(int(math.ceil( strucC_new.lat._lengths[d]/l_n)))
        if(  rank == 0 ):
            log_line = " Adding %d into volume %f \n"%(n_i,strucC_new.volume)
            log_line += "  %d along v_1 \n"%(n_lat[0])
            log_line += "  %d along v_2 \n"%(n_lat[1])
            log_line += "  %d along v_3 \n"%(n_lat[2])
            log_line += "  %d lattice positions \n"%(n_lat[0]*n_lat[1]*n_lat[2])
            logger.debug(log_line)
        #
        # Start adding molecules to the system
        add_strucC = True
        while ( add_strucC ):
            add_strucC = True
            struc_add_cnt = 0                    # Number of structures add to the new structure
            for indx_i in range(n_lat[0]):
                if( not add_strucC): break 
                for indx_j in range(n_lat[1]):
                    if( not add_strucC): break 
                    for indx_k in range(n_lat[2]):
                        # Calculate fractional coord
                        f_i = float(indx_i)/float(n_lat[0])
                        f_j = float(indx_j)/float(n_lat[1])
                        f_k = float(indx_k)/float(n_lat[2])
                        frac_o  = np.array([f_i,f_j,f_k])
                        pos_o = strucC_new.lat.fractoreal(frac_o)
                        #
                        if ( rank == 0 ):
                            #
                            #  Get random rotation angles from single processor  
                            #
                            rot_angle_i_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                            rot_angle_j_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                        if( size > 0 ):
                            # Broadcast random rotation angles and translations to all processors 
                            p.barrier() # Barrier for MPI_COMM_WORLD
                            rot_angle_i = p.bcast(rot_angle_i_o)
                            rot_angle_j = p.bcast(rot_angle_j_o)
                            r_random = p.bcast(pos_o)
                            p.barrier() # Barrier for MPI_COMM_WORLD
                        else:
                            rot_angle_i = rot_angle_i_o
                            rot_angle_j = rot_angle_j_o
                            r_random = pos_o
                        #
                        # Get coordinates of randomly rotated and shifted
                        other.calc_center_mass()
                        other.shift_pos(-1.0*other.center_mass)  # Place center of mass at origin
                        other.rotate_xy(rot_angle_i)
                        other.rotate_xz(rot_angle_j)
                        other.shift_pos(r_random)  # Place center of mass at random position in lattice 
                        #
                        # Calculate inter particle seperation between particles in the structure and structure to be added
                        overlap = 0                        
                        if( strucC_new.n_particles > 0 and calc_overlap ):
                            npos_i = other.positions
                            npos_j = strucC_new.positions
                            npos_ij,nd_ij = strucC_new.lat.delta_npos(npos_i,npos_j)
                            # If the particles are further than the cut off
                            for key_i in other.particles.keys():
                                for key_j in strucC_new.particles.keys():
                                    if( nd_ij[key_i][key_j] <= pos_cut ):
                                        overlap = 1
                        if( size > 0 ):
                            #
                            # Reduce sum the overlap variable from all the processors
                            #   if it is zero everywhere there was no overlap detected 
                            overlap_sum = p.allReduceSum(overlap)
                            p.barrier() # Barrier for MPI_COMM_WORLD
                        else:
                            overlap_sum = overlap
                        #
                        # If no overlap detected add molecule to the system
                        if( overlap_sum == 0 ):
                            # Create copy of structure to add to preserve original state 
                            structoadd= copy.deepcopy(other)
                            # Update mol number
                            for pkey_i, particle_i  in structoadd.particles.iteritems():
                                particle_i.properties["mol"] = struc_add_cnt
                            struc_add_cnt += 1
                            strucC_new += structoadd
                            logger.info( "Molecule %d/%d added "%(struc_add_cnt,n_i))
                        if( struc_add_cnt == n_i ):
                            add_strucC = False

                            # Record replication
                            strucC_new.tag = tag 
                            replication_i = Replication(self.tag,other.tag,strucC_new.tag,"grid",n_i)
                            strucC_new.replications.append(replication_i)        
                            
                            return strucC_new
            #
            # If attempts to place structure into the new structure failed
            exlat_cnt += 1
            #
            # If the new structure has been reset over max_sys times expand the box size by lc_expand
            matrix_i = copy.deepcopy(strucC_new.lat._matrix) # Save state of lattice
            # Reset new structure to original state 
            strucC_new = copy.deepcopy(self) # Create copy to add new structures to 
            strucC_new.lat.set_matrix(matrix_i)  # Set lattice to saved state
            strucC_new.lat.expand_matrix(exlat_frac)
            #
            # Add one row to a single dimension of the lattice grid 
            n_lat[exlat_indx] += 1
            exlat_indx += 1
            if( exlat_indx > 2 ):
                    exlat_indx = 0
            #
            if( verbose ):
                log_line =  "  %d  added after "%(struc_add_cnt)
                log_line += "%d  exlat_indx "%(exlat_indx)
                log_line += "  n_lat 0 %f "%(n_lat[0])
                log_line += "  n_lat 1 %f "%(n_lat[1])
                log_line += "  n_lat 2 %f "%(n_lat[2])
                log_line += "%d  expantions "%(exlat_cnt)
                print log_line
                            

        return None
                
    def add_struc(self,other,n_i,seed,p=None,tag="blank",verbose=False):
        """
        Add structure other to self n times via random placement
        """
        #
        # MPI setup
        #
        if( p == None ):
            rank = 0
            size = 0
        else:
            rank = p.getRank()
            size = p.getCommSize()
        #
        # Initialize random
        
        if( rank == 0 ):            
            logging.info('Starting add_struc %s to %s with seed %d '%(other.tag,self.tag,seed))
            
        random.seed(seed)
        #
        pos_cut = 10.0 # minimum distance between particles of added structure and current structure
        ang_acc = 1000  # number of digets in random angle
        max_sys = 3
        max_mol_place = 10 
        exlat_frac = 0.10 # Fraction to expand the lattice by if max_sys is exceeded 
        #
        # Initilize
        exlat_cnt = 0                   # Number of lattice expantions 
        reset_cnt = 0                   # Number of times the new structure  has been reset
        struc_add_cnt = 0                # Total number of structures added to the new structure
        strucC_new = copy.deepcopy(self) # Create copy to reset to if addition of other fails
        strucC_new.tag = self.tag 
        strucC_new.lat = self.lat 
        # Set lattice to be recoppied over the strucC_new during initial add
        other.lat = self.lat         
        other.calc_mass()
        # 
        # Create a list of atomic indices for each processor 
        # particle_keys = self.particles.keys()  
        # particle_keys_p  = p.splitListOnProcs(particle_keys)
        #
        #
        # Start adding molecules to the system
        #   
        add_strucC = True
        while ( add_strucC ):
            #
            # Initialize 
            #
            add_strucC = True
            poxpass  = False 
            placement_cnt = 0 
            # 
            # For each structure add to 
            # while ( overlap_sum != 0  ):
            while ( not poxpass ):
                placement_cnt += 1
                if( rank == 0 ):
                    logging.debug('-Placement %d %d '%(placement_cnt,rank))

                rot_angle_i_o = 0.0 
                rot_angle_j_o = 0.0 
                r_random_o  = np.zeros(strucC_new.lat.n_dim) 

                if ( rank == 0 ):
                    #
                    #  Get random rotation angles from single processor  
                    #
                    rot_angle_i_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                    rot_angle_j_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                    #
                    #  Get random translation from single processor 
                    r_random_o = strucC_new.lat.random_pos()
                if( size > 0 ):
                    # Broadcast random rotation angles and translations to all processors 
                    p.barrier() # Barrier for MPI_COMM_WORLD
                    rot_angle_i = p.bcast(rot_angle_i_o)
                    rot_angle_j = p.bcast(rot_angle_j_o)
                    r_random = p.bcast(r_random_o)
                    p.barrier() # Barrier for MPI_COMM_WORLD
                else:
                    rot_angle_i = rot_angle_i_o
                    rot_angle_j = rot_angle_j_o
                    r_random = r_random_o

                if( rank == 0 ):
                    logging.debug('Angle i %f angle j %f on %d '%(rot_angle_i,rot_angle_j,rank))
                #
                # Get coordinates of randomly rotated and shifted
                #
                other.calc_center_mass()
                other.shift_pos(-1.0*other.center_mass)  # Place center of mass at origin
                other.rotate_xy(rot_angle_i)
                other.rotate_xz(rot_angle_j)
                other.shift_pos(r_random)  # Place center of mass at random position in lattice 
                #
                # Calculate inter particle seperation between particles in the structure and structure to be added
                #
                poxpass = True
                if( strucC_new.n_particles > 0 ):
                    npos_i = strucC_new.positions
                    npos_j = other.positions
                    poxpass = strucC_new.lat.proximitycheck(npos_i,npos_j,pos_cut,p=p)
                    
                # If no overlap detected add molecule to the system
                if( poxpass ):
                    if ( rank == 0 ):
                        # logging.debug('No overlap found adding structure %d'%(struc_add_cnt))
                        print 'No overlap found adding structure %s'%(struc_add_cnt)
                    # Create copy of structure to add to preserve original state 
                    structoadd= copy.deepcopy(other)
                    # Update mol number
                    for pkey_i, particle_i  in structoadd.particles.iteritems():
                        particle_i.properties["mol"] = struc_add_cnt
                    struc_add_cnt += 1
                    strucC_new += structoadd
                if( placement_cnt >= max_mol_place ):

                    # If attempts to place molecule into the system exceed max set by max_mol_place
                    #   reset system and star over 
                    if(  rank == 0  ):
                        logging.debug('Max placments %d exceeded resetting to original system '%(max_mol_place))
                    reset_cnt += 1                 # Number of times the new structure  has been reset
                    struc_add_cnt = 0                    # Number of structures add to the new structure 
                    placement_cnt = 0
                    matrix_i = copy.deepcopy(strucC_new.lat._matrix) # Save state of lattice
                    # Reset new structure to original state 
                    strucC_new = copy.deepcopy(self) # Create copy to add new structures to 
                    strucC_new.lat.set_matrix(matrix_i)  # Set lattice to saved state 

                        
                if( reset_cnt >= max_sys  ):
                    if(  rank == 0  ):
                        logging.debug('Max resets %d exceeded expanding box size by %f '%(max_sys,exlat_frac))

                    exlat_cnt += 1
                    # If the new structure has been reset over max_sys times expand the box size by lc_expand
                    matrix_i = copy.deepcopy(strucC_new.lat._matrix) # Save state of lattice
                    # Reset new structure to original state 
                    strucC_new = copy.deepcopy(self) # Create copy to add new structures to 
                    strucC_new.lat.set_matrix(matrix_i)  # Set lattice to saved state
                    strucC_new.lat.expand_matrix(exlat_frac)

                    reset_cnt = 0
                    struc_add_cnt = 0                    # Number of structures add to the new structure 
                    placement_cnt = 0



            if( verbose and   rank == 0  ):

                    log_line =  "  %d  added after "%(struc_add_cnt)
                    log_line += "%d  expantions "%(exlat_cnt)
                    log_line += "%d resets "%(reset_cnt)
                    log_line += "%d placements "%(placement_cnt)
                    logging.info(log_line)


            if( struc_add_cnt ==  n_i  ):
                if(  rank == 0  ):
                    logging.info(" All structures added with %d particles "%(strucC_new.n_particles))
                # If all the molecule have been added exit while loop and print system 
                add_strucC = False


        # Record replication
        strucC_new.tag = tag
        replication_i = Replication(self.tag,other.tag,strucC_new.tag,"random",n_i)
        strucC_new.replications.append(replication_i)        

        if( rank == 0 ):
            logging.info("Finished %s "%(datetime.now()))
           
        return strucC_new



    def calc_elcnt(self,key_i,nb_list):
        '''
        Calculate the number of each element type in the neighbors of particle key_i

        Return:
            el_cnt (list) count of each element at it's atomic number index

            example: el_cnt[6] is the number of carbon neighbors 
        '''

        el_cnt = np.zeros(periodictable.n_elements, dtype =int )
        for key_j in nb_list.getnbs(key_i):
            el_j = int( self.particles[key_j].properties["number"] )
            if( el_j >= 0 ):
                el_cnt[el_j] += 1
                
        return el_cnt

    def change_mass(self,symbol_o,mass_o):
        '''
        Change the mass of the particles with a certain symbol to a certain number
        
        Args:
            symbol_o (str) symbol of particle to be changed 
            mass_o (float) new mass of particle with symbol symbol_o
        
        Reutrns:
            null
            
        '''

        for pkey_i, particle_i  in self.particles.iteritems():
            if( particle_i.properties["symbol"] == symbol_o ):
                particle_i.properties["mass"] = mass_o

    def guess_oplsa(self):
        """
        Guess OPLS-aa atom types based on coordination 
        """
        for pkey_i,particle_i in self.particles.iteritems():
            nb_cnt_i = self.bonded_nblist.calc_nnab(pkey_i)
            el_cnt_i = self.calc_elcnt(pkey_i,self.bonded_nblist)
            #
            # label carbons 
            #
            if particle_i.properties["number"] == 6 :
                if int(nb_cnt_i) == 4 :
                    particle_i.properties["fftype"] = 'CT' # Alkane
                if  int(nb_cnt_i) == 3 :
                    particle_i.properties["fftype"] = 'CA'  # Conjugated 
                if int(nb_cnt_i) == 2 :
                    particle_i.properties["fftype"] = 'C:'   # Allene


                if int(nb_cnt_i) == 1 :
                    particle_i.properties["fftype"] = '' # Aromatic C
                    error_line =  " WARNING!!! carbon index ",pkey_i," bonded to single atom "
                    sys.exit(error_line)
            #
            # label oxygens
            #
            if( particle_i.properties["number"] == 8 ):
                if int(nb_cnt_i) == 1 :
                    particle_i.properties["fftype"] = 'O' # double bonded
                if int(nb_cnt_i) == 2 :
                    particle_i.properties["fftype"] = 'OS' # ether

            #
            # label nitrogens 
            #
            if particle_i.properties["number"] == 7 :
                if int(nb_cnt_i) == 3 :      # amide
                    particle_i.properties["fftype"] = 'N' 

            #
            # label sulfurs
            #
            if( particle_i.properties["number"] == 16 ):
                if int(nb_cnt_i) == 2 :
                    particle_i.properties["fftype"] = 'S'   #  Thioether RSR (UA)


        #
        # label hydrogens
        #
        for pkey_i,particle_i in self.particles.iteritems():
            nb_cnt_i = self.bonded_nblist.calc_nnab(pkey_i)
            el_cnt_i = self.calc_elcnt(pkey_i,self.bonded_nblist)
            if( particle_i.properties["number"] == 1 ):
                if ( nb_cnt_i > 1  ):
                    sys.exit(' over coordinated H')
                if ( nb_cnt_i < 1  ):
                    error_line = ' unbonded H %d '%(pkey_i)
                    sys.exit(error_line)
                key_j = self.bonded_nblist.getnbs(pkey_i)[0]
                particle_j = self.particles[key_j]
                el_cnt_j = self.calc_elcnt(key_j,self.bonded_nblist)

                if ( particle_j.properties["fftype"]== 'CA' ):
                    particle_i.properties["fftype"] = 'HA' #
                if ( particle_j.properties["fftype"]== 'CT' ):
                    particle_i.properties["fftype"] = 'HC' #
        return 




                

    def findbond_key(self,pid_i,pid_j):
        '''
        Find key for bond with particles pid_i and pid_j
        '''
        
        for b_indx,bond_i in self.bonds.iteritems():
            print bond_i.pkey1, bond_i.pkey2,bond_i.border
            
            if( pid_i == bond_i.pkey1 and pid_j == bond_i.pkey2 ):
                return b_indx
            if( pid_i == bond_i.pkey2 and pid_j == bond_i.pkey1 ):
                return b_indx
            
        return None 

            
    def find_pairs(self,list_i,list_j,mol_inter=False,mol_intra=False):
        '''
        Find pairs based on criteria
        '''
        N_i = len(list_i)
        N_j = len(list_j)

        if( N_i == 0 or N_j == 0 ):
            logger.warning(" Empty list passed to structure.find_pairs ")
            return 
        
        # probabilityperpair = 1.0     # Probability per pair i-j 
        logger.info("Finding %d x %d  pairs  "%(N_i,N_j))
        
        pairvalue_ij =  np.zeros((N_i,N_j), dtype=np.float64)   # value assigned to each pair 
        # 
        for indx_i in range(N_i):
            pid_i = list_i[indx_i]
            for indx_j in range(N_j):
                pid_j = list_j[indx_j]
                if( pid_i != pid_j ):
                    pairvalue_ij[indx_i][indx_j] = 1.0
                    if( mol_inter and self.particles[pid_i].properties["mol"] == self.particles[pid_j].properties["mol"] ):
                        pairvalue_ij[indx_i][indx_j] = 0.0
                    elif( mol_intra and self.particles[pid_i].properties["mol"] != self.particles[pid_j].properties["mol"] ):
                        pairvalue_ij[indx_i][indx_j] = 0.0
                    logger.debug(" keyi %d keyj %d has probility value of %f "%(pid_i,pid_j,pairvalue_ij[indx_i][indx_j]))
                        
        return pairvalue_ij


    def get_subdihedrals(self,list_k,list_i,list_j,list_l):
        """
        Find sets of dihedrals in system

        k-i-j-l

        Arguments
            list_k (list) of atom indexes in of the first bonded atom in the dihedral 
            list_i (list) of atom indexes in of the second bonded atom in the dihedral 
            list_j (list) of atom indexes in of the third bonded atom in the dihedral 
            list_l (list) of atom indexes in of the fourth bonded atom in the dihedral
        Return
            sub_dihedrals (list) of four aotms in each dihedral 
        """
        #
        # Find  atom groups k-i-j-l
        #
        subdihedrals = {}
        #
        # Find atom indices  of group i and j
        #
        #for p_k, ptcl_k  in self.ptclC(list_k):
        for pkey_k  in list_k:
            for pkey_i in self.bonded_nblist.getnbs(pkey_k):
                if( pkey_i in list_i ):
                    for pkey_j in self.bonded_nblist.getnbs(pkey_i):
                        if( pkey_j in list_j ):
                            for pkey_l in self.bonded_nblist.getnbs(pkey_j):
                                if( pkey_l in list_l ):
                                    # Check to make sure not in list already
                                    add_dih = True
                                    for dkey_i,dih_i in subdihedrals.iteritems():
                                        a_k = dih_i.pkey1
                                        a_i = dih_i.pkey2
                                        a_j = dih_i.pkey3
                                        a_l = dih_i.pkey4
                                        if( pkey_k == a_k and pkey_i == a_i and pkey_j == a_j and pkey_l == a_l ):
                                            add_dih = False 
                                        if( pkey_k == a_l and pkey_i == a_j and pkey_j == a_i and pkey_l == a_k ):
                                            add_dih = False 
                                    if( add_dih ):
                                        dkey_i =  len(subdihedrals)
                                        dih_i = Dihedral(pkey_k,pkey_i,pkey_j,pkey_l)
                                        subdihedrals[dkey_i] = dih_i
        #
        # Return dictionary to allow multiple instances
        #
        return subdihedrals

# Position manipulations
    def get_pos(self,list_i=[]):

        if( len(list_i) == 0 ):
            list_i = self.particles.keys()
            
        npos_i = []
        for pkey_i in list_i:
            npos_i.append(self.positions[pkey_i])
        return npos_i


# Particle manipulations

    def propcompile_particles(self):
        '''
        Compile all the properties of each particle into a single dictionary
        '''
        if( self.n_particles > 0 ):
            # Create dict of properties 
            self.prop_particles = {}
            # Use first particle's keys 
            part_i = self.particles[0]
            columns = part_i.properties.keys()
            for prop_key in columns:
                self.prop_particles[prop_key] = []

            for pkey_i, particle_i  in self.particles.iteritems():
                for prop_i in columns:
                    self.prop_particles[prop_i].append(particle_i.properties[prop_i])
                    
        return 

    def write_particles(self,particle_file,keys=[]):
        '''
        Write out bonds in csv file 
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.particles.keys()
        
        if( len(keys) > 0 ):
            fout = open(particle_file,'wb')
            writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            header = ['key']
            part_i = self.particles[0]
            for prop_key in part_i.properties.keys():
                header.append(prop_key)
            writer.writerow(header)
            for key_i in keys:
                bond_i = self.bonds[key_i]
                row_i = [key_i,bond_i.pkey1,bond_i.pkey2]
                for prop_key,prop_val in bond_i.properties.iteritems():
                    row_i.append(prop_val)                    
                writer.writerow(row_i)
            fout.close()
        else:
            logger.warning("Empty bond dictionary passed to write_bonds no file will be writen")
        return
        
# Bond manipulation 
    def bonded_bonds(self):
        """
        Generate bonds from bonded neighbor list 
        """
        self.bonds = dict()                                   # Creates empty dict struc
        self.n_bonds = 0    
        for pkey_i in self.particles.keys():
            for pkey_j in self.bonded_nblist.getnbs(pkey_i):
                if( pkey_i < pkey_j ):
                    bond_i = Bond( pkey_i, pkey_j )            
                    self.add_bond(bond_i)


    def calc_bond(self,bond_i):
        """
        Calculate the distacne between two vectors 
        """
        if isinstance(bond_i, Bond):
            r_i = self.positions[bond_i.pkey1]
            r_j = self.positions[bond_i.pkey2]
            r_ij,bond_l = self.lat.delta_pos_c(r_i, r_j)
            bond_i.properties['length'] = bond_l
            
            return bond_l
        else:
            print "Attempting to calculate non-Bond type"
            raise TypeError
        
        return None

    def find_bonds(self,list_i=[],list_j=[]):
        '''
        Find bonds according to list of particles
        '''
        if( len(list_i) == 0 ):
            # If list is not specified use all particles
            list_i = self.particles.keys()
        if( len(list_j) == 0 ):
            # If list is not specified use all particles
            list_j = self.particles.keys()
        keys =[]
        for key_i,bond_i in self.bonds.iteritems():
            calc = False
            if( bond_i.pkey1 in list_i and  bond_i.pkey2 in list_j ):
                calc = True
            if( bond_i.pkey2 in list_i and  bond_i.pkey1 in list_j ):
                calc = True
            if( calc ):
                keys.append(key_i)
                
        return keys
                                    
                    
    def calc_bonds(self,keys=[]):
        '''
        Calculate bond lengths all bond according to list of keys
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.bonds.keys()
        
        for key_i in keys:
            bond_i = self.bonds[key_i]
            self.calc_bond(bond_i)

    def write_bonds(self,bonds_file,keys=[]):
        '''
        Write out bonds in csv file 
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.bonds.keys()
        
        if( len(keys) > 0 ):
            fout = open(bonds_file,'wb')
            writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            header = ['dkey','pkey_i','pkey_j']
            bond_i = self.bonds[0]
            for prop_key in bond_i.properties.keys():
                header.append(prop_key)
            writer.writerow(header)
            for key_i in keys:
                bond_i = self.bonds[key_i]
                row_i = [key_i,bond_i.pkey1,bond_i.pkey2]
                for prop_key,prop_val in bond_i.properties.iteritems():
                    row_i.append(prop_val)                    
                writer.writerow(row_i)
            fout.close()
        else:
            logger.warning("Empty bond dictionary passed to write_bonds no file will be writen")
        return
              
# Angle manipulation
    def bonded_angles(self):
        """
        Generate angles from bonded neighbor list 
        """
        self.angles = dict()
        self.n_angles = 0    
        for pkey_i in self.particles.keys():
            nb_cnt_i = self.bonded_nblist.calc_nnab(pkey_i)
            if( nb_cnt_i >= 2 ):            
                for pkey_j in self.bonded_nblist.getnbs(pkey_i):
                    #if( pkey_i < pkey_j ):
                    for pkey_k in self.bonded_nblist.getnbs(pkey_i):
                        if ( pkey_j < pkey_k ):
                            a_i = Angle( pkey_k,pkey_i, pkey_j )            
                            self.add_angle(a_i)
           
    def calc_angle(self,angle_i):
        """
        Calculate the cosine theta between two vectors 
        """
        if isinstance(angle_i, Angle):
            r_k = self.positions[angle_i.pkey1]
            r_i = self.positions[angle_i.pkey2]
            r_j = self.positions[angle_i.pkey3]
            r_ki = self.lat.norm_delta_pos_c(r_k, r_i)
            r_ij = self.lat.norm_delta_pos_c(r_i, r_j)
            
            r_i_norm = r_ki/np.linalg.norm(r_ki)
            r_j_norm = r_ij/np.linalg.norm(r_ij) 
            cos_kij = np.dot(r_i_norm,r_j_norm)
            
            angle_i.properties['cosine'] = cos_kij
        
            return cos_kij
        else:
            print "Attempting to calculate non-Angle type"
            raise TypeError
        
        return None 

    def find_angles(self,list_k=[],list_i=[],list_j=[]):
        '''
        Find angles according to list of particles
        '''
        if( len(list_k) == 0 ):
            # If list is not specified use all particles
            list_k = self.particles.keys()
        if( len(list_i) == 0 ):
            # If list is not specified use all particles
            list_i = self.particles.keys()
        if( len(list_j) == 0 ):
            # If list is not specified use all particles
            list_j = self.particles.keys()
        keys =[]
        for key_i,angle_i in self.angles.iteritems():
            calc = False
            if( angle_i.pkey1 in list_k and  angle_i.pkey2 in list_i and  angle_i.pkey3 in list_j ):
                calc = True
            if( angle_i.pkey3 in list_k and  angle_i.pkey2 in list_i and  angle_i.pkey1 in list_j ):
                calc = True
            if( calc ):
                keys.append(key_i)
                
        return keys
                                     
    def calc_angles(self,keys=[]):
        '''
        Calculate the cosine theta between two vectors according to list_i and list_j
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.angles.keys()

        for key_i in keys:
            angle_i = self.angles[key_i]
            self.calc_angle(angle_i)
            
    def write_angles(self,angles_file,keys=[]):
        '''
        Write out angles in csv file 
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.angles.keys()
        
        if( len(keys) > 0 ):
        
            fout = open(angles_file,'wb')
            writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            header = ['dkey','pkey_k','pkey_i','pkey_j']
            angle_i = self.angles[0]
            for prop_key in angle_i.properties.keys():
                header.append(prop_key)
            writer.writerow(header)
            for key_i in keys:
                angle_i = self.angles[key_i]
                row_i = [key_i,angle_i.pkey1,angle_i.pkey2,angle_i.pkey3]
                for prop_key,prop_val in angle_i.properties.iteritems():
                    row_i.append(prop_val)                    
                writer.writerow(row_i)
            fout.close()
        else:
            logger.warning("Empty angles dictionary passed to write_angles no file will be writen")

        return


# Dihedral manipulation
                            
    def bonded_dih(self):
        """
        Generate dihedrals from nbonded eighbor list 
        """
        self.dihedrals =dict()                                # Creates empty dict struc
        self.n_dihedrals = 0    
        for pkey_i in self.particles.keys():
            for pkey_j in self.bonded_nblist.getnbs(pkey_i):
                if( pkey_i < pkey_j ):
                    for pkey_k in self.bonded_nblist.getnbs(pkey_i):
                        if ( pkey_k != pkey_j ):
                            for pkey_l in self.bonded_nblist.getnbs(pkey_j):
                                if ( pkey_l != pkey_i and pkey_l != pkey_k ):
                                    d_i = Dihedral( pkey_k, pkey_i, pkey_j, pkey_l )            
                                    self.add_dihedral(d_i)
                                    
    def calc_dihedral(self,dihedral_i):
        '''
        Calculate the dihedral angel 
        '''
        if isinstance(dihedral_i, Dihedral):
            r_k = self.positions[dihedral_i.pkey1]
            r_i = self.positions[dihedral_i.pkey2]
            r_j = self.positions[dihedral_i.pkey3]
            r_l = self.positions[dihedral_i.pkey4]

            r_ki = self.lat.norm_delta_pos_c(r_k, r_i)
            r_ij = self.lat.norm_delta_pos_c(r_i, r_j)
            r_jl = self.lat.norm_delta_pos_c(r_j, r_l)

            v1v2 = np.cross(r_ki,r_ij)
            v2v3 = np.cross(r_ij,r_jl)

            r_i_norm = v1v2/np.linalg.norm(v1v2)
            r_j_norm = v2v3/np.linalg.norm(v2v3) 
            cos_kijl = np.dot(r_i_norm,r_j_norm)
            
            dihedral_i.properties['cosine'] = cos_kijl

            return cos_kijl

            cos_ang = np.arccos(cos_kijl )
            ang_deg = np.rad2deg( cos_ang )

            logger.debug(" Cosine %f angle %f (deg) "%(cos_kijl,ang_deg))

        else:
            print "Attempting to calculate non-dihedral type"
            raise TypeError
        
        return None

    def find_dihedrals(self,list_k=[],list_i=[],list_j=[],list_l=[]):
        '''
        Find dihedrals according to list of particles
        '''
        if( len(list_k) == 0 ):
            # If list is not specified use all particles
            list_k = self.particles.keys()
        if( len(list_i) == 0 ):
            # If list is not specified use all particles
            list_i = self.particles.keys()
        if( len(list_j) == 0 ):
            # If list is not specified use all particles
            list_j = self.particles.keys()
        if( len(list_l) == 0 ):
            # If list is not specified use all particles
            list_l = self.particles.keys()
        keys =[]
        
        for key_i,dihedral_i in self.dihedrals.iteritems():
            calc = False
            if( dihedral_i.pkey1 in list_k and  dihedral_i.pkey2 in list_i and  dihedral_i.pkey3 in list_j and  dihedral_i.pkey4 in list_l ):
                calc = True
            if( dihedral_i.pkey4 in list_k and  dihedral_i.pkey3 in list_i and  dihedral_i.pkey2 in list_j and  dihedral_i.pkey1 in list_l ):
                calc = True
            if( calc ):
                keys.append(key_i)
                
                
        return keys
                            
    def calc_dihedrals(self,keys=[]):
        '''
        Calculate the cosine theta of a dihedral angle according to keys
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.dihedrals.keys()

        for key_i in keys:
            dihedral_i = self.dihedrals[key_i]
            cos_kijl = self.calc_dihedral(dihedral_i)
        
    def write_dihedrals(self,dih_file,keys=[]):
        '''
        Calculate the dihedral angels from angle_list
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.dihedrals.keys()
        
        if( len(keys) > 0 ):
            fout = open(dih_file,'wb')
            writer = csv.writer(fout,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            header = ['dkey','pkey_k','pkey_i','pkey_j','pkey_l']
            dih_i = self.dihedrals[0]
            for prop_key in dih_i.properties.keys():
                header.append(prop_key)
            writer.writerow(header)
            for key_i in keys:
                dih_i = self.dihedrals[key_i]
                row_i = [key_i,dih_i.pkey1,dih_i.pkey2,dih_i.pkey3,dih_i.pkey4]
                for prop_key,prop_val in dih_i.properties.iteritems():
                    row_i.append(prop_val)                    
                writer.writerow(row_i)
            fout.close()
        else:
            logger.warning("Empty dihedrals dictionary passed to write_dihedrals no file will be writen")

        return
          
# Pair manipulations 
    def distbin_pairs(self,list_i,list_j,pairvalue_ij,bin_size,r_cut,rank):
        '''
        Bin distances between particle pairs
        '''
        # Add size to cutoff
        # Assumes cut off is evenly divisable by bin_size 
        #         |                      |
        #    -bin_size/2.0   r_cut   +bin_size/2.0 
        r_cut += bin_size/2.0 

        N_i = len(list_i)
        N_j = len(list_j)


        probabilityperpair = 1
        # Calculate rdf relate values
        n_bins = int(r_cut/bin_size) + 1 
        bin_r = np.zeros(n_bins, dtype=np.float64)
        bin_r_nn = np.zeros(n_bins, dtype=np.float64)    # Nearest neighbor count 
        volumes  = []
        rdf_frames = 0
        coor_i = self.get_pos(list_i)
        coor_j = self.get_pos(list_j)
        
        npos_ij,dist = self.lat.delta_npos(coor_i,coor_j)
        
        self.calc_volume
        volumes.append(self.volume)
        
        for ref_i in range(N_i):
            a_i_hasnieghbor = False
            r_ij_nn = r_cut   # Nearest Neighbor distance  
            for ref_j in range(N_j):
                logger.debug(" Checking pair %d - %d dr %f  pp %f "%(ref_i,ref_j,dist[ref_i,ref_j],pairvalue_ij[ref_i][ref_j] ))
                if(  pairvalue_ij[ref_i][ref_j] > 0.0 ):
                    if(  dist[ref_i,ref_j] <= r_cut ):
                        # bin distance =
                        bin_index = int( round( dist[ref_i,ref_j] / bin_size) )
                        # print " dist / bin / bin_sit", dist[ref_i,ref_j],bin_index,bin_size*float(bin_index)
                            
                        bin_r[bin_index] += pairvalue_ij[ref_i][ref_j]

                        logger.debug(" %f %f "%(bin_r[bin_index] , pairvalue_ij[ref_i][ref_j]))
                        # Find nearest neighbor distance 
                        a_i_hasnieghbor = True
                        if( dist[ref_i,ref_j] < r_ij_nn ):
                            r_ij_nn = dist[ref_i,ref_j]
                            p_ij_nn = pairvalue_ij[ref_i][ref_j]
                            
            # Record nearest neighbor distance 
            if( a_i_hasnieghbor ):
                bin_nn_index = int( round( r_ij_nn /bin_size) )
                bin_r_nn[bin_nn_index] += p_ij_nn

        return bin_r,bin_r_nn,volumes

    

'''
Extra functions  
'''
                      
def prop_list(propkey,keys,dic):
    '''
    Calculate bond lengths all bond according to list of keys
    '''
    if( len(keys) == 0 ):
        # If list is not specified use all bonds
        keys = dic.keys()
    prop_list = []
    for key_i in keys:
        dic_entry = dic[key_i]
        prop_list.append(dic_entry.properties[propkey])

    return prop_list


def hterm_Csp3(hb_length,r_i,r_ij_array):
    """
    Hydrogen terminate segment 

     (j)     (l)
        \   /
         (i)   
        /   \
     (k)     (m)



    """

    add_jk = np.zeros(3)
    for r_ij in r_ij_array:
        add_jk +=  r_ij
    add_jk  = -1.0*add_jk
    add_scale = hb_length*add_jk/np.linalg.norm(add_jk)

    r_l = r_i + add_scale

    return r_l 

def hterm_Csp2(hb_length,r_i,r_ij_array):
    """
    Hydrogen terminate conjugated atom

    (j)  
        \  
         (i)  - (l)
        /  
    (k)
    Hydrogens will be added to sp2 carbons in a plainer configuration

    """
    debug = False

    add_jk = -1.0*( r_ij_array[0] + r_ij_array[1] )
    add_scale = hb_length*add_jk/np.linalg.norm(add_jk)
    r_l = r_i + add_scale

    return r_l


