#! /usr/bin/env python
"""
This module defines the classes relating to a periodic lattice
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"


class Lattice(object):
    
    def __init__(self):
        """        
        Create a lattice 
        """
        self.n_dim = 3 # Number of spatial dimensions is set to 3
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
        #    
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
        #
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

        random.seed(seed) needs to be initialized
        
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
        # 
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
                  