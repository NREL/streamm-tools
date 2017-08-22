# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

"""
This module defines the classes relating to containters of particles 
"""
try:
    from pymatgen.core.lattice import Lattice
except:
    logger.warning("pymatgen import error for Lattice object")
    exit()
    
class Container(object):
    """
    Data structure for describing a collection of Particles that have associated
    positions within a Lattice, and consistent set keys corresponding to
    Bond, Angle, Dihedral and Improper descriptions 
    """

    def __init__(self,tag="blank",matrix=[100.0,0.0,0.0,100.0,0.0,0.0,100.0,0.0,0.0],verbose=False):
        """
        Constructor for a composite structure. 
        """
        if isinstance(tag, str):
            self.tag = tag
        else:
            raise TypeError("1st arg (tag) in %s Container initialization should be string"%(__name__))
        
        self.lat = Lattice(matrix)                                  # Creates lattice object for structure
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

    def write_xyz_prop(self,prop_key,xyz_file=''):
        '''
        Write a list of certain particles of the structure  to an xyz file
        taggd by a property type 

        Args:
            xyz_file    (str) xyz file tag
            
        Reutrns:
            null
            
        '''
        if( len(xyz_file) == 0 ):
            xyz_file = "%s.xyz"%(prop_key)
        # Loop over structures
        F = open(xyz_file,"w")
        F.write(" %d \n" % self.n_particles )
        F.write(" %s \n"%" structure.Container  %s "%(prop_key))
        for pkey_i,particle_i  in self.particles.iteritems():
            pos_i = self.positions[pkey_i]
            tag_i = particle_i.properties[prop_key]
            F.write(" %5s %16.8f %16.8f %16.8f \n"%(tag_i,pos_i[0],pos_i[1],pos_i[2] ))
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
            bond_i.properties = bond_other.properties
            self.add_bond(bond_i)
        #
        # Add angles 
        for akey_other,angle_other in new_strcC.angles.iteritems():
            pkey1 = self.keyupdate[angle_other.pkey1]
            pkey2 = self.keyupdate[angle_other.pkey2]
            pkey3 = self.keyupdate[angle_other.pkey3]
            angle_i = Angle(pkey1,pkey2,pkey3)
            angle_i.properties = angle_i.properties
            self.add_angle(angle_i)
        #
        # Add dihedrals 
        for dkey_other,dihedral_other in new_strcC.dihedrals.iteritems():
            pkey1 = self.keyupdate[dihedral_other.pkey1]
            pkey2 = self.keyupdate[dihedral_other.pkey2]
            pkey3 = self.keyupdate[dihedral_other.pkey3]
            pkey4 = self.keyupdate[dihedral_other.pkey4]
            dihedral_i = Dihedral(pkey1,pkey2,pkey3,pkey4)
            dihedral_i.properties = dihedral_i.properties
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
                logger.debug(log_line)
                            

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

        Example: el_cnt[6] is the number of carbon neighbors
        
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
        
        Args:
            bond_i (Bond) Bond object
        """
        if isinstance(bond_i, Bond):
            r_i = self.positions[bond_i.pkey1]
            r_j = self.positions[bond_i.pkey2]
            r_ij,bond_l = self.lat.delta_pos_c(r_i, r_j)
            bond_i.properties['length'] = bond_l
            
            return bond_l
        else:
            raise TypeError("Attempting to calculate non-Bond type")
        
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
        
    def proc_dihedrals(self,dih_file,keys=[]):
        '''
        Calculate the dihedral angels from angle_list
        '''
        self.dih_dic = {}
        self.dih_dic['cosine'] = []
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.dihedrals.keys()
        if( len(keys) > 0 ):
            for key_i in keys:
                dih_i = self.dihedrals[key_i]
                self.dih_dic['cosine'] .append(dih_i.properties['cosine'] )
                
                
                
                
                
    def write_dihedrals(self,dih_file,keys=[]):
        '''
        Calculate the dihedral angels from angle_list
        '''
        if( len(keys) == 0 ):
            # If list is not specified use all bonds
            keys = self.dihedrals.keys()
        
        if( len(keys) > 0 ):
            fout = open(dih_file,'wb')
            writer = csv.writer(fout,delimiter=',')  #QUOTE_ALL)
            header = ['dkey','pkey_k','pkey_i','pkey_j','pkey_l','mol_i','g_i','g_j','cosine']
            dih_i = self.dihedrals[0]
            for prop_key in dih_i.properties.keys():
                header.append(prop_key)
            writer.writerow(header)
            for key_i in keys:
                dih_i = self.dihedrals[key_i]
                row_i = [key_i,dih_i.pkey1,dih_i.pkey2,dih_i.pkey3,dih_i.pkey4]
                row_i.append(self.particles[dih_i.pkey2].properties['mol'])                    
                row_i.append(self.particles[dih_i.pkey2].properties['group'])                    
                row_i.append(self.particles[dih_i.pkey3].properties['group'])                    
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
        
        Add size to cutoff
        Assumes cut off is evenly divisible by bin_size 
                 |                      |
            -bin_size/2.0   r_cut   +bin_size/2.0
        
        '''
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

    