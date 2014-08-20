"""
Class data structures for atomic data
"""

from particles import Particle, ParticleContainer
from bonds     import Bond,     BondContainer
from angles    import Angle,    AngleContainer
from dihedrals import Dihedral, DihedralContainer

import pbcs
import copy
import numpy as np 
import json
import sys
import os
import groups 


class StructureContainer:
    """
    Data structure for describing a collection of Particles that have associated
    Bond, Angle, Dihedral descriptions

    The contained ParticleContainer object has a special role... this object is a
    dictionary that tracks the globalIDs for all particles in this Structure.

    GlobalIDs for a structure should be unique. All gid-s referred to by the
    BondContainer, AngleContainer etc should be consistent.

    This class will also contain data related to global parameters for a simulation
    (eg box length, pressure, temperature.....)
    """

    def __init__(self, ptclC=ParticleContainer(), bondC=BondContainer(), angleC=AngleContainer(), dihC=DihedralContainer(), verbose=True):
        """
        Constructor for a composite structure. Deepcopy of containers is used
        so this is the main container that manages memory for all sim objects

        Args:
            ptclC  (ParticleContainer)  
            bondC  (BondContainer)  
            angleC (AngleContainer)
            dihC   (DihedralContainer)
        """
        
        if isinstance(ptclC, ParticleContainer):
            self.ptclC = copy.deepcopy(ptclC)
        else:
            print "1st arg should be a ParticleContainer object"
            raise TypeError

        if isinstance(bondC, BondContainer):
            self.bondC = copy.deepcopy(bondC)
        else:
            print "2nd arg should be a BondContainer object"
            raise TypeError

        if isinstance(angleC, AngleContainer):
            self.angleC = copy.deepcopy(angleC)
        else:
            print "3rd arg should be an AngleContainer object"
            raise TypeError

        if isinstance(dihC, DihedralContainer):
            self.dihC = copy.deepcopy(dihC)
        else:
            print "4rd arg should be an DihedralContainer object"
            raise TypeError


        # Debug print flag
        self.verbose = verbose

        # Length of cartesian box size [ [xmin, xmax], [ymin, ymax]...]
        self.boxLengths = [ [0.0, 1.0], [0.0, 1.0], [0.0, 1.0] ]

        # Lattice vectors 
        self.latvec = [  np.array([100.0,0.0,0.0]) ,  np.array( [0.0,100.0,0.0]),  np.array( [0.0,0.0,100.0]) ]


    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        if self.verbose:
            print "Cleaning structureContainer"
            
        del self.ptclC
        del self.bondC
        del self.angleC
        del self.dihC


    def __len__(self):
        """
        'Magic' method for returning size of container. This is defined as the
        size of the particle container
        """
        return len(self.ptclC)


    def dump(self, filePrefix):
        """
        Dump a pickled version of this object

        Args:
            filePrefix (str): name of pickle file. will dump filePrefix.pkl
        """

        if self.verbose:
            print "Dumping structure container to pickle file"

        import pickle 
        fileObj = open(filePrefix + '.pkl', 'w')
        pickle.dump(self, fileObj)
        fileObj.close()


    def restore(self, pickleName):
        """
        Restore state of object from a pickled file. Sets information
        in an instance of this class

        Args:
            pickleName (str): Complete name of pickle file.
        """
        
        import pickle 
        fileObj = open(pickleName, 'r')
        struc = pickle.load(fileObj)
        fileObj.close()

        # Custom restore for each data member
        self.ptclC  = copy.deepcopy(struc.ptclC)
        self.bondC  = copy.deepcopy(struc.bondC)
        self.angleC = copy.deepcopy(struc.angleC)
        self.dihC   = copy.deepcopy(struc.dihC)

        self.boxLengths = copy.deepcopy(struc.boxLengths)
        self.latvec     = copy.deepcopy(struc.latvec)


    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """

        strucStr =  "\n"
        strucStr += "************************************* \n"
        strucStr += " Global parameters: \n"
        strucStr += "************************************* \n"
        strucStr += "  Box lengths: \n"
        strucStr += "      Lx = " + str(self.boxLengths[0]) + "\n"
        strucStr += "      Ly = " + str(self.boxLengths[1]) + "\n"
        strucStr += "      Lz = " + str(self.boxLengths[2]) + "\n"
        strucStr += "\n"                                        
        strucStr += "************************************* \n"
        strucStr += " Structure contains: \n"
        strucStr += "************************************* \n"
        strucStr += str(self.ptclC)
        strucStr += str(self.bondC)
        strucStr += str(self.angleC)
        strucStr += str(self.dihC)
        return strucStr


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator (eg struc1 += struc2)
        
        Compare global IDs of particles and reassign globalIDs for particle
        container using the max ID between the two lists. Tracks these changes
        for all other (bond, angle, dihedral) containers that reference particleIDs
        """

        # Empty container checks
        if len(other) == 0:             # If struc2 is empty (has no particles)
            return self                 #   simply return unchanged current container
        if len(self) == 0:              # If struc1 (this struc) is empty (has no particles)
            return copy.deepcopy(other) #    full copy other and return

        idFromToDict = dict()  # Need to keep track of all ptcl ID changes at once
                               # {fromID1:toID1, fromID2:toID2...}
                               # eg {1:3, 3:5, 2:20...}
        
        bondC  = BondContainer()              # Local bond container copy so ptclIDs
        angleC = AngleContainer()             # Local angle container copy so ptclIDs
        dihC   = DihedralContainer()          # Local dihedral container copy so ptclIDs
        bondC  = copy.deepcopy(other.bondC)   #  inside can be changed (for adding below)
        angleC = copy.deepcopy(other.angleC)  #  inside can be changed (for adding below)
        dihC   = copy.deepcopy(other.dihC)    #  inside can be changed (for adding below)

        keys1 = self.ptclC.particles.keys()    # global IDs of particles in this object
        keys2 = other.ptclC.particles.keys()   # global IDs in object being added
        self.ptclC.maxgid= max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        for ptclkey2 in other.ptclC.particles:
            self.ptclC.put(other.ptclC.particles[ptclkey2]) # Pushes ptcl to this struc's ptcl container
            fromPtclID = ptclkey2                           # Track IDs from--->to
            toPtclID   = self.ptclC.maxgid                  #  --> toID (to is the maxid of this ptclC)
            idFromToDict[fromPtclID]=toPtclID               # Store ID changes


        bondC.replacePtclIDs(idFromToDict)  # Use tracked list of ID changes
        self.bondC += bondC                 # Now add bondC with 'corrected' IDs

        angleC.replacePtclIDs(idFromToDict) # Now add angleC with 'corrected' IDs
        self.angleC += angleC               # Use tracked list of ID changes

        dihC.replacePtclIDs(idFromToDict)   # Use tracked list of ID changes
        self.dihC += dihC                   # Now add dihC with 'corrected' IDs

        return self



    def compressPtclIDs(self):
        """
        Replace all particle IDs such that if there are N particles in structure
        the particle ID's run from 1...N. Tracks these changes
        for all other (bond, angle, dihedral) containers that reference particleIDs
        """

        idFromToDict = dict()  # Need to keep track of all ptcl ID changes at once
                               # {fromID1:toID1, fromID2:toID2...}
                               # eg {1:3, 3:5, 2:20...}

        localDict = dict() # Local copy for reordering
        
        for toPtclID, ptclTuple in enumerate(self.ptclC): # Enumerate returns (ID, obj) tuple for ptclTuple
            toPtclID +=1                                  # Sets reordering index correctly
            fromPtclID = ptclTuple[0]                     # Picks out ID from ptclTuple
            idFromToDict[fromPtclID]=toPtclID             # Store ID changes
            ptclObj = self.ptclC[fromPtclID]              # Get particle object
            localDict[toPtclID] = ptclObj                 # Set local dictionary with reordered index

        del self.ptclC.particles              # Ensure memory is free
        self.ptclC.particles = localDict      # Assign reordered local copy to ptcl container
        del localDict                         # Clear local memory

        # ----------- Redo bonds -------------
        self.bondC.replacePtclIDs(idFromToDict)         # Use tracked list of ID changes bonds

        localDict = dict()
        for toBondID, bondTuple in enumerate(self.bondC):   # Enumerate returns (ID, obj) tuple for ptclTuple
            toBondID +=1                                    # Sets reordering index correctly
            fromBondID = bondTuple[0]                       # Picks out ID from ptclTuple
            bondObj = self.bondC[fromBondID]                # Get bond object
            localDict[toBondID] = bondObj                   # Set local dict with reordered index

        del self.bondC.bonds                       # Ensure memory is free
        self.bondC.bonds = localDict               # Assign reordered local copy to bond container
        del localDict

        # ------------ Redo angles -------------
        self.angleC.replacePtclIDs(idFromToDict)        # Use tracked list of ID changes angles

        localDict = dict()
        for toAngleID, angleTuple in enumerate(self.angleC): # Enumerate returns (ID, obj) tuple for ptclTuple
            toAngleID +=1                                    # Sets reordering index correctly
            fromAngleID = angleTuple[0]                      # Picks out ID from ptclTuple
            angleObj = self.angleC[fromAngleID]              # Get angle object
            localDict[toAngleID] = angleObj                  # Set local dict with reordered index

        del self.angleC.angles                    # Ensure memory is free
        self.angleC.angles = localDict            # Assign reordered local copy to bond container
        del localDict

        # ------------ Redo dihedrals-----------
        self.dihC.replacePtclIDs(idFromToDict)  # Use tracked list of ID changes angles

        localDict = dict()
        for toDihedralID, dihedralTuple in enumerate(self.dihC): # Enumerate returns (ID, obj) tuple for ptclTuple
            toDihedralID +=1                                     # Sets reordering index correctly
            fromDihedralID = dihedralTuple[0]                    # Picks out ID from ptclTuple
            dihedralObj = self.dihC[fromDihedralID]              # Get dihedral object
            localDict[toDihedralID] = dihedralObj                # Set local dict with reordered index

        del self.dihC.dihedrals                    # Ensure memory is free
        self.dihC.dihedrals = localDict            # Assign reordered local copy to bond container
        del localDict



    def getSubStructure(self, ptclIDList):
        """
        Return a new Structure object with partcleID's in input list
        Preserves IDs of particles (and any bonds, angles, dihedrals...)
        Bonds, angles, dihedrals... etc are included if and ONLY if all the
        particles of which it consists is in the ptclIDList

        Args:
            ptclIDList (list) global particles ID's for which to return structure

        Return:
            New Structure() object. IDs in new object are unique
        """

        if (len(self)==0 and len(ptclIDList)>0):
            print "Error: getSubStructure using non-zero ptcl list on empty container"
            sys.exit(0)

        subAtoms = ParticleContainer(ptclIDList) # Initial ptcl container w/input IDs

        bondIDList = self.bondC.keys()           # Get keys of bond container
        subBonds = BondContainer(bondIDList)     # Intitialize subbond container

        angleIDList = self.angleC.keys()         # Get keys of angle container
        subAngles = AngleContainer(angleIDList)  # Initialize subangle container

        dihedralIDList = self.dihC.keys()                 # Get keys of dihedral container
        subDihedrals = DihedralContainer(dihedralIDList)  # Initialize subdihedral container        

        # Grab particles from IDlist and put into sub-particle container
        for pgid in ptclIDList:
            atom = self.ptclC[pgid]
            subAtoms[pgid] = atom

        # For each bond object in container check that both
        # particles in bond are in ptcl search list
        for gid, bondObj in self.bondC:
            if ( (bondObj.pgid1 in ptclIDList) and \
                 (bondObj.pgid2 in ptclIDList) ):
                subBonds[gid] = bondObj
            else:
                # Need to remove empty key generated above
                del subBonds[gid]

        # For each angle object in container check that both
        # particles in angle are in ptcl search list
        for gid, angleObj in self.angleC:
            if ( (angleObj.pgid1 in ptclIDList) and \
                 (angleObj.pgid2 in ptclIDList) and \
                 (angleObj.pgid3 in ptclIDList) ):
                subAngles[gid] = angleObj
            else:
                # Need to remove empty key generated above
                del subAngles[gid]

        # For each dihedral object in container check that both
        # particles in dihedral are in ptcl search list
        for gid, dObj in self.dihC:
            if ( (dObj.pgid1 in ptclIDList) and \
                 (dObj.pgid2 in ptclIDList) and \
                 (dObj.pgid3 in ptclIDList) and \
                 (dObj.pgid4 in ptclIDList) ):
                subDihedrals[gid] = dObj
            else:
                # Need to remove empty key generated above
                del subDihedrals[gid]

        return StructureContainer(subAtoms, subBonds, subAngles, subDihedrals)



    def setPtclPositions(self, ptclPosList):
        """
        Reset all particle positions from an indexed list

        Args:
            ptclPosList (list) [ [1, 0.5, 0.1, 12.0],
                                 [2, 0.4, 33.3, -0.1] .....]
        """

        for ptclPos in ptclPosList:
            index = ptclPos[0]
            pos = ptclPos[1:]
            ptclObj = self.ptclC.particles[index]
            ptclObj.position = pos


    def getPtclPositions(self):
        """
        Return list of indexed particle positions
        
        Return: (list) [ [1, 0.5, 0.1, 12.0],
                         [2, 0.4, 33.3, -0.1] .....]
        """

        idxPosList = list()
        
        for pid, ptclObj in self.ptclC:
            pos = ptclObj.position       # Get pos directly from ptcl object
            idxPos = copy.deepcopy(pos)  # Local copy for editing
            idxPos.insert(0, pid)        # Prepend index to pos eg  [1, 0.5, 0.1, 12.0],
            idxPosList.append(idxPos)    # Collect list of index positions
            
        return idxPosList


    def setBoxLengths(self, bLs):
        """
        Set length of cartesian box size... [ [xmin, xmax], [ymin, ymax]...]

        Args:
            bLs (list) box lengths (x,y,z) values
        """

        if isinstance(bLs, list):
            self.boxLengths = bLs
        else:
            print "box lengths should be a list"
            raise TypeError


    def getBoxLengths(self):
        """
        Return: list of cartesian box lengths [units ?]
        """
        return self.boxLengths

    def setLatVec(self, latvec_list ):
        """
        Set length of lattice vector 
        """
        self.latvec = latvec_list


    def readOutput(self, fileName):
        """
        This is the 'effective' base class interface for a method
        that reads in an external output file and populates an instance
        of this class object

        This method should be redefined for each kind of file types
        (typically defined by simulation version eg LAMMPS, Gaussian etc)
        The derived classes must implement the following:
        
        def readOutput(self, fileName):
          ...
          ...
        return None

        Args:
            fileName (str) string of filename to input
        """

        print "No StructureContainer:readOutput method defined for pure base class"
        sys.exit(0)


    def writeInput(self, fileName):
        """
        This is the 'effective' base class interface for a method
        that writes an input file based on the internal attributes of an instance
        of the StructureContainer

        This method should be redefined for each kind of file types
        (typically defined by simulation version eg LAMMPS, Gaussian etc)
        The derived classes must implement the following:
        
        def writeOutput(self, fileName):
          ...
          ...
        return None

        Args:
            fileName (str) string of filename to input
        """

        print "No StructureContainer:writeInput method defined for pure base class"
        sys.exit(0)



    #########################################################

    def getpartnumb(self):
        """
        Return number of particles in a structure 
        """
        NP = 0
        for pid, ptclObj in self.ptclC :
            NP += 1

        return NP


    def expandLatVec(self, precent ):
        """
        Expand lattice vector by a certain percent 
        """
        self.latvec[0] = self.latvec[0]*(1.0+precent)
        self.latvec[1] = self.latvec[1]*(1.0+precent)
        self.latvec[2] = self.latvec[2]*(1.0+precent)

    def getLatVec(self):
        """
        Get lattice vector 
        """
        
        return self.latvec

    def getVolume(self):
        """
        Calculate volume

        Method:
            Volume = ( v_i x v_j ) \dot v_k
        """

        br1 = np.cross(self.latvec[0],self.latvec[1])
        vol = np.dot(br1,self.latvec[2])

        return vol

    def getTotMass(self):
        """
        Calculate total mass of system 

        """
        # Sum mass, charges
        total_mass = 0.0
        for pid, ptclObj in self.ptclC :
            total_mass += ptclObj.mass

        return float(total_mass)

    def getDen(self):
        """
        Calculate density of system 

        """
	# const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

	volume_i = self.getVolume()    
	total_mass_i = self.getTotMass()

	density_i = total_mass_i/volume_i #/const_avo*10.0
	
	return density_i

    def getlength(self):
        """
        Calculate length from maximum seperation between particles

        Return
          struc_len (float) 
        
        """
        
        sq_maxdr = -1000000.0 
        for p_i, ptclObj_i in self.ptclC :
            r_i = np.array( ptclObj_i.position )
            for p_j, ptclObj_j in self.ptclC :
                r_j = np.array( ptclObj_j.position )
                r_ij_sq = pbcs.sq_drij_c(r_i,r_j,self.latvec)
                if( r_ij_sq > sq_maxdr):
                    sq_maxdr = r_ij_sq


        struc_len = np.sqrt(sq_maxdr)

        return struc_len

    def vec_shift(self,r_shift):
        """
        Shift structure by vector

        Arguments
          r_shift (numpy vector) to shift all the cordinates by
          
        """
    
        for pid, ptclObj in self.ptclC :
             r_i = np.array( ptclObj.position ) + r_shift
             ptclObj.position = [r_i[0],r_i[1],r_i[2]]

            
    def center_mass(self):
        """
        Find center of mass of a structure

        Return
          r_mass (numpy array) position of the center of mass
          
        """
        import numpy as np

        total_mass_i = self.getTotMass()

        r_mass = np.array( [0.0,0.0,0.0] )
    
        
        for pid, ptclObj in self.ptclC :
             r_mass += ptclObj.mass*np.array( ptclObj.position )

        r_mass = r_mass/total_mass_i

        return r_mass

            
    def shift_center_mass(self,r_shift):
        """
        Translate center of mass of a structure to a location 

        Return
          r_shift (numpy array) position of the center of mass
          
        """
        
        r_mass = self.center_mass()
        r_m_s = r_shift - r_mass
        self.vec_shift(r_m_s)
        
    def rotate(self,rot_angle_i,rot_angle_j):
        """
        Rotate particles in particle container

        Arguments
          rot_angle_i (float)  0 - pi 
          rot_angle_j (float)  0 - pi 

        """
        import numpy as np 
        import math

        # set variables of rotation matrix
        #   for rotation i around y axis 
        cy = math.cos(rot_angle_i)
        sy = math.sin(rot_angle_i)
        #   for rotation j around z axis 
        cz = math.cos(rot_angle_j)
        sz = math.sin(rot_angle_j)
        # loop over each particle 
        for pid, ptclObj in self.ptclC :
            xd = ptclObj.position[0]
            yd = ptclObj.position[1]
            zd = ptclObj.position[2]
            # Apply rotation matrix i and j to get new postion 
            r_x =  cy*cz*xd - sz*cy*yd + sy*zd 
            r_y =  sz*xd    + cz*yd            
            r_z = -sy*cz*xd + sy*sz*yd + cy*zd

            # ptclObj.position = numpy.array( [r_x,r_y,r_z] )
            ptclObj.position = [r_x,r_y,r_z]
        
    def printprop(self):
        """
        Print properties of a structure 
        """
        print "  Particles %d "%(len(self.ptclC))
        print "    Volume %f "%self.getVolume()
        print "    Mass %f "%self.getTotMass()
        print "    Density %f "%self.getDen()
        print "  Lattice vectors "
        print "    v_i ",self.latvec[0]
        print "    v_j ",self.latvec[1]
        print "    v_k ",self.latvec[2]
        print "  Bonds %d "%(len(self.bondC))


    def maxminLatVec(self):
        """
        Set lattice vector based on the max min postion of the particles 
        """
        x_max = -100000.0
        y_max = -100000.0
        z_max = -100000.0
        
        x_min = 100000.0
        y_min = 100000.0
        z_min = 100000.0
        l_max = -100000.0
        
        for pid, ptclObj in self.ptclC :

            r_x = float(ptclObj.position[0])
            r_y = float(ptclObj.position[1])
            r_z = float(ptclObj.position[2])

            if( r_x > x_max): x_max = r_x
            if( r_y > y_max): y_max = r_y
            if( r_z > z_max): z_max = r_z
                
            if( r_x < x_min): x_min = r_x
            if( r_y < y_min): y_min = r_y
            if( r_z < z_min): z_min = r_z


        if( (x_max - x_min) > l_max ): l_max = (x_max - x_min)
        if( (y_max - y_min) > l_max ): l_max = (y_max - y_min)
        if( (z_max - z_min) > l_max ): l_max = (z_max - z_min)
            
        self.latvec[0][0] = l_max
        self.latvec[1][1] = l_max
        self.latvec[2][2] = l_max
        

    def getchainnumb(self):   # Move out of class
        """
        Return number of chains in a structure 
        """
        n_chains = 0
        for pid, ptclObj in self.ptclC :
            if( ptclObj.tagsDict["chain"] > n_chains): n_chains = ptclObj.tagsDict["chain"]

        return n_chains

    def putstruc_json(self, json_data ):  # Move out of class
        """
        Write a structure into json file
                
        Args:
            json_data (json) json data structure 
            json_file (srt) name of json file 

        """        
        # Initialize json data 
        #   
        struc_data = {}        # Data for entire structure  
        particle_data = {}     # Data for particles and positions 
        twobody_data = {}   # Data for connections between particles (two body interactions) 
        threebody_data = {}        # Data for angles between particles (three body interactions) 
        fourbody_data = {}     # Data for dihedrals between particles  (four body interactions) 
        
        json_data["structure"] = struc_data
        struc_data["particle"] = particle_data
        struc_data["twobody"] = twobody_data
        struc_data["threebody"] = threebody_data
        struc_data["fourbody"] = fourbody_data
    
	# Structure data
        lv_string = str( "%f %f %f %f %f %f %f %f %f " % ( self.latvec[0][0], self.latvec[0][1], self.latvec[0][2], self.latvec[1][0], self.latvec[1][1], self.latvec[1][2], self.latvec[2][0], self.latvec[2][1], self.latvec[2][2]))

	struc_data["latvector"] = lv_string

	# Particle data 
        particle_data["number_id"] = []
        particle_data["type"] = []
        particle_data["position"] = []
        particle_data["mass"] = []
        particle_data["charge"] = []
        particle_data["chain"] = []
        particle_data["ring"] = []
        particle_data["resname"] = []
        particle_data["residue"] = []
        particle_data["linkid"] = []
        particle_data["fftype"] = []

	# Loop over particles and ad them to json data 
        for  pid, ptclObj in self.ptclC:
            particle_data["number_id"].append(pid )
            particle_data["type"].append( ptclObj.type )
            particle_data["position"].append( ptclObj.position )
            particle_data["mass"].append( ptclObj.mass )
            particle_data["charge"].append( ptclObj.charge )
            # Dictionary items 
            particle_data["chain"].append( ptclObj.tagsDict["chain"] )
            particle_data["ring"].append( ptclObj.tagsDict["ring"] )
            particle_data["resname"].append( ptclObj.tagsDict["resname"] )
            particle_data["residue"].append( ptclObj.tagsDict["residue"] )
            particle_data["linkid"].append( ptclObj.tagsDict["linkid"] )        
            particle_data["fftype"].append( ptclObj.tagsDict["fftype"] )        

        twobody_data["bonds"] = []
        for b_i,bondObj in  self.bondC:
            pt_i = bondObj.pgid1
            pt_j = bondObj.pgid2
            twobody_data["bonds"].append( [pt_i,pt_j])
            
        return json_data


    def write_json(self,dir_id,output_id ):
	"""
	Write structure information into a json file
	"""
	
	json_data = {}
	json_data = self.putstruc_json(json_data)


	json_file = dir_id+"/"+output_id + ".json"
	f = open(json_file, 'w')
        json.dump(json_data,f, indent=2)
        f.close()

    def getsys_json(self, json_file):
        """
        Read in structure information from json file 

        Args:
            json_file (json) file with structure information

        Return
          json_data (json structure data )
          
        """

        f = open(json_file, 'r')
        json_data = json.load(f)
        f.close()

        # Place paticle data in sperate data structure 
        struc_data = json_data["structure"]
        particle_data = json_data["structure"]["particle"]

        # Create structure container for particles
        

        for p_i in range( len( particle_data["number_id"])):
            r_i = particle_data["position"][p_i]
            atomic_symb = str( particle_data["type"][p_i] )
            m_i = float(particle_data["mass"][p_i])
            q_i = float(particle_data["charge"][p_i])
            # Create particle
            pt_i = Particle( r_i,atomic_symb,q_i,m_i )
            # Find needed tags
            chain_i = int( particle_data["chain"][p_i] )
            ring_i = particle_data["ring"][p_i]
            resname_i = particle_data["resname"][p_i]
            residue_i = particle_data["residue"][p_i]
            linkid_i = particle_data["linkid"][p_i]
            fftype_i = particle_data["fftype"][p_i]
            # _i = particle_data[""][p_i]
            # Add particle to structure 
            tagsD = {"chain":chain_i,"ring":ring_i,"resname":resname_i,"residue":residue_i,"linkid":linkid_i,"fftype":fftype_i}
            pt_i.setTagsDict(tagsD)
            self.ptclC.put(pt_i)

        # Read in bonds
        twobody_data =  json_data["structure"]["twobody"]
        bondC_i = BondContainer()
        for b_indx in range( len(twobody_data["bonds"] )):
            a_i = int(twobody_data["bonds"][b_indx][0] )
            a_j = int(twobody_data["bonds"][b_indx][1] )
            b_i = Bond( a_i, a_j )            
            self.bondC.put(b_i)

        # Read in lattice vectors
        self.latvec = []
        lv_array = struc_data["latvector"].split()
        self.latvec.append(  np.array( [float(lv_array[0]),float(lv_array[1]),float(lv_array[2])] ) )
        self.latvec.append(  np.array( [float(lv_array[3]),float(lv_array[4]),float(lv_array[5])] ) )
        self.latvec.append(  np.array( [float(lv_array[6]),float(lv_array[7]),float(lv_array[8])] ) )

        return json_data
    
    def create_top(self,ff_charges): # Move out of class (or derived class)
        """
        Find topology information for force-field input files 
        """

        # Version 1 will be dependent on Atomicpy
        import elements , lammps ,gromacs , atom_types, top 

        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        GTYPE = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( [ float( ptclObj.position[0] ),float( ptclObj.position[1] ),float( ptclObj.position[2] )]  )  )
            AMASS.append( float( ptclObj.mass)  )
            CHARGES.append( float( ptclObj.charge ) )
            MOLNUMB.append( int( ptclObj.tagsDict["chain"] ) )
            RESID.append( int( ptclObj.tagsDict["residue"] ) )
            RESN.append(  ptclObj.tagsDict["resname"]  )
            GTYPE.append( ptclObj.tagsDict["gtype"] )
            
            
        # Direct copy of top.print_ff_files
            

        # Read in ff file
        #itp_file = 'ff.itp'
        ##print "   Read in parameters from ",itp_file
        #FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(itp_file)
                 
        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        verbose = True 

        # Find atomic number based on atomic symbol 
        ELN = elements.asymb_eln(ASYMB)
        AMASS = elements.eln_amass(ELN)

        #   Initialize  topology values
        #GTYPE = top.initialize_gtype( ELN )  # need to do this in frag read in 
        CHARN = top.initialize_charn( ELN )        

        NA = len(ELN)
        
        find_bonds = False 
        if( find_bonds ):

            #   Build covalent nieghbor list for bonded information 
            NBLIST, NBINDEX = top.build_covnablist(ELN,R)
            BONDS = top.nblist_bonds(NA,NBLIST, NBINDEX)

            print_bonds = True
            if( print_bonds ):
                for b_indx in range(len(BONDS)):
                    print " bond ",BONDS[b_indx][0]+1,BONDS[b_indx][1]+1

                sys.exit(" print bonds from proximity ")
        else:
            BONDS = []
            for b_i,bondObj in  self.bondC:
                pt_i = bondObj.pgid1
                pt_j = bondObj.pgid2
                BONDS.append( [pt_i-1,pt_j-1])
                print "create_top  bond ",pt_i,pt_j

            NBLIST,NBINDEX = groups.build_nablist_bonds(ELN,BONDS)
            #sys.exit(" passing bonds checking ")
            

        print_bonds = True
        if( print_bonds ):
            for b_indx in range(len(BONDS)):
                print " bond ",BONDS[b_indx][0]+1,BONDS[b_indx][1]+1
            
        
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
        #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)

        #
        # Set charge groups
        #
        CG_SET = []
        one = 1
        for i in range( len(ELN) ):
            CG_SET.append(one)
        #

        if( verbose ):
            print "      Finding atom types  "
            if( ff_charges ):
                print "      Using ff charges "

        d_mass = 0
        d_charge = 0

        find_rings = False 
        NA = len(ELN)
        if( find_rings ):
            RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
        else:
            zero = 0.0

            RINGLIST = []
            RINGINDEX = []
            RING_NUMB = []

            # relabel based on neighbors
            for i in range(NA):
                RINGLIST.append(zero)
                RING_NUMB.append(zero)
                RINGINDEX.append(zero)

            RINGLIST.append(zero)
            RINGINDEX.append(zero)

        # Asign oplsaa atom types
        ATYPE, CHARGES = atom_types.oplsaa( ff_charges,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB)

        #ATYPE , CHARGES = atom_types.biaryl_types( ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
        ATYPE ,RESID, CHARGES = atom_types.set_pmmatypes(ff_charges, ELN, ATYPE,GTYPE,RESID,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB )

        
        ATYPE,RESID,CHARGES,CG_SET,CHARN = atom_types.set_ptmatypes( ff_charges, ELN,ASYMB, ATYPE,GTYPE,RESID,CHARGES,AMASS,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB,CG_SET,CHARN )

        debug = False 
        if(debug):
            # relabel based on neighbors
            NA = len(ELN)
            for i in range(NA):
                print  i+1,ATYPE[i] ,RESID[i], CHARGES[i]
            sys.exit(" atom chagre check 1 ")

        #Refind inter ring types
        ATYPE , CHARGES  = atom_types.interring_types(ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )

        # Update structure information
        pt_cnt = 0
        for pid, ptclObj  in self.ptclC:
             ptclObj.tagsDict["fftype"] = ATYPE[pt_cnt]
             ptclObj.mass = AMASS[pt_cnt]
             ptclObj.charge = CHARGES[pt_cnt]
             ptclObj.tagsDict["ring"] = RING_NUMB[pt_cnt]
             pt_cnt += 1 

        # Add bonds to system
        #for i in range( len(BONDS) ):
        #    #
        #    a_i = BONDS[i][0] 
        #    a_j = BONDS[i][1]
        #    b_i = Bond( a_i+1, a_j+1 )
        #    self.bondC.put(b_i)

            #  Sudo code
            #   Read in parameter file
            #      gromacs itp file
            #      tinker parameter file
            #      lammps parameters in data file
            #   Find bonds
            #      from gaussian output optimization
            #      distance cut-off
            #         system_i = system_i.bonds()
            #   Find Rings
            #   Guess atom types
            #      amber
            #      oplsaa
            #         oligomer = oligomer.guess_oplsaatypes()



    def write_gro(self,dir_id,output_id ): # Move out of class
        """
        Write out gromacs gro file
        """
        # Version 1 will be dependent on Atomicpy
        import gromacs 

        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        
        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        ATYPE = []
        RING_NUMB = []
        GTYPE = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( ptclObj.position)  )
            AMASS.append( float(ptclObj.mass)  )
            CHARGES.append( float(ptclObj.charge)  )
            MOLNUMB.append( int(ptclObj.tagsDict["chain"])  )
            RESID.append( ptclObj.tagsDict["resname"][0:5]  )
            RESN.append( int(ptclObj.tagsDict["residue"])  )
            ATYPE.append( ptclObj.tagsDict["fftype"]  )
            RING_NUMB.append( int(ptclObj.tagsDict["ring"])  )
            GTYPE.append( ptclObj.tagsDict["gtype"]  )


            print " GTYPE ", ptclObj.tagsDict["gtype"] 
            print " RESID ", ptclObj.tagsDict["resname"] 
            print " RESN ", ptclObj.tagsDict["residue"] 
       
        # Set cubic lattice constant to 5 nm arbitrary 
        LV = np.zeros( (3,3) )
            
        LV[0][0] = self.latvec[0][0]
        LV[1][1] = self.latvec[1][1]
        LV[2][2] = self.latvec[2][2]
        
        out_gro = dir_id+"/"+output_id + ".gro"
        gromacs.print_gro(out_gro,GTYPE,RESID,RESN,R,LV)

    def write_top(self,dir_id,output_id,norm_dihparam,itp_file ): # Move out of class
        """
        Write out gromacs gro file
        """
        # Version 1 will be dependent on Atomicpy
        import gromacs , elements, top , lammps, groups , atom_types

        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        
        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        ATYPE = []
        RING_NUMB = []
        GTYPE = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( ptclObj.position)  )
            AMASS.append( float(ptclObj.mass)  )
            CHARGES.append( float(ptclObj.charge)  )
            MOLNUMB.append( int(ptclObj.tagsDict["chain"])  )
            RESID.append( ptclObj.tagsDict["resname"][0:5]  )
            RESN.append( int(ptclObj.tagsDict["residue"])  )
            ATYPE.append( ptclObj.tagsDict["fftype"]  )
            RING_NUMB.append( int(ptclObj.tagsDict["ring"])  )
            GTYPE.append( ptclObj.tagsDict["gtype"]  )
       
        BONDS = []
        for b_i,bondObj in  self.bondC:
            BONDS.append( [bondObj.pgid1 - 1, bondObj.pgid2 -1])

            print " make_top bonds ",bondObj.pgid1 , bondObj.pgid2

        # Set cubic lattice constant to 5 nm arbitrary 
        LV = np.zeros( (3,3) )
            
        LV[0][0] = self.latvec[0][0]
        LV[1][1] = self.latvec[1][1]
        LV[2][2] = self.latvec[2][2]
        
        # Find atomic number based on atomic symbol 
        ELN = elements.asymb_eln(ASYMB)
        NA = len(ELN)
        
        # Create neighbor list form bonds
        NBLIST,NBINDEX = groups.build_nablist_bonds(ELN,BONDS)
        #NBLIST,NBINDEX = self.bonded_nblist() #groups.build_nablist_bonds(ELN,BONDS)

        print_nb = True
        # Make compatable with 0-(N-1) index of atoms 
        #for n_indx in range(len(NBLIST)):
        #    if( print_nb):
        #        print " changing NBLIST ",NBLIST[n_indx] ," to ",NBLIST[n_indx] -1 
        #    NBLIST[n_indx] =NBLIST[n_indx] -1
        #for n_indx in range(len(NBINDEX)-1):
        #    if( print_nb):
        #        print " changing NBINDEX ",NBINDEX[n_indx] ," to ",NBINDEX[n_indx+1]
        #    NBINDEX[n_indx] =NBINDEX[n_indx+1]
        #
        if( print_nb):

            for p_i in range(len(self.ptclC)):
                N_i_o = NBINDEX[p_i]
                N_i_f = NBINDEX[p_i+1]
                print " atom ",p_i+1, ELN[p_i]," has ",N_i_f - N_i_o
                for indx_j in range( N_i_o,N_i_f):
                    atom_j = NBLIST[indx_j]
                    print "      nb ",atom_j+1," atomic # ", ELN[atom_j]
                    
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)

        for a_indx in range(len(ANGLES)):
            print " angles ",ANGLES[a_indx][0]+1,ANGLES[a_indx][1]+1,ANGLES[a_indx][2]+1
        
        #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)


        #
        # Set charge groups
        #
        verbose = False 
        CG_SET = []
        CHARN = []
        one = 1
        for i in range( len(ELN) ):
            CG_SET.append(one)
            CHARN.append(one)
        CHARN = top.set_chargegroups(verbose,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LV)

        # Read in parameter files 
        
        FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(itp_file)

        # Identify total number of atom types for lammps output 
        ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH)

        # Check atom types to be sure each atom of the same type has the same number of neighbors 
        ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)
    
        ATYPE_EP, ATYPE_SIG = top.atom_parameters(itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
        BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(itp_file,BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
        ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(itp_file,ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
        DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(itp_file, norm_dihparam, DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )

        IMPTYPE_F  = top.imp_parameters(itp_file)

        IMPS,IMPTYPE_F = atom_types.set_ptma_imps(NA,NBLIST, NBINDEX,ELN,ASYMB,IMPS,IMPTYPE_F)


        top_file = dir_id+"/"+output_id + ".top"
        DIH_CONST = []
        DIH_CONST_ANGLE = []
        gromacs.print_top( top_file,ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
               ,DIH_CONST,DIH_CONST_ANGLE
               ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
               ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LV)
        
    def lmp_writedata(self,data_file,norm_dihparam,itp_file):
        """
        Write out lammps data file
        """

        check_bonds = True 

        # Version 1 will be dependent on Atomicpy
        import elements , lammps ,gromacs , atom_types, top , groups

        # New options that need to be passed 
        limdih =  0
        limitdih_n = 1
        
        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        ATYPE = []
        RING_NUMB = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( ptclObj.position)  )
            AMASS.append( float(ptclObj.mass)  )
            CHARGES.append( float(ptclObj.charge)  )
            MOLNUMB.append( int(ptclObj.tagsDict["chain"])  )
            RESID.append( ptclObj.tagsDict["resname"]  )
            RESN.append( int(ptclObj.tagsDict["residue"])  )
            ATYPE.append( ptclObj.tagsDict["fftype"]  )
            RING_NUMB.append( int(ptclObj.tagsDict["ring"])  )

        BONDS = []
        for b_i,bondObj in  self.bondC:
            BONDS.append( [bondObj.pgid1 - 1, bondObj.pgid2 -1])

                    
        # Set cubic lattice constant to 5 nm arbitrary 
        LV = np.zeros( (3,3) )
            
        LV[0][0] = self.latvec[0][0]
        LV[1][1] = self.latvec[1][1]
        LV[2][2] = self.latvec[2][2]
        
        # Find atomic number based on atomic symbol 
        ELN = elements.asymb_eln(ASYMB)
        GTYPE = top.initialize_gtype( ELN )
        NA = len(ELN)
        
        # Create neighbor list form bonds
        # NBLIST,NBINDEX = groups.build_nablist_bonds(ELN,BONDS)
        NBLIST,NBINDEX = self.bonded_nblist() #groups.build_nablist_bonds(ELN,BONDS)

        print_nb = False
        # Make compatable with 0-(N-1) index of atoms 
        for n_indx in range(len(NBLIST)):
            if( print_nb):
                print " changing NBLIST ",NBLIST[n_indx] ," to ",NBLIST[n_indx] -1 
            NBLIST[n_indx] =NBLIST[n_indx] -1
        for n_indx in range(len(NBINDEX)-1):
            if( print_nb):
                print " changing NBINDEX ",NBINDEX[n_indx] ," to ",NBINDEX[n_indx+1]
            NBINDEX[n_indx] =NBINDEX[n_indx+1]

        if( print_nb):

            for p_i in range(len(self.ptclC)):
                N_i_o = NBINDEX[p_i]
                N_i_f = NBINDEX[p_i+1]
                print " atom ",p_i, ELN[p_i]," has ",N_i_f - N_i_o
                for indx_j in range( N_i_o,N_i_f):
                    atom_j = NBLIST[indx_j]
                    print "      nb ",atom_j," atomic # ", ELN[atom_j]
        #
        # Check that the altered neighbor list was done correctly 
        #  
        BONDS_check = top.nblist_bonds(NA,NBLIST, NBINDEX)
        if( len(BONDS_check) != len(BONDS)):
            sys.exit("error in bonds ")
        else:
            for b_i in range(len(BONDS)):                
                if( BONDS[b_i][0] !=  BONDS_check[b_i][0] ):
                    print BONDS[b_i][0], BONDS[b_i][1], BONDS_check[b_i][0], BONDS_check[b_i][1]
                    sys.exit("error in bonds ")
                if( BONDS[b_i][1] !=  BONDS_check[b_i][1] ):
                    print BONDS[b_i][0], BONDS[b_i][1], BONDS_check[b_i][0], BONDS_check[b_i][1]
                    sys.exit("error in bonds ")
                    
            #sys.exit(" checking bonds ")

        
        
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
        #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)

        # Read in parameter files 
        
        FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(itp_file)

        # Identify total number of atom types for lammps output 
        ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH)

        # Check atom types to be sure each atom of the same type has the same number of neighbors 
        ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)
    
        ATYPE_EP, ATYPE_SIG = top.atom_parameters(itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
        BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(itp_file,BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
        ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(itp_file,ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
        DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(itp_file, norm_dihparam, DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )
    
        IMPTYPE_F  = top.imp_parameters(itp_file)

	
	lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
	      BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
	      ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
	      DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
	      RESN,ATYPE_IND,CHARGES,R , ATYPE,
	      BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LV)


    def write_xmol(self, xmol_file,comment,append):
        """
        Write a structure  to an xmol file
        
        Args:
          xmol_file    (str) xmol file name
          comment  (str) for comment line 
          append  (boolean) to append or create a new file 
        Reutrns
          null
        """
        # Open xmol file 
        if(append):
            F = open(xmol_file,"a")
        else:
            F = open(xmol_file,"w")
            
        # Loop over structures
        NP = len( self.ptclC )
        F.write(" %d \n" % NP )
        F.write(" %s \n"%comment)
        for pid, ptclObj  in self.ptclC:
            r_i = ptclObj.position
            atomic_symb = ptclObj.type
            F.write( " %5s %16.8f %16.8f %16.8f \n"  % (atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]) ) )   
        F.close()

    

    def bonded_nblist(self):
        """
        Create neighbor list of bonded particles


        """
        #sys.exit("bonded_nblist not working !!! ")
        
        debug = False 
        NNAB  = 0

        maxnnab = len(self.bondC)*2 + 1
            
        if(debug):
            print " maxnnab",maxnnab
            
        #NBLIST = numpy.empty( maxnnab,  dtype=int )
        #NBINDEX = numpy.empty( maxnnab,  dtype=int )

        # python style nieghbor list
        nblist_py = [] #numpy.empty( maxnnab,  dtype=int )

        NBLIST = []
        NBINDEX = []
        NBLIST.append( 0 )
        
        # First create an N diminsional list of index lists for each particle

        for p_i, prtclC in self.ptclC:
            nblist_py.append( [  ] )
        nblist_py.append( [  ] )

        # bassed on bonds add index of neighbros to particle index of nblist_py
        for b_i,bondObj in  self.bondC:
            bnd_i = bondObj.pgid1 
            bnd_j = bondObj.pgid2 
            
            nblist_py[bnd_i].append( bnd_j )
            nblist_py[bnd_j].append( bnd_i )

            if(debug): print " adding bond bonded_nblist",bnd_i,bnd_j

        # Translate 2D into 1D array
        #   mostly to match perviously writen fortran code
        for p_i in range( len(nblist_py)):
            # loop over  each particle p_i and get list of neighbors nlist_i
            nlist_i = nblist_py[p_i]
            NBINDEX.append( NNAB + 1 )
            # Loop over neighbor list of each particle nlist_i and get neighbor p_j
            for p_j in  nlist_i:

                if(debug): print " p_i ",p_i," p_j ",p_j
                
                #if( p_j > p_i):
                # remove redundent neighbors 
                NNAB +=  1
                # add to neighbor list 
                NBLIST.append( p_j )

        NBINDEX.append( NNAB + 1 )

        if ( debug ):
            print ' total nbs ',NNAB

            for p_i, prtclC in self.ptclC:
                N_i_o = NBINDEX[p_i]
                N_i_f = NBINDEX[p_i+1]
                print " atom ",p_i,prtclC.type, " has ",N_i_f - N_i_o
                					
                for indx_j in range( N_i_o,N_i_f):
                    atom_j = NBLIST[indx_j]
                    print "      nb ",atom_j, self.ptclC[atom_j].type
                
            sys.exit('bonded_nblist debug')
            
        return (NBLIST,NBINDEX)


        
    def calc_rdf(self, rdf_cnt_ij,bin_size,list_i,list_j,sq_r_cut):
        """
        Calculate RDF for a group of particles
        """
        print_vmd = True 
        #
        # Loop over list i
        #
        for p_i, ptcl_i in self.ptclC(list_i):
            r_i = np.array( [float(ptcl_i.position[0]),float(ptcl_i.position[1]),float(ptcl_i.position[2] )] )
            #
            # Loop over list j
            #
            for p_j, ptcl_j in self.ptclC(list_j):
                
                if( p_j > p_i):
                    r_j =  np.array( [float(ptcl_j.position[0]),float(ptcl_j.position[1]),float(ptcl_j.position[2])] )
                    r_ij_sq = pbcs.sq_drij_c(r_i,r_j,self.latvec)
                    if( r_ij_sq <= sq_r_cut ):
                        m_ij = np.sqrt(r_ij_sq)
                        bin_index = int( round( m_ij/bin_size) )
                        rdf_cnt_ij[bin_index] += 2

                        if( print_vmd ):
                            print " index %d or %d "%(p_i,p_j)

        return rdf_cnt_ij

    def get_dihatoms(self,list_k,list_i,list_j,list_l):
        """
        Find sets of dihedrals in system

        k-i-j-l

        Arguments
          list_k (list) of atom indexes in of the first bonded atom in the dihedral 
          list_i (list) of atom indexes in of the second bonded atom in the dihedral 
          list_j (list) of atom indexes in of the third bonded atom in the dihedral 
          list_l (list) of atom indexes in of the fourth bonded atom in the dihedral
        Return
          angle_list (list) of four aotms in each dihedral 
        """

        import datetime

        debug = False

        if(debug):
            print " list_k ",list_k
            print " list_i ",list_i
            print " list_j ",list_j
            print " list_l ",list_l

        
        # Create neighbor list form bonds
        NBLIST,NBINDEX = self.bonded_nblist()

        
        if(debug):
            print " found     NBLIST,NBINDEX "
           
        #
        # Find  atom groups k-i-j
        #
        angle_list = []
        #
        sum_angles = 0
        #
        # Find atom indices  of group i and j
        #
        #for p_k, ptcl_k  in self.ptclC(list_k):
        for atom_k  in list_k:
            if(debug):
               t_i = datetime.datetime.now()

	    N_k_o = NBINDEX[atom_k]
	    N_k_f = NBINDEX[atom_k+1]

            if(debug): print  "checking k ",atom_k,self.ptclC[atom_k].type," with ",N_k_f - N_k_o," nbs"
	    
	    for indx_i in range( N_k_o,N_k_f):
                atom_i = NBLIST[indx_i]
                add_i = False

                if( debug ): print " checking i ",atom_i,self.ptclC[atom_i].type
                
                for  p_i in list_i:
                    if( atom_i == p_i ): #and atom_i > atom_k ):
                        add_i = True
                if( add_i ): #atom_i in list_i ):

                    
		    N_i_o = NBINDEX[atom_i]
		    N_i_f = NBINDEX[atom_i+1] 
					
		    for indx_j in range( N_i_o,N_i_f):
			atom_j = NBLIST[indx_j]
			add_j = False

                        if( debug ): print " checking j ",atom_j,self.ptclC[atom_j].type
                            
                        for  p_j  in list_j:
                            if( atom_j == p_j ): #and atom_j > atom_i ):
                                add_j = True
                        if( add_j ): # atom_j  in list_j  ):

                            N_j_o = NBINDEX[atom_j]
                            N_j_f = NBINDEX[atom_j+1] 

                            for indx_l in range( N_j_o,N_j_f):
                                atom_l = NBLIST[indx_l]
                                add_l = False

                                if( debug ):
                                    print " checking l ",atom_l,self.ptclC[atom_l].type                                    
                                for  p_l  in list_l:
                                    if( atom_l == p_l ): #and atom_l > atom_j ):
                                        add_l = True
                                if( add_l ): #atom_l  in list_l ):
                                    # Check to make sure not in list already
                                    add_dih = True
                                    for indx_kij in angle_list:
                                        a_k = indx_kij[0]                
                                        a_i = indx_kij[1]
                                        a_j = indx_kij[2]
                                        a_l = indx_kij[3]
                                        if( atom_k == a_k and atom_i == a_i and atom_j == a_j and atom_l == a_l ):
                                            add_dih = False 
                                        if( atom_k == a_l and atom_i == a_j and atom_j == a_i and atom_l == a_k ):
                                            add_dih = False 
                                    if(add_dih ):
                                        angle_list.append( [atom_k,atom_i,atom_j,atom_l] )

                                        if(debug):

                                            t_f = datetime.datetime.now()
                                            dt_sec  = t_f.second - t_i.second
                                            dt_min  = t_f.minute - t_i.minute
                                            if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                            if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0         
                                            print " found angle ",atom_k,atom_i,atom_j,atom_l
                                            print "  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "
                                        

        if ( debug ):
            dih_cnt = 0 
            for indx_kij in angle_list:
                dih_cnt += 1 
                a_k = indx_kij[0]                
                a_i = indx_kij[1]
                a_j = indx_kij[2]
                a_l = indx_kij[3]
                print dih_cnt," found angle ",a_k,a_i,a_j,a_l
                                    
            sys.exit('get_dihatoms debug')
            
	return angle_list

    def getDihedral(self,a_k,a_i,a_j,a_l):
        """
        Calculate dihedral angle of set of particles

        k
           \
             i  - j
                    \
                       l

                       
        
        """
        debug = False
        
        r_k = np.array( [ float( self.ptclC[a_k].position[0] ),float( self.ptclC[a_k].position[1] ),float( self.ptclC[a_k].position[2] ) ] )
        r_i = np.array( [ float( self.ptclC[a_i].position[0] ),float( self.ptclC[a_i].position[1] ),float( self.ptclC[a_i].position[2] ) ] )
        r_j = np.array( [ float( self.ptclC[a_j].position[0] ),float( self.ptclC[a_j].position[1] ),float( self.ptclC[a_j].position[2] ) ] )
        r_l = np.array( [ float( self.ptclC[a_l].position[0] ),float( self.ptclC[a_l].position[1] ),float( self.ptclC[a_l].position[2] ) ] )

        if( debug ):
            print " lat vec ",self.latvec
            print "r_k ",r_k
            print "r_i ",r_i
            print "r_j ",r_j
            print "r_l ",r_l

        v1 = pbcs.norm_r_ij(r_k, r_i,self.latvec)
        v2 = pbcs.norm_r_ij(r_i, r_j,self.latvec)
        v3 = pbcs.norm_r_ij(r_j, r_l,self.latvec)

        v1v2 = np.cross(v1,v2)
        v2v3 = np.cross(v2,v3)


        angle_i = pbcs.getAngle(v1v2,v2v3)

        
        if( debug):
            print "v1 ",v1
            print "v2 ",v2
            print "v3 ",v3
            print "v1v2 ",v1v2
            print "v2v3 ",v2v3
            print "angle_i ",angle_i


        #
        # Find sign of angle 
        #
        v1v3 = np.cross(v1,v3)
        sign_v = np.dot(v2,v1v3)

        if( sign_v > 0.0  ):
            angle_i = -1.0*angle_i

        return angle_i

    def get_gromacs(self, gro_file,top_file):  # Move out of class
        """
        Read in structure information from gromacs files

        Args:
          gro_file (str) gro file
          top_file (str) top file
      
          
        """
        import gromacs , elements

        GTYPE,R,VEL,LV = gromacs.read_gro(gro_file)
        ATYPE,RESN,RESID,GTYPE,CHARN,CHARGES,AMASS,BONDS,ANGLES,DIH,MOLNUMB,MOLPNT,MOLLIST = gromacs.read_top(top_file)
        ASYMB,ELN  = elements.mass_asymb(AMASS)
        

        for p_i in range( len( ASYMB)):
            r_i =  [ float(R[p_i][0]), float(R[p_i][1]), float(R[p_i][2]) ] 
            atomic_symb = str( ASYMB[p_i] )
            m_i = float(AMASS[p_i])
            q_i = float(CHARGES[p_i])
            # Create particle
            pt_i = Particle( r_i,atomic_symb,q_i,m_i )
            # Find needed tags
            chain_i = int( MOLNUMB[p_i] )
            ring_i = 0 #particle_data["ring"][p_i]
            resname_i = RESID[p_i]
            residue_i = RESN[p_i]
            linkid_i = "UNKNOWN" #particle_data["linkid"][p_i]
            fftype_i =ATYPE[p_i]
            gtype_i = GTYPE[p_i]
            # _i = particle_data[""][p_i]
            # Add particle to structure 
            tagsD = {"chain":chain_i,"ring":ring_i,"resname":resname_i,"residue":residue_i,"linkid":linkid_i,"fftype":fftype_i,"gtype":gtype_i}
            pt_i.setTagsDict(tagsD)
            self.ptclC.put(pt_i)

        # Read in bonds
        for b_indx in range( len(BONDS )):
            a_i = int(BONDS[b_indx][0] )
            a_j = int(BONDS[b_indx][1] )
            b_i = Bond( a_i + 1 , a_j + 1  )            
            self.bondC.put(b_i)

        # Read in lattice vectors
        self.latvec = []
        self.latvec.append(  np.array( [float(LV[0][0]),float(LV[0][1]),float(LV[0][2])] ) )
        self.latvec.append(  np.array( [float(LV[1][0]),float(LV[1][1]),float(LV[1][2])] ) )
        self.latvec.append(  np.array( [float(LV[2][0]),float(LV[2][1]),float(LV[2][2])] ) )


    def get_gaussian(self,fchk_file,json_data):
        """
        Read in structure information from gaussian fchk file

        Args:
          fchk_file (str) grochkfile
        
        """
        import numpy as np
        # Energy conversion
        # http://physics.nist.gov/cgi-bin/cuu/Value?threv
        HtoeV = 27.211385

        bohr2angstrom = 0.5291772086

        
        particle_data = json_data["structure"]["particle"]
        
        F = open(fchk_file,'r')
        Lines = F.readlines()
        F.close()

        read_r = False
        read_eln = False
        read_esp = False
        for line in Lines :
            col = line.split()

            if( read_r ):

                if (  col[0] == "Force" and col[1] == "Field" ):
                    read_r = False
                    p_i = 0 
                    for p_indx in range(NA):
                        #print atom_i ,atom_i*3,atom_i*3+2,R_all[atom_i*3:atom_i*3+3]
                        vec_r_i =  R_all[p_indx*3:p_indx*3+3]
                        p_i += 1                         
                        self.ptclC[p_i].position = vec_r_i
                        particle_data["position"][p_indx] = vec_r_i
                else:
                    for r_i in  map(float,col) :
                        R_all.append( r_i*bohr2angstrom )

            if( read_eln ):
                if ( eln_p_cnt == NA ):
                    read_eln = False
                else:                    
                    for eln_i in  map(int,col):
                        eln_p_cnt += 1                        
                        if( eln_i != self.ptclC[eln_p_cnt].type ):
                            print "  Particle ",eln_p_cnt, self.ptclC[eln_p_cnt].type," != ",eln_i

            if( read_esp ):
                if ( esp_p_cnt == NA ):
                    read_esp = False
                else:
                    
                    for q_i in  map(float,col):
                        esp_p_cnt += 1
                        p_indx = esp_p_cnt -1 
                        self.ptclC[esp_p_cnt].charge = q_i
                        particle_data["charge"][p_indx] = q_i
                                                

            if( len(col) > 2 ):
                if( col[0] == "Total" and col[1] == "Energy" ):
                    TOTAL_ENERGY = float( col[3] )*HtoeV

            if( len(col) == 5 ):
                if( col[0] == "Number" and col[1] == "of"  and col[2] == "atoms" ):
                    NA = int(col[4])
                    
                    if( NA != len(self.ptclC) ):
                        print " json file contains %d atoms and fchk file contains %d "%(self.ptclC,NA)
                        sys.exit("inconsistent files ")

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

        return json_data


    def set_cply_tags(self):  # Move outside of class
        """
        Set tags for new cply file 
        Use ctype tag and bonding enviroment
        """

        # Version 1 will be dependent on Atomicpy
        import elements , top 

        # Create list to pass to Atomicpy
        ASYMB = []
        R = []
        AMASS = []
        CHARGES = []
        MOLNUMB = []
        RESID = []
        RESN = []
        CTYPE = []
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( [ float( ptclObj.position[0] ),float( ptclObj.position[1] ),float( ptclObj.position[2] )]  )  )
            AMASS.append( float( ptclObj.mass)  )
            CHARGES.append( float( ptclObj.charge ) )
            MOLNUMB.append( int( ptclObj.tagsDict["chain"] ) )
            RESID.append( int( ptclObj.tagsDict["residue"] ) )
            RESN.append(  ptclObj.tagsDict["resname"]  )
            CTYPE.append( ptclObj.tagsDict["linkid"]  )
            
        # Find atomic number based on atomic symbol 
        ELN = elements.asymb_eln(ASYMB)
        AMASS = elements.eln_amass(ELN)

        #   Build covalent nieghbor list for bonded information 
        NBLIST, NBINDEX = top.build_covnablist(ELN,R)

        verbose = True 
        cply_tag = top.set_cply_tags(  verbose, ELN, CTYPE,RESN ,NBLIST, NBINDEX )
        
        atom_i = 0
        for pid, ptclObj  in self.ptclC:
            ptclObj.tagsDict["cply_tag"] = cply_tag[atom_i]

            print " setting cply tag ",cply_tag[atom_i]
            
            atom_i +=1 

    def replicate(self,p,options): # Move outside of class

        """
        Replicate structures

        Arguments
            p (object) mpirNREL 

            options
                verbose
                ptime
                output_id
                dir_id
                json
                top
                gro
                calc_overlap 

                sol_json
                sol_top
                sol_gro

                sol_buf
                atomic_cut
                den_target
                atoms_target
                max_mol_place
                max_sys
                lc_expand
                perc_sol
                itp_file
                norm_dihparam


        Returns: None

        """
        import mpiNREL
        import file_io
        from structureContainer import StructureContainer
        import sys , datetime, random, math 

        debug = False 
        const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

        # If options set for fixed then set random seed (for regression testing)
        if options.fixed_rnd_seed:
            random.seed(0)
        
        # MPI setup
        rank = p.getRank()
        size = p.getCommSize()

        # Open log file 
        if( rank == 0  ):
            log_file = options.output_id + ".log"
            log_out = open(log_file,"w")

            log_lines = " Running on %d processors  \n"%(size)
            print log_lines
            log_out.write(log_lines)
            
            t_i = datetime.datetime.now()

        # Read in oligomers
        #   from json files
        oligo_array = []
        if( len(options.json) > 0 ):
            oligo_array = file_io.struc_array_json(oligo_array,options.json)
        #    from gromacs  files
        if( len(options.gro) > 0 ):
            oligo_array = file_io.struc_array_gromacs(oligo_array,options.gro,options.top)

        # Read in solvents
        #   from json files 
        sol_array = []
        if( len(options.sol_json) > 0 ):
            sol_array = file_io.struc_array_json(sol_array,options.sol_json)
        #   from gromacs  files
        if( len(options.gro) > 0 ):
            sol_array = file_io.struc_array_gromacs(sol_array,options.sol_gro,options.sol_top)

        # Set location of origin
        org = np.array( [0.0,0.0,0.0] )

        oligo_cnt     = 0      # Number of oligomers to replicate 
        oligo_nprt    = 0      # Total number of particles in all the oligomers 
        oligo_mass    = 0.0    # Total mass of particles in all the oligomers
        oligo_maxlength = 0.0 
        for oligo_i in oligo_array:
            oligo_cnt += 1
            nprt = len( oligo_i.ptclC )

            print "nprt",nprt

            oligo_nprt += nprt 
            tot_mass = oligo_i.getTotMass()
            oligo_mass += tot_mass
            nbonds = len( oligo_i.bondC )
            # Shift center of mass to origin
            oligo_i.shift_center_mass(org)

            # Record oligomer molecule length
            oligo_length = oligo_i.getlength()

            # Record solvent molecule length for grid spacing 
            oligo_length = oligo_i.getlength()
            if( oligo_length > oligo_maxlength):
                oligo_maxlength = oligo_length

            if( options.verbose and rank == 0  ):
                log_lines = ""
                log_lines += "  Oligomers %d \n"%oligo_cnt 
                log_lines +=  "    Particles  %d \n"%nprt
                log_lines +=  "    Total mass   %f \n"%tot_mass
                log_lines +=  "    Length   %f \n"%oligo_length
                log_lines +=  "    Bonds  %d \n"%nbonds
                print log_lines
                log_out.write(log_lines)


        sol_cnt     = 0   # Number of solvent to replicate 
        sol_nprt    = 0   # Total number of solvent in all the solture 
        sol_mass    = 0.0    # Total mass of solvent in all the solture
        sol_maxlength = -100000.0 
        #

        for sol_i in sol_array:
            sol_cnt += 1
            nprt = len( sol_i.ptclC )
            sol_nprt += nprt 
            tot_mass = sol_i.getTotMass()
            sol_mass += tot_mass
            nbonds = len( sol_i.bondC )

            # Shift center of mass to origin
            sol_i.shift_center_mass(org)

            # Record solvent molecule length for grid spacing 
            sol_length = sol_i.getlength()
            if( sol_length > sol_maxlength):
                sol_maxlength = sol_length

            if( options.verbose and rank == 0  ):
                log_lines = ""
                log_lines += "  Solvents %d \n"%sol_cnt 
                log_lines +=  "    Particles  %d \n"%nprt
                log_lines +=  "    Total mass   %f \n"%tot_mass
                log_lines +=  "    Length   %f \n"%sol_length
                log_lines +=  "    Bonds  %d \n"%nbonds
                print log_lines
                log_out.write(log_lines)

        #
        # Calculate the number of oligomers and  solvent molecules
        #
        if( options.perc_sol > 0.0 ):
            # Variables
            #   atoms_target - target number of atoms # options.atoms_target
            #   perc_sol - perecent solvent by mass
            #   frac_sol - fraction solvent by mass
            #   sol_nprt - number atoms in the list of solvent molecules 
            #   oligo_nprt - number atoms in the list of oligomer molecules 
            #   n_sol_l - number of solvent list replications
            #   n_olgio_l - number of oligomer list replications
            #   sol_mass - mass of all the solvents in the solvent list
            #   oligo_mass - mass of all the oligomers in the oligomer list
            # Equations
            #   Equ 1 : perc_sol = n_sol_l*sol_mass/( n_sol_l*sol_mass + n_olgio_l*oligo_mass )
            #   Equ 2 : atoms_target =  n_sol_l*sol_nprt + n_olgio_l*oligo_nprt
            # Solutions
            frac_sol = options.perc_sol/100.0
            n_sol_l= int((-1.0*frac_sol*oligo_mass*float(options.atoms_target))/(frac_sol*sol_mass*float(oligo_nprt)  - frac_sol*oligo_mass*float(sol_nprt) - sol_mass*float(oligo_nprt) ))
            n_olgio_l = int( (float(options.atoms_target) - float(sol_nprt)*float(n_sol_l))/float(oligo_nprt))
        else:
            n_sol_l = 0
            n_olgio_l = int(float(options.atoms_target)/float(oligo_nprt))
            solv_box_l = 0.0


            print " atoms_target oligo_bnprt ",float(options.atoms_target),float(oligo_nprt)

        #
        # Calculate the box size for a target density 
        #
        target_density_amuang = options.den_target*const_avo/10.0 # densit in AMU/Angstrom^3


        total_n = n_olgio_l*oligo_nprt + n_sol_l*sol_nprt
        total_mass = oligo_mass*n_olgio_l + n_sol_l*sol_mass
        volume_target_ang = total_mass/target_density_amuang
        len_target_ang = volume_target_ang**(1.0/3.0)

        if( options.perc_sol > 0.0 ):
            # Check to be sure the grid is large enough
            sol_box_side = int(math.ceil(n_sol_l**(1.0/3.0) ) )         # Number of solvents per box side 
            sol_length = sol_maxlength + options.sol_buf           # length of solvent 
            sol_length_sq = sol_length*sol_length                  # length of solvent squared for overlap calculation 
            vol_olgio = float(n_olgio_l)*(oligo_maxlength**3.0)    # Volume occupied by oligomer 
            len_sol_box =  float(sol_box_side)*sol_length              # Length of pure solvent box
            vol_sol = len_sol_box**3.0                                 # Volume occupied by solvent 
            tot_mol_vol = vol_olgio + vol_sol                      # total volume for oligomers and solvents
            len_tot_box = tot_mol_vol**(1.0/3.0)
                
            if( len_target_ang < len_tot_box ):
                if( options.verbose and rank == 0 ):
                    log_lines = "  Expanding box length from %f to %f based on the length of the solvent  "%(len_target_ang,len_tot_box)
                    print log_lines
                    log_out.write(log_lines)
                    
                len_target_ang = len_tot_box
                volume_target_ang = len_tot_box**3.0

            sol_grid_len = len_target_ang/float(sol_box_side)
            
            print "len_target_ang ",len_target_ang 
            print "sol_grid_len ",sol_grid_len 
            
        else:
            len_target_ang = volume_target_ang**(1.0/3.0)

        cut_ij_sq = options.atomic_cut* options.atomic_cut
        # Calculate actual final structure properties
        vol_f = len_target_ang**3.0
        den_AMU_f = total_mass/vol_f
        den_f = den_AMU_f/const_avo*10.0
        perc_sol_f = (n_sol_l*sol_mass)/(n_sol_l*sol_mass + n_olgio_l*oligo_mass)

        # Recalculate solvent molecules along box length 
        
        print "target_density_amuang",target_density_amuang, options.den_target
        print "total_n",total_n
        print "total_mass",total_mass
        print "len_target_ang",len_target_ang
        print n_sol_l,sol_mass,n_olgio_l,oligo_mass

        # oligomer molecule container
        oligomer_rep = StructureContainer() 
        # Set lattice vector to new box size
        latvec_list = [np.array([len_target_ang,0.0,0.0]),np.array( [0.0,len_target_ang,0.0]),np.array( [0.0,0.0,len_target_ang]) ]
        print "latvec_list" , latvec_list
        oligomer_rep.setLatVec(latvec_list)

        print " oligomer_rep.getLatVec() ",oligomer_rep.getLatVec()

                                

        # Solvent molecule container
        sol_rep = StructureContainer()  
        sol_rep.setLatVec(latvec_list)
        self.setLatVec(latvec_list)

        print " s1 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s1 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s1 self. .getLatVec() ",self.getLatVec()
        

        # Print script information and settings 
        if( rank == 0  ):

            log_lines =  "  Replication settings \n"
            log_lines += "   - Tragets \n"
            log_lines += "       Total atoms %d \n"%(options.atoms_target)
            log_lines += "       Density %f g/cm^3 %f AMU/Angstrom^3 \n"%(options.den_target,target_density_amuang)
            log_lines += "       Solvent mass percentage %f \n"%(options.perc_sol)
            log_lines += "   - Oligomers \n"
            log_lines += "       Total atoms in set %d \n"%(oligo_nprt)
            log_lines += "       Set of structures will replicated %d times \n"%(n_olgio_l)
            if( options.perc_sol > 0.0 ):
                log_lines += "   - Solvents \n"
                log_lines += "       Total atoms in set %d \n"%(sol_nprt)
                log_lines += "       With a length of  %f Angstrom \n"%(sol_length)
                log_lines += "       Set of structures will replicated %d times \n"%(n_sol_l)
                log_lines += "       With an initial grid spacing of %f  \n"%(sol_grid_len)

            log_lines += "   - Final porperties  \n"
            log_lines += "       Total atoms %d \n"%(total_n)
            log_lines += "       Volume %f Angstrom^3 \n"%(volume_target_ang)
            log_lines += "       Density %f g/cm^3 %f AMU/Angstrom^3 \n"%(den_f,den_AMU_f)
            log_lines += "       Solvent mass percentage %f \n"%(perc_sol_f)
            log_lines += "       Cubic cell with length of %f Angstroms \n"%(len_target_ang)
            log_lines += "   - Placement \n"
            log_lines += "       Maximum structure placements %d \n"%(options.max_mol_place)
            log_lines += "       Maximum number of restarts before box is expanded %d \n"%(options.max_sys)
            log_lines += "       Percent of box size to add during expantion %8.2f \n"%(100*options.lc_expand)
            log_lines += "   - id's \n"
            log_lines += "       Directory %s \n"%(options.dir_id)	
            log_lines += "       Output id %s \n"%(options.output_id)
            print log_lines
            log_out.write(log_lines)

        p.barrier()
        


        #sys.exit(" debug 2 ")
        
        # Record initial time
        if( rank == 0  ): 
            t_i = datetime.datetime.now()

        sys_oligo_n = 0     # Number of oligomers add to the system
        sys_attempts = 0    # Number of times the system has been reset
        struc_add_cnt = 0   # Total number of structures added to the final structure 
        #
        # Start adding molecules to the system
        #
        add_oligo = True

        print " Adding %d  oligomers  "%n_olgio_l
        while ( add_oligo ):
            #
            # Initialize 
            #
            add_oligo = True
            overlap_found = True
            strucadd_atempts = 0


            # Record intial time for pereformance testing 
            if( rank == 0  ): 
                tadd_i = datetime.datetime.now()
            for oligo_l in range( n_olgio_l ):
                # loop over the number of times each oligomer in the oligomer list needs to be replicated
                for struc_i in oligo_array:


                    print " s12 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
                    print " s12 sol_rep.getLatVec() ",sol_rep.getLatVec()
                    print " s12 self. .getLatVec() ",self.getLatVec()

        
                    #
                    # Place the atomic indices into list 
                    # 
                    pointIndices = range( len(struc_i.ptclC)  )
                    if( debug ):
                        print rank, size," splitOnProcs "
                    # Create a list of atomic indices for each processor 
                    myChunk  = p.splitListOnProcs(pointIndices)
                    p_debug = False 
                    if(p_debug):                
                        print " cpu ",rank ," has atoms ",myChunk[0]," - ",myChunk[len(myChunk)-1],"  \n"

                        sys.exit("P debug 1 ")
                    
                    # For each structure add to 
                    while ( overlap_found ):
                        strucadd_atempts += 1

                        n_dim = 3
                        rot_angle_i_o = 0.0 
                        rot_angle_j_o = 0.0 
                        r_random_o  = np.zeros(n_dim)
                        
                        if ( rank == 0 ):
                            #
                            #  Get random rotation angles from single processor  
                            #
                            ang_acc = 1000  # number of digets in random angle 
                            rot_angle_i_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 1
                            rot_angle_j_o = float(random.randrange(0,ang_acc))*np.pi/float(ang_acc) # Random angle 2
                            #
                            #  Get random translation from single processor 
                            #
                            r_random_o = np.zeros(n_dim)
                            for x_indx in range( n_dim ):
                                r_random_o[x_indx] = random.randrange(0,ang_acc*int(oligomer_rep.latvec[x_indx][x_indx]) )/float(ang_acc)

                            debug = 0
                            if( debug ):
                                print " ran ",x_indx,(oligomer_rep.latvec[x_indx][x_indx])
                                print rot_angle_i_o,rot_angle_j_o,r_random_o
                                #sys.exit(" Random # 's test 1")
                                
                        p.barrier() # Barrier for MPI_COMM_WORLD
                        #
                        # Broadcast random rotation angles and translations to all processors 
                        #
                        rot_angle_i = p.bcast(rot_angle_i_o)
                        rot_angle_j = p.bcast(rot_angle_j_o)
                        r_random = p.bcast(r_random_o)
                        p.barrier() # Barrier for MPI_COMM_WORLD
                        #
                        # Get coordinates of randomly rotated and shifted 
                        #
                        struc_i.shift_center_mass(org)
                        struc_i.rotate(rot_angle_i,rot_angle_j)
                        struc_i.vec_shift(r_random)

                        overlap = 0
                        if( len(oligomer_rep.ptclC) > 0 ):
                            #
                            # If there are particles in the system check atoms do not overlap
                            #                            
                            if( options.calc_overlap ):

                                for p_i, ptclObj_i in struc_i.ptclC(myChunk):
                                    r_i = np.array( ptclObj_i.position )
                                    for p_sys, ptclObj_sys in oligomer_rep.ptclC :
                                        r_sys = np.array( ptclObj_sys.position )
                                        r_ij_sq = pbcs.sq_drij_c(r_i,r_sys,oligomer_rep.getLatVec() )
                                        if( r_ij_sq < cut_ij_sq ):
                                            overlap = 1

                        p.barrier() # Barrier for MPI_COMM_WORLD
                        #
                        # Reduce sum the overlap variable from all the processors
                        #   if it is zero everywhere there was no overlap detected 
                        #
                        overlap_sum = p.allReduceSum(overlap)
                        p.barrier() # Barrier for MPI_COMM_WORLD
                        
                        print " s13 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
                        print " s13 sol_rep.getLatVec() ",sol_rep.getLatVec()
                        print " s13 self. .getLatVec() ",self.getLatVec()

        
                        if( overlap_sum ==  0 ):
                            # If no overlap detected add molecule to the system 
                            sys_oligo_n += 1
                            struc_add_cnt += 1
                            # Rest molecule numbers
                            for pid, ptclObj in struc_i.ptclC :
                                ptclObj.tagsDict["chain"] = struc_add_cnt

                            # add oligomer structure to system structure
                            struc_i.setLatVec(oligomer_rep.getLatVec())
                            oligomer_rep += struc_i

                            if( options.verbose ):
                                if( rank == 0  ):
                                    print "      -  Molecule ",sys_oligo_n," has been added to the system after ",strucadd_atempts," placment attempts "
                                    print "         system has %d atoms and %d bonds "%(len(oligomer_rep.ptclC),len(oligomer_rep.bondC))
                                    #print " Printing  oligomer_rep bonds "
                                    #oligomer_rep.printbondlengths()


                                    print " s11 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
                                    print " s11 sol_rep.getLatVec() ",sol_rep.getLatVec()
                                    print " s11 self. .getLatVec() ",self.getLatVec()


                                    if( options.ptime ):
                                        t_f = datetime.datetime.now()
                                        dt_sec  = t_f.second - t_i.second
                                        dt_min  = t_f.minute - t_i.minute
                                        if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                        print "        - with placement time ",dt_min," min ",dt_sec," seconds "

                            overlap_found = False
                        else:
                            overlap_found = True

                        if( strucadd_atempts >= options.max_mol_place ):
                            # If attempts to place molecule into the system exceed max set by max_mol_place

                            #   reset system and star over 
                            if(  rank == 0  ):
                                if(  options.verbose ):


                                    print " add0 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()


                                    print "        -  Attempts to add molecule ",sys_oligo_n+1," has exceeded max attempts ",options.max_mol_place," system will be reset for the ",sys_attempts," time "

                            sys_oligo_n = 0
                            struc_add_cnt = 0 
                            strucadd_atempts = 0
                            sys_attempts += 1

                            # Save lattice vectors as to no loose any expansions 
                            latvec_i = oligomer_rep.getLatVec()

                            print " saving s1 latvec_i ",latvec_i
                            
                            # Delete system 
                            del oligomer_rep
                            oligomer_rep = StructureContainer()  # Output replicated structure
                            # Set lattice vectors 
                            oligomer_rep.setLatVec(latvec_i) 


                        if( sys_attempts >= options.max_sys  ):


                            print " exp0 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()


                            # If the system has been reset over max_sys times expand the box size by lc_expand
                            oligomer_rep.expandLatVec(options.lc_expand)

                            # Save lattice vectors as to no loose any expansions 
                            latvec_i = oligomer_rep.getLatVec()


                            print " exp1 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
                            
                            # Delete system 
                            del oligomer_rep
                            oligomer_rep = StructureContainer()  # Output replicated structure
                            # Set lattice vectors 
                            oligomer_rep.setLatVec(latvec_i) 

                            sys_attempts = 0

                            if( options.verbose ):                
                                if( rank == 0  ):
                                    print '          - Number of system resets has exceeded the maximum  (option max_sys) ',options.max_sys
                                    print '          - Lattice vectors will be expanded by (option lc_expand)',options.lc_expand
                                    print '             v_1 ',latvec_i[0]
                                    print '             v_2 ',latvec_i[1]
                                    print '             v_3 ',latvec_i[2]

                p.barrier() # Barrier for MPI_COMM_WORLD


            if( sys_oligo_n ==  n_olgio_l  ):
                # If all the molecule have been added exit while loop and print system 
                add_oligo = False
                latvec_oligo = oligomer_rep.getLatVec()
                p.barrier() # Barrier for MPI_COMM_WORLD
                if( options.verbose and rank == 0  ):
                    print " All oligomers  have been added "

        # Add replicated oligmers to final structure
        # self = StructureContainer()

        debgu_n = False
        if( debgu_n ):

            print "self pre add "
            print self.ptclC
            print self.bondC

            print "oligomer_rep"
            print oligomer_rep.ptclC
            print oligomer_rep.bondC


        print " s21 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s21 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s21 self. .getLatVec() ",self.getLatVec()
                

        self =  oligomer_rep
        self.setLatVec(latvec_oligo)
        # Rest solvent lattice vectors to match oligo 
        sol_rep.setLatVec(latvec_oligo)

        if( debgu_n ):
            print "self"
            print self.ptclC
            print self.bondC

            sys.exit("debug ")


        print " s2 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s2 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s2 self. .getLatVec() ",self.getLatVec()
                

        if( options.perc_sol > 0.0 ):

            sys_sol_n = 0    # Number of solvent add to the system
            sys_attempts = 0  # Number of times the system has been reset
            #
            # Start adding molecules to the system
            #
            print " Adding %d  solvent  "%n_sol_l
            add_sol = True
            while ( add_sol ):
                #
                # Initialize 
                #
                strucadd_atempts = 0

                # Record intial time for pereformance testing 
                if( rank == 0  ): 
                    tadd_i = datetime.datetime.now()

                bad_grid_point = False 
                for sol_l in range( n_sol_l ):

                    #
                    # 
                    #
                    print " starting %d of %d "%(sol_l,n_sol_l)

                    # loop over the number of times each solvent in the solvent list needs to be replicated
                    if( bad_grid_point ): break
                    if( not add_sol): break 
                    for struc_i in sol_array:
                        overlap_found = True
                        if( bad_grid_point ): break
                        if( not add_sol): break 
                        while ( overlap_found ):
                            strucadd_atempts += 1
                            for x_indx in range(sol_box_side):
                                if( not add_sol): break 
                                for y_indx in range(sol_box_side):
                                    if( not add_sol): break 
                                    for z_indx in range(sol_box_side):
                                        
                                        if( sys_sol_n == n_sol_l*sol_cnt ):
                                            add_sol = False
                                            break
                                            
                                
                                        #l_x =  float(lat_indx[0])*sol_length
                                        #l_y =  float(lat_indx[1])*sol_length
                                        #l_z =  float(lat_indx[2])*sol_length

                                        l_x =  float(x_indx)*sol_grid_len
                                        l_y =  float(y_indx)*sol_grid_len
                                        l_z =  float(z_indx)*sol_grid_len

                                        print " Checking overlap for solvent %d at lattice point %f %f %f "%(sys_sol_n,l_x,l_y,l_z)

                                        if( l_x > oligomer_rep.latvec[0][0] or l_y > oligomer_rep.latvec[1][1] or l_z > oligomer_rep.latvec[2][2] ):
                                            print " Lattic point beyond box %f %f %f "%(oligomer_rep.latvec[0][0],oligomer_rep.latvec[1][1], oligomer_rep.latvec[2][2])
                                            bad_grid_point = True
                                            break 

                                        lat_pos = np.array( [l_x,l_y,l_z] )

                                        # Make sure there is no overlap with the added molecules
                                        overlap = 0
                                        for p_i, ptclObj_i in oligomer_rep.ptclC :
                                            r_i = np.array( ptclObj_i.position )
                                            r_ij_sq = pbcs.sq_drij_c(r_i,lat_pos,oligomer_rep.getLatVec() )
                                            if( r_ij_sq < sol_length_sq ):
                                                overlap = 1

                                        p.barrier() # Barrier for MPI_COMM_WORLD
                                        #
                                        # Reduce sum the overlap variable from all the processors
                                        #   if it is zero everywhere there was no overlap detected 
                                        #
                                        overlap_sum = p.allReduceSum(overlap)

                                        if( overlap_sum ==  0 ):
                                            # If no overlap detected add molecule to the system 
                                            sys_sol_n += 1
                                            overlap_found = False 
                                            # Shift molecule to lattice point
                                            struc_i.shift_center_mass(org)
                                            struc_i.vec_shift(lat_pos)

                                            struc_i.setLatVec(sol_rep.getLatVec())
                                            sol_rep += struc_i

                                            if( options.verbose ):
                                                if( rank == 0  ):
                                                    print "      -  Molecule %d  has been added to the system at lattice point %f %f %f  "%(sys_sol_n,l_x,l_y,l_z)

                                                    if( options.ptime ):
                                                        t_f = datetime.datetime.now()
                                                        dt_sec  = t_f.second - t_i.second
                                                        dt_min  = t_f.minute - t_i.minute
                                                        if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
                                                        if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0
                                                        print "        - with placement time ",dt_min," min ",dt_sec," seconds "
                                        else:
                                            overlap_found = True

                                            print "  lattice point was found to overlap "



                if( sys_sol_n == n_sol_l*sol_cnt ):
                    add_sol = False
                    p.barrier() # Barrier for MPI_COMM_WORLD
                    if( options.verbose and rank == 0  ):
                        print " All solvents  have been added "

                else:
                    #
                    # If attempts to place solvent molecule into the system failed
                    #
                    sys_sol_n = 0                
                    sys_attempts += 1
                    #
                    # Save lattice vectors as to not loose any expansions
                    #
                    latvec_i = sol_rep.getLatVec()
                    #
                    # Delete system
                    #
                    del sol_rep
                    sol_rep = StructureContainer()  # Output replicated structure
                    #
                    # Set lattice vectors
                    #
                    sol_rep.setLatVec(latvec_i)
                    #
                    #
                    if( (sol_grid_len*0.90) > sol_length):
                        # If grid spacing is larger than the solvent length shrink grid
                        sol_grid_len = sol_grid_len*0.90
                    else:
                        # Otherwise increase volume and set grid spacing to solvent length
                        sol_grid_len = sol_length
                        #
                        # Expand the box size by a single solvent length 
                        #
                        #sol_rep.expandLatVec(options.lc_expand)
                        sol_rep.latvec[0] = sol_rep.latvec[0] + sol_length
                        sol_rep.latvec[1] = sol_rep.latvec[1] + sol_length
                        sol_rep.latvec[2] = sol_rep.latvec[2] + sol_length
                        #
                        # Increase the number of solvent molecules along the box by 1
                        #
                        sol_box_side += 1


                    if( options.verbose ):                
                        if( rank == 0  ):
                            print '          - Lattice vectors will be expanded by solvent length %f ',sol_length

                p.barrier() # Barrier for MPI_COMM_WORLD

            # Add replicated solvents to final structure
            self += sol_rep
            self.setLatVec(sol_rep.latvec)

        

        print " s3 oligomer_rep. .getLatVec() ",oligomer_rep.getLatVec()
        print " s3 sol_rep.getLatVec() ",sol_rep.getLatVec()
        print " s3 self. .getLatVec() ",self.getLatVec()
                

        self.compressPtclIDs()
        
        print "         f_rep has %d atoms and %d bonds "%(len(self.ptclC),len(self.bondC))
        print "             lat vec ",sol_rep.latvec
        #self.printbondlengths()


        if( rank == 0 ):


            t_f = datetime.datetime.now()
            dt_sec  = t_f.second - t_i.second
            dt_min  = t_f.minute - t_i.minute
            if ( dt_sec < 0 ): dt_sec = 60.0 - dt_sec
            if ( dt_sec > 60.0 ): dt_sec = dt_sec - 60.0 
            log_line="\n  Finished time  " + str(t_f)
            log_out.write(log_line)
            log_line="\n  Computation time "+str(dt_min) + " min "+str(dt_sec)+" seconds "
            log_out.write(log_line)

            log_out.close()


        debgu_n = False 
        if( debgu_n ):
            print "self"
            print self.ptclC
            print self.bondC

            sys.exit("debug 2")

        return self
