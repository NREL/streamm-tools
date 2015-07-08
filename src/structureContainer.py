"""
Class data structures for atomic data
"""

import copy
import numpy as np 
import json
import sys
import os

from particles     import Particle, ParticleContainer
from bonds         import Bond,     BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer
from impropers     import Improper, ImproperContainer
# from periodictable import periodictable

import pbcs, units

sperator_line = "---------------------------------------------------------------------\n"

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

    def __init__(self, ptclC=ParticleContainer(), bondC=BondContainer(), angleC=AngleContainer(), dihC=DihedralContainer(), impC=ImproperContainer(), verbose=True):
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

        if isinstance(impC, ImproperContainer):
            self.impC = copy.deepcopy(impC)
        else:
            print "5rd arg should be an ImproperContainer object"
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
        del self.boxLengths
        del self.latvec


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
        self.impC   = copy.deepcopy(struc.impC)

        self.boxLengths = copy.deepcopy(struc.boxLengths)
        self.latvec     = copy.deepcopy(struc.latvec)


    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  "\n"
        strucStr += sperator_line
        strucStr += "    Structure properties \n"
        strucStr += sperator_line
        strucStr += "      Box lengths: \n"
        strucStr += "        Lx (A) = " + str(self.boxLengths[0]) + "\n"
        strucStr += "        Ly (A) = " + str(self.boxLengths[1]) + "\n"
        strucStr += "        Lz (A) = " + str(self.boxLengths[2]) + "\n"
        strucStr += "      Volume %f  A^3 \n"%self.getVolume()
        strucStr += "      Mass %f  AMU \n"%self.getTotMass()     
        strucStr += "      Density %f g/cm^3 \n"%self._getDensity()     
        strucStr += "      Lattice vectors \n"
        latvec_i = self.getLatVec()
        strucStr += "        v_i (A)  ( %f , %f , %f ) \n"%(latvec_i[0][0],latvec_i[0][1],latvec_i[0][2])
        strucStr += "        v_j (A)  ( %f , %f , %f ) \n"%(latvec_i[1][0],latvec_i[1][1],latvec_i[1][2])
        strucStr += "        v_k (A)  ( %f , %f , %f ) \n"%(latvec_i[2][0],latvec_i[2][1],latvec_i[2][2])
        strucStr += "\n"                                        
        strucStr += "      Particles %d \n"%(len(self.ptclC))
        strucStr += "      Bonds  %d \n"%(len(self.bondC))
        strucStr += "      Angles %d \n"%(len(self.angleC))
        strucStr += "      Dihedrals %d \n"%(len(self.dihC))
        strucStr += "      Impropers %d \n"%(len(self.impC))

        if( self.verbose ):
            # Print
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
        impC   = copy.deepcopy(other.impC)    #  inside can be changed (for adding below)
        
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

        impC.replacePtclIDs(idFromToDict)   # Use tracked list of ID changes
        self.impC += impC                   # Now add impC with 'corrected' IDs

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

        # ------------ Redo improper dihedrals-----------
        self.impC.replacePtclIDs(idFromToDict)  # Use tracked list of ID changes angles

        localDict = dict()
        for toDihedralID, dihedralTuple in enumerate(self.impC): # Enumerate returns (ID, obj) tuple for ptclTuple
            toDihedralID +=1                                     # Sets reordering index correctly
            fromDihedralID = dihedralTuple[0]                    # Picks out ID from ptclTuple
            dihedralObj = self.impC[fromDihedralID]              # Get dihedral object
            localDict[toDihedralID] = dihedralObj                # Set local dict with reordered index

        del self.impC.impropers                    # Ensure memory is free
        self.impC.impropers = localDict            # Assign reordered local copy to bond container
        del localDict



    def getSubStructure(self, ptclIDList, particlesOnly=False):
        """
        Return a new Structure object with partcleID's in input list
        Preserves IDs of particles (and any bonds, angles, dihedrals...)
        Bonds, angles, dihedrals... etc are included if and ONLY if all the
        particles of which it consists is in the ptclIDList

        Args:
            ptclIDList    (list): global particles ID's for which to return structure
            particlesOnly (bool): Flag for including ParticleContainer only (executes faster) Default includes any angles, bonds, dihedrals
        Return:
            New StructureContainer() object. IDs in new object are unique
        """

        if (len(self)==0 and len(ptclIDList)>0):
            print "Error: getSubStructure using non-zero ptcl list on empty container"
            sys.exit(0)

        # Grab particles from IDlist and put into sub-particle container
        subAtoms = ParticleContainer(ptclIDList) # Initial ptcl container w/input IDs
        for pgid in ptclIDList:
            atom = self.ptclC[pgid]
            subAtoms[pgid] = atom


        if (not particlesOnly):

            # For each bond object in container check that both
            # particles in bond are in ptcl search list
            bondIDList = self.bondC.keys()             # Get keys of bond container
            subBonds   = BondContainer(bondIDList)     # Intitialize subbond container
            for gid, bondObj in self.bondC:
                if ( (bondObj.pgid1 in ptclIDList) and \
                     (bondObj.pgid2 in ptclIDList) ):
                    subBonds[gid] = bondObj
                else:
                    # Need to remove empty key generated above
                    del subBonds[gid]

            # For each angle object in container check that both
            # particles in angle are in ptcl search list
            angleIDList = self.angleC.keys()           # Get keys of angle container
            subAngles   = AngleContainer(angleIDList)  # Initialize subangle container
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
            dihedralIDList = self.dihC.keys()                   # Get keys of dihedral container
            subDihedrals   = DihedralContainer(dihedralIDList)  # Initialize subdihedral container        
            for gid, dObj in self.dihC:
                if ( (dObj.pgid1 in ptclIDList) and \
                     (dObj.pgid2 in ptclIDList) and \
                     (dObj.pgid3 in ptclIDList) and \
                     (dObj.pgid4 in ptclIDList) ):
                    subDihedrals[gid] = dObj
                else:
                    # Need to remove empty key generated above
                    del subDihedrals[gid]


            # For each imp dihedral object in container check that both
            # particles in imp dihedral are in ptcl search list
            impdihedralIDList = self.impC.keys()                 # Get keys of imp dihedral container
            subimpDihedrals = ImproperContainer(impdihedralIDList)  # Initialize subdihedral container        
            for gid, dObj in self.impC:
            	if ( (dObj.pgid1 in ptclIDList) and \
                 (dObj.pgid2 in ptclIDList) and \
                 (dObj.pgid3 in ptclIDList) and \
                 (dObj.pgid4 in ptclIDList) ):
                     subimpDihedrals[gid] = dObj
                else:
                     # Need to remove empty key generated above
                     del subimpDihedrals[gid]


        else:
            if self.verbose:
                print "getSubStructure including only ParticleContainer"


        if (not particlesOnly):
            return StructureContainer(subAtoms, subBonds, subAngles, subDihedrals, subimpDihedrals)
        else:
            return StructureContainer(subAtoms)


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


    def setLatVec(self, latvec):
        """
        Set lattice vector

        Args:
            latvec (list): 3-element lattice vector
        """
        self.latvec = latvec


    def getLatVec(self):
        """
        Get lattice vector 

        Returns:
               (list) 3-element lattice vector list
        """
        return self.latvec


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


    def getPtclNum(self):
        """
        Return number of particles in a structure
        """
        return len(self.ptclC)


    def expandLatVec(self, percent ):
        """
        Expand lattice vector by a certain percent
        NOTE: needs to consider units
        """
        self.latvec[0] = self.latvec[0]*(1.0+percent)
        self.latvec[1] = self.latvec[1]*(1.0+percent)
        self.latvec[2] = self.latvec[2]*(1.0+percent)

    def getLatVec(self):
        """
        Get lattice vector
        NOTE: needs to consider units
        """
        return self.latvec


    def getVolume(self):
        """
        Calculate volume
        NOTE: needs to consider units

        Method:
            Volume = ( v_i x v_j ) . v_k
        """
        br1 = np.cross(self.latvec[0],self.latvec[1])
        vol = np.dot(br1,self.latvec[2])
        return vol


    def getTotMass(self):
        """
        Calculate total mass of system 
        NOTE: needs to consider units
        """
        # Sum mass, charges
        total_mass = 0.0
        for pid, ptclObj in self.ptclC :
            total_mass += ptclObj.mass

        return float(total_mass)


    #
    # Private class methods
    #

    def _getDensity(self):
        """
        Calculate density of system in AMU/A^3 and convert to g/cm^3
        NOTE: mass units contained in PtclConatiner
        """

	volume_i = self.getVolume()    
	total_mass_i = self.getTotMass()
	density_i = units.convert_AMUA3_gcm3(total_mass_i/volume_i) 
	
	return density_i
    #########################################################


    """
    TWK: Testing removing method so pbcs module can be erased from release.
    Travis will either reimplement elsewhere or remove permanently
        
    def getlength(self):

        Calculate length from maximum seperation between particles
        NOTE: pbcs also calls methods from StructureContainer and this should be moved.
        NOTE: pbcs functionality should be moved here since lattice vectors are here

        Return
          struc_len (float) 

        
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
    """


    def vec_shift(self,r_shift):
        """
        Shift structure by vector
        NOTE: only called in pbcs.py

        Args:
          r_shift (numpy vector) to shift all the cordinates by
        """
    
        for pid, ptclObj in self.ptclC :
             r_i = np.array( ptclObj.position ) + r_shift
             ptclObj.position = [r_i[0],r_i[1],r_i[2]]

            
    def center_mass(self):
        """
        Find center of mass of a structure
        NOTE: needs to consider units

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
        NOTE: keep here (maybe re-name method) only in pbcs.py

        Return
          r_shift (numpy array) position of the center of mass

        """
        r_mass = self.center_mass()
        r_m_s = r_shift - r_mass
        self.vec_shift(r_m_s)

        
    def rotate(self, rot_angle_i, rot_angle_j):
        """
        Rotate particles in particle container. Rotation angle i around y axis
        and rotation angle j around z-axis

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
        

    def maxminLatVec(self):
        """
        Set lattice vector based on the max/min position of the particles
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
        NOTE: only in replicate
        """
        n_chains = 0
        for pid, ptclObj in self.ptclC :
            if( ptclObj.tagsDict["chain"] > n_chains): n_chains = ptclObj.tagsDict["chain"]

        return n_chains

