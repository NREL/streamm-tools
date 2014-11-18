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
from periodictable import periodictable

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
        strucStr += "      Density %f g/cm^3 \n"%self.getDen()     
        strucStr += "      Lattice vectors \n"
        latvec_i = self.get_latvec()
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

        impdihedralIDList = self.impC.keys()                 # Get keys of dihedral container
        subimpDihedrals = ImproperContainer(impdihedralIDList)  # Initialize subdihedral container        

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


        # For each dihedral object in container check that both
        # particles in dihedral are in ptcl search list
        for gid, dObj in self.impC:
            if ( (dObj.pgid1 in ptclIDList) and \
                 (dObj.pgid2 in ptclIDList) and \
                 (dObj.pgid3 in ptclIDList) and \
                 (dObj.pgid4 in ptclIDList) ):
                subimpDihedrals[gid] = dObj
            else:
                # Need to remove empty key generated above
                del subimpDihedrals[gid]

        return StructureContainer(subAtoms, subBonds, subAngles, subDihedrals, subimpDihedrals)



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

    def get_latvec(self ):
        """
        Get lattice vector 
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
            Volume = ( v_i x v_j ) . v_k
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
        Calculate density of system in AMU/A^3 and convert to g/cm^3

        """

	volume_i = self.getVolume()    
	total_mass_i = self.getTotMass()

	density_i = units.convert_AMUA3_gcm3(total_mass_i/volume_i) 
	
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
        import topology

        debug = False

        if(debug):
            print " list_k ",list_k
            print " list_i ",list_i
            print " list_j ",list_j
            print " list_l ",list_l

        
        # Create neighbor list form bonds
        cov_nblist, cov_nbindx = topology.build_covnablist(self)
        #NBLIST,NBINDEX = self.bonded_nblist()

        
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

	    N_k_o = cov_nbindx[atom_k]
	    N_k_f = cov_nbindx[atom_k+1]

            if(debug): print  "checking k ",atom_k,self.ptclC[atom_k].type," with ",N_k_f - N_k_o," nbs"
	    
	    for indx_i in range( N_k_o,N_k_f):
                atom_i = cov_nblist[indx_i]
                add_i = False

                if( debug ): print " checking i ",atom_i,self.ptclC[atom_i].type
                
                for  p_i in list_i:
                    if( atom_i == p_i ): #and atom_i > atom_k ):
                        add_i = True
                if( add_i ): #atom_i in list_i ):

                    
		    N_i_o = cov_nbindx[atom_i]
		    N_i_f = cov_nbindx[atom_i+1] 
					
		    for indx_j in range( N_i_o,N_i_f):
			atom_j = cov_nblist[indx_j]
			add_j = False

                        if( debug ): print " checking j ",atom_j,self.ptclC[atom_j].type
                            
                        for  p_j  in list_j:
                            if( atom_j == p_j ): #and atom_j > atom_i ):
                                add_j = True
                        if( add_j ): # atom_j  in list_j  ):

                            N_j_o = cov_nbindx[atom_j]
                            N_j_f = cov_nbindx[atom_j+1] 

                            for indx_l in range( N_j_o,N_j_f):
                                atom_l = cov_nblist[indx_l]
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
                        #if( eln_i != self.ptclC[eln_p_cnt].type ):
                        #    print "  Particle ",eln_p_cnt, self.ptclC[eln_p_cnt].type," != ",eln_i

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
                        print " json file contains %d atoms and fchk file contains %d "%(len(self.ptclC),NA)
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
