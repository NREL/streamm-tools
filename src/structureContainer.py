"""
Class data structures for atomic data
"""

import copy
import json
import sys
import os
import numpy as np
import math 

from particles     import Particle, ParticleContainer
from bonds         import Bond,     BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer
from impropers     import Improper, ImproperContainer
# from periodictable import periodictable

# TWK: import pbcs, units
import units
from periodictable import periodictable

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

    def __init__(self, ptclC=None, bondC=None, angleC=None, dihC=None, impC=None, verbose=True):
        """
        Constructor for a composite structure. Deepcopy of containers is used
        so this is the main container that manages memory for all sim objects

        Args:
            ptclC  (ParticleContainer)  
            bondC  (BondContainer)  
            angleC (AngleContainer)
            dihC   (DihedralContainer)
            verbose (bool) -- Flag for debug/status messages. Default=True
        """
        
        # Creating default containers here so constructor does not
        # instantiate at definition time
        if ptclC is None:
            ptclC = ParticleContainer()
        if bondC is None:
            bondC = BondContainer()
        if angleC is None:
            angleC = AngleContainer()
        if dihC is None:
            dihC = DihedralContainer()
        if impC is None:
            impC = ImproperContainer()            
            

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
        self.latvec = [  [100.0,0.0,0.0] ,  [0.0,100.0,0.0] , [0.0,0.0,100.0] ]


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
        strucStr += "      Volume %f  A^3 \n"%self.getVolume_c()
        strucStr += "      Mass %f  AMU \n"%self.getTotMass()     
        strucStr += "      Density %f g/cm^3 \n"%self.getDensity()     
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


    def getVolume_c(self):
        """
        Calculate volume of Orthorhombic unit cell 
        Method:
        none cubic Volume = ( v_i x v_j ) . v_k / cubic Volume =v_i  v_j  v_k
        """
        vol = self.latvec[0][0]*self.latvec[1][1]*self.latvec[2][2]
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
    def getDensity(self):
        """
        Calculate density of system in AMU/A^3 and convert to g/cm^3
        NOTE: mass units contained in PtclConatiner
        """

	volume_i = self.getVolume_c()    
	total_mass_i = self.getTotMass()
	density_i = units.convert_AMUA3_gcm3(total_mass_i/volume_i) 
	
	return density_i


    def bondC_nblist(self):
        """
        Create neighbor list of bonded particles
        """

        debug = False

        NNAB  = 0

        maxnnab = len(self.bondC)*2 + 1

        # python style nieghbor list
        nblist_py = [] #numpy.empty( maxnnab,  dtype=int )

        self.bonded_nblist = []
        self.bonded_nbindx = []
        self.bonded_nblist.append( 0 )

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
        for nblis_indx in range( len(nblist_py)):
            # loop over  each particle p_i and get list of neighbors nlist_i
            nlist_i = nblist_py[nblis_indx]
            self.bonded_nbindx.append( NNAB + 1 )
            # Loop over neighbor list of each particle nlist_i and get neighbor p_j
            p_i = nblis_indx #+ 1
            for p_j in  nlist_i:

                if(debug):
                    print "nblis_indx",nblis_indx
                    print " p_i ",p_i," p_j ",p_j

                #if( p_j > p_i):
                # remove redundent neighbors 
                NNAB +=  1
                # add to neighbor list 
                self.bonded_nblist.append( p_j )

        self.bonded_nbindx.append( NNAB + 1 )

        if ( debug ):
            print ' total nbs ',NNAB

            for p_i, prtclC in self.ptclC:
                N_i_o = self.bonded_nbindx[p_i]
                N_i_f = self.bonded_nbindx[p_i+1]
                print " atom ",p_i,prtclC.type, " has ",N_i_f - N_i_o

                for indx_j in range( N_i_o,N_i_f):
                    atom_j = self.bonded_nblist[indx_j]
                    print "      nb ",atom_j, self.ptclC[atom_j].type

            sys.exit('bonded_nblist debug')

        return 

    def build_bonded_nblist(self,max_nn=12.0,radii_buffer=1.25):
        """
        Build neighbor list from particle radii

        Arguments:
           max_nn        (int)  maximum nearest neighbors to set the size of the nieghbor list 
           radii_buffer (float) buffer to modify radii 
        """

        maxnnab = float(len( self.ptclC ))*float(max_nn)
        
        self.bonded_nblist = np.empty( maxnnab,  dtype=int )
        self.bonded_nbindx = np.empty( maxnnab,  dtype=int )

        NNAB = 0
        if( len( self.ptclC) <= 0 ):
            error_line = " empyt particle container passed to  build_covnablist "
            sys.exit(error_line)

        for pid_i, ptclObj_i  in self.ptclC:
            self.bonded_nbindx[pid_i] = NNAB + 1
            #r_i = np.array( [ float( ptclObj_i.position[0] ),float( ptclObj_i.position[1] ),float( ptclObj_i.position[2] )] )
            r_i = np.array( ptclObj_i.position  )
            rc_i = ptclObj_i.radii*radii_buffer

            NNAB_i = NNAB 
            for pid_j, ptclObj_j  in self.ptclC:
                if( pid_j != pid_i ):
                    r_j = np.array(  ptclObj_j.position )
                    r_ij = r_j - r_i
                    mag_dr =  np.linalg.norm(r_ij)
                    #r_ij = delta_r(r_i,r_j)
                    rc_j = ptclObj_j.radii*radii_buffer
                    r_cov = rc_i + rc_j

                    if( mag_dr <= r_cov ):
                        NNAB = NNAB + 1
                        self.bonded_nblist[NNAB] =  pid_j
                        
        # Account for final atom position
        self.bonded_nbindx[pid_i+1] =  NNAB + 1

        return 


    def nblist_bonds(self,debug=False):
        """
        Generate bonds from bonded neighbor list
        """

        if( debug):
            print " Initial # of bonds ",len(self.bondC)

        for pid_i, ptclObj_i  in self.ptclC:    
            N_o = self.bonded_nbindx[pid_i ]
            N_f = self.bonded_nbindx[pid_i+ 1 ] 
            for pid_j in self.bonded_nblist[N_o:N_f]:
                if( pid_j > pid_i):
                    b_i = Bond( pid_i, pid_j ) 
                    self.bondC.put(b_i)
                    if( debug ):
                        print pid_i, pid_j

        if( debug):
            print " Final # of bonds ",len(self.bondC)
            sys.exit(" debug nblist_bonds")

        return 


    def nblist_angles(self,debug=False):
        """
        Generate angles from bonded neighbor list 
        """
        if( debug):
            print " Initial # of angles ",len(self.angleC)

        for pid_i, ptclObj_i  in self.ptclC:    
            N_o = self.bonded_nbindx[pid_i ]
            N_f = self.bonded_nbindx[pid_i+ 1 ]
            NNAB = N_f - N_o 
            if( NNAB >= 2 ):
                for indx_j in range( N_o,N_f-1):
                    pid_j = self.bonded_nblist[indx_j]
                    for indx_k in range( indx_j+1,N_f):
                        pid_k = self.bonded_nblist[indx_k]
                        if ( pid_j != pid_k ):
                            a_i = Angle( pid_k,pid_i, pid_j )            
                            self.angleC.put(a_i)
                            if( debug ):
                                print pid_k,pid_i, pid_j
                            
        if( debug):
            print " Final # of angles ",len(self.angleC)
            sys.exit(" debug nblist_bonds")

        return
        
    def nblist_dih(self,debug=False):
        """
        Generate dihedrals from bonded neighbor list 
        """
        if( debug):
            print " Initial # of dih angles ",len(self.dihC)

        for pid_i, ptclObj_i  in self.ptclC:    
            N_o = self.bonded_nbindx[pid_i ]
            N_f = self.bonded_nbindx[pid_i+ 1 ] - 1
            NNAB = N_f - N_o + 1
            for indx_j in range( N_o,N_f+1):
                pid_j = self.bonded_nblist[indx_j]
                if( pid_j > pid_i):
                    dih_ij_cnt = 0
                    atom_k_i = -1
                    atom_l_i = -1
                    for indx_k in range( N_o,N_f+1 ):
                        pid_k = self.bonded_nblist[indx_k]
                        if ( pid_k != pid_j ):
                            No_j =  self.bonded_nbindx[ pid_j  ]
                            Nf_j = self.bonded_nbindx[pid_j+ 1 ] - 1
                            for indx_l in range( No_j,Nf_j+1):
                                pid_l = self.bonded_nblist[indx_l]
                                if ( pid_l != pid_i and pid_l != pid_k ):
                                    d_i = Dihedral( pid_k, pid_i, pid_j, pid_l )            
                                    self.dihC.put(d_i)
                                    if( debug ):
                                        print pid_k, pid_i, pid_j, pid_l


        if( debug):
            print self.dihC
            print " Final # of dih ",len(self.dihC)
            sys.exit(" debug nblist_bonds")
            
        return

    def find_max_qgroup_id(self,debug = False):
        """
        Find maximum qgroup value 
        """        
        self.max_qgroup_id = 0
        for pid, ptclObj  in self.ptclC:
            if( debug ): print "qgroup", ptclObj.tagsDict["qgroup"]
            if(ptclObj.tagsDict["qgroup"] > self.max_qgroup_id ): self. max_qgroup_id = ptclObj.tagsDict["qgroup"]
                
        return
    
    def set_label(self,debug = False):
        """
        Set particle label 
        """
        pt = periodictable()
        self.composition = np.zeros(pt.maxgid,dtype=np.int)
        
        for pid, ptclObj  in self.ptclC:
            atomic_number  = ptclObj.tagsDict["number"]
            self.composition[atomic_number] += 1 
            ptclObj.tagsDict["label"] = "{}{}".format(ptclObj.tagsDict["symbol"],self.composition[atomic_number])

        return


    def calc_nnab(self,pid_i,nbindx):
        """
        Return number of nieghbors for a given particle  
        """
        return nbindx[pid_i+1] - nbindx[pid_i]


    def calc_elcnt(self,pid_i,nblist,nbindx):
        """
        Return 
        """
        #
        # Find number of elements 
        #
        pt = periodictable()
        ELCNT = np.zeros(pt.maxgid,dtype=np.int)
        N_o = nbindx[pid_i]
        N_f = nbindx[  pid_i +1  ] - 1 
        for indx_j in range( N_o,N_f+1):
            j = nblist[indx_j]
            el_j = int( self.ptclC[j].tagsDict["number"] )
            if( el_j >= 0 ):
                ELCNT[el_j] = ELCNT[el_j] + 1

        return ELCNT 


    def find_ring_type(self,pid_i,debug = False):
        """
        Criteria for a conjugated particle 
        """

        number_j = self.ptclC[pid_i].tagsDict["number"]
        NNAB_j = self.calc_nnab(pid_i,self.bonded_nbindx)
        ELCNT_j = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)

        if(debug):
            print '            atomic # %d  NN  %d '%(number_j, NNAB_j)

        # Test if atom j is conjugated 
        if ( number_j == 6 and NNAB_j == 3 ):                       return True
        if ( number_j == 16 and NNAB_j == 2 and ELCNT_j[6] == 2 ):  return True
        if ( number_j == 7 and NNAB_j >= 2  ):                      return True

        return False


    def id_ring3(self,pid_o,debug = False):
        """
        Find atoms in conjugated rings 
        """
        ptclObj_o = self.ptclC[pid_o]
        
        max_paths = 1000

        R_SET = []
        RING_ATOMS = []
        BAD_PATH = []

        #a_i = pid_o -1 
        one = 1
        zero = 0
        atoms_in_ring = 0
        # relabel based on neighbors
        n_ptcl = len(self.ptclC)

        if(debug):
            print " for a system of %d partiles checking particle %d "%(n_ptcl,pid_o)

        for pid_i in range(n_ptcl+1): 
            R_SET.append(one)

        R_SET[pid_o] = 0 
        r_term = 1
        cnt = 0
        p_cnt = 0
        if(debug): print ' initializing ring ',pid_o," w NN ", self.calc_nnab(pid_o,self.bonded_nbindx)

        last_i = pid_o
        NNAB_last = self.calc_nnab(pid_i,self.bonded_nbindx)
        while ( r_term ):
            N_o = self.bonded_nbindx[last_i]
            N_f = self.bonded_nbindx[last_i+1] - 1

            if( debug): print " checking neighbors ",N_o,N_f

            for n_indx in range( N_o,N_f+1):
                j = self.bonded_nblist[n_indx]
                cnt = cnt + 1
                # If non of the neighboring atoms of atom j are conjugated
                #    make the path as bad 
                if( cnt > NNAB_last+1 ):
                    p_cnt += 1
                    if(debug): print ' bad path found resetting at atom ',j,self.ptclC[j].tagsDict["number"]
                    for ring_a in range(len(RING_ATOMS)):
                        j = RING_ATOMS[ring_a]

                    BAD_PATH.append(j)
                    # Reset lists and counts
                    RING_ATOMS = []
                    for pid_k in range(n_ptcl):
                        R_SET[pid_k] = 1
                    atoms_in_ring = 0
                    R_SET[pid_o] = 0 
                    cnt = 0
                    last_i = pid_o

                    if(debug): print '  resetting last_i to ',pid_o,ptclObj_o.tagsDict["number"]

                    for bad_i in range( len(BAD_PATH)):
                        j = BAD_PATH[bad_i]
                        R_SET[j] = 0
                        if(debug): print '  bad path atoms ',self.ptclC[j].tagsDict["number"]
                    break

                # If path traces back to original atom ring has been found 
                if( atoms_in_ring  > 1 and j == pid_o ):
                    r_term = 0
                    if(debug): print ' ring found with ',atoms_in_ring 
                    break

                number_j = self.ptclC[j].tagsDict["number"]
                NNAB_j = self.calc_nnab(j,self.bonded_nbindx)
                ELCNT_j = self.calc_elcnt(j,self.bonded_nblist,self.bonded_nbindx)
                ring_type = False 

                if(debug):
                    print '           with ', number_j
                    print '           with ', NNAB_j,
                    print '           with ', R_SET[j] 

                # Test if atom j is conjugated 
                if ( number_j == 6 and NNAB_j == 3 and R_SET[j] == 1 ):                      ring_type = True
                if ( number_j == 16 and NNAB_j == 2 and R_SET[j] == 1 and ELCNT_j[6] == 2 ): ring_type = True
                if ( number_j == 7 and NNAB_j >= 2 and R_SET[j] == 1 ):                      ring_type = True
                if( ring_type ):
                    atoms_in_ring = atoms_in_ring + 1
                    RING_ATOMS.append(j)
                    R_SET[j] = 0
                    last_i = j
                    cnt = 0
                    NNAB_last = NNAB_j 
                    if(debug): print '   atom ',j,number_j,' added '
                    break

            if (p_cnt > max_paths ):
                if(debug): print '     max paths considered ',max_paths
                r_term = 0

        return RING_ATOMS

                     
    def find_rings(self,debug = False):
        """
        Find conjugate rings
        """

        self.ring_nblist = []
        self.ring_nbindex = []
        RING_NUMB = []

        IN_RING = []

        n_ptcl = len(self.ptclC)
        one = 0 
        zero = 0 
        ring_cnt = 0
        a_cnt = 0
        # relabel based on neighbors
        for pid_i, ptclObj_i  in self.ptclC:
            self.ring_nblist.append(zero)
            self.ring_nbindex.append(zero)
            IN_RING.append(one)
            RING_NUMB.append(zero)

        RING_NUMB.append(zero)
        self.ring_nblist.append(zero)
        self.ring_nbindex.append(zero)


        # relabel based on neighbors
        for pid_i, ptclObj_i  in self.ptclC:
            RING_ATOMS = []
            # If particle i is conjugated 
            if( self.find_ring_type(pid_i) ):
                RING_ATOMS = self.id_ring3(pid_i)
            if( debug ):
                print pid_i,RING_ATOMS

            if ( len(RING_ATOMS) > 1 ):
                # If member of ring alread part of another ring add new ring to existing ring
                r_numb = 0 
                for ring_a in range(len(RING_ATOMS)):
                    j = RING_ATOMS[ring_a]
                    if( RING_NUMB[j] != 0 ):
                        r_numb = RING_NUMB[j]
                if( r_numb == 0 ):
                    ring_cnt = ring_cnt + 1
                    r_numb = ring_cnt 
                    #print ' new ring ',ring_cnt

                for ring_a in range(len(RING_ATOMS)):
                    j = RING_ATOMS[ring_a]
                    if( RING_NUMB[j] == 0 ): 
                        a_cnt = a_cnt + 1
                        RING_NUMB[j] = r_numb

        a_cnt = 0
        for r_numb in range(1,ring_cnt+1):
            self.ring_nbindex[r_numb] = a_cnt + 1
            for i in range(n_ptcl):
                if( RING_NUMB[i] == r_numb ):
                    a_cnt = a_cnt + 1
                    self.ring_nblist[a_cnt] = i

        self.ring_nbindex[ring_cnt+1] = a_cnt + 1


        # Set attached H to have same ring #
        attach_H = False
        if( attach_H ):

            for r_numb in range(1,ring_cnt+1):
                Nr_o = self.ring_nbindex[r_numb]
                Nr_f = self.ring_nbindex[r_numb+1] - 1
                print ' Ring ', r_numb
                for r_indx in range(Nr_o,Nr_f+1):
                    i = self.ring_nblist[r_indx]
                    N_o = self.bonded_nbindx[i]
                    N_f = self.bonded_nbindx[i+1] - 1
                    for n_indx in range( N_o,N_f+1):
                        j = self.bonded_nblist[n_indx]
                        if( self.ptclC[j].tagsDict["number"] == 1 ):
                            RING_NUMB[j] = r_numb


        # Update particles
        for pid_i, ptclObj_i  in self.ptclC:
            ptclObj_i.tagsDict["ring"] =  RING_NUMB[pid_i]
            if( debug ):
                print ' atom  ',pid_i,ptclObj_i.tagsDict["number"],ptclObj_i.tagsDict["ring"]

        if(debug):
            print ring_cnt , " rings found "
            for r_numb in range(1,ring_cnt+1):
                N_o = self.ring_nbindex[r_numb]
                N_f = self.ring_nbindex[r_numb+1] - 1
                print ' Ring ', r_numb
                for r_indx in range(N_o,N_f+1):
                    i = self.ring_nblist[r_indx]
                    print ' type ',self.ptclC[i].tagsDict["number"],' or ',i,RING_NUMB[i]
                    #sys.exit('find_rings')

            for pid_i, ptclObj_i  in self.ptclC:    
                print ' atom  ',pid_i,ptclObj_i.tagsDict["number"],ptclObj_i.tagsDict["ring"]

            sys.exit('top.find_rings')

        return  #(  self.ring_nblist, self.ring_nbindex  )


    def atomtypes(self, update_chr = False ):
        """
        Set OPLSaa atom types 
        """
    

        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            if ptclObj_i.tagsDict["number"] == 6 :

                # simple guess based on coordination 
                if int(NNAB) == 4 :
                    ptclObj_i.tagsDict["fftype"] = 'CT' # Alkane
                    # refine guess based on nieghbors 
                    if( ELCNT[1] == 4 ):                              # Methane 
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.24
                    elif( ELCNT[6] == 1 and ELCNT[1] == 3):         # Methyl 
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.18
                    elif(   ELCNT[6] == 2 and ELCNT[1] == 2 ):      
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.12
                    elif( ELCNT[6] == 3 and ELCNT[1] == 1 ):         # Alkane 
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.06
                    elif(  ELCNT[6] == 4 ):
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  0.0
                    elif(  ELCNT[1] == 1 and ELCNT[8] == 2 ):       #acetal
                        ptclObj_i.tagsDict["fftype"] = 'CO'
                        if( update_chr ): ptclObj_i.charge =  -0.4

                if  int(NNAB) == 3 :
                    r_numb = ptclObj_i.tagsDict["ring"]
                    if( r_numb != 0 ):
                        nring = self.calc_nnab(r_numb,self.ring_nbindex)
                        nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                        ptclObj_i.tagsDict["fftype"] = 'CA' # Aromatic C
                        if( ELCNT[6] == 2 and ELCNT[1] == 1 ):           # "Aromatic C"  
                            ptclObj_i.tagsDict["fftype"] = 'CA'
                            if( update_chr ): ptclObj_i.charge =  -0.1150
                        if( ELCNT[6] == 3  and nring_mod == 0 ):            # "Naphthalene Fusion C"
                            ptclObj_i.tagsDict["fftype"] = 'CA'
                            if( update_chr ): ptclObj_i.charge =  0.0

                        debug = 0
                        if( debug  ):
                            print " r_numb ",nring,r_numb, nring_mod, ELCNT[6], ELCNT[7], ptclObj_i.tagsDict["fftype"]
                            # sys.exit(' debug ')

                    else: # not aromatic
                        if( ELCNT[6] == 2 and ELCNT[1] == 1 ):          # diene 
                            ptclObj_i.tagsDict["fftype"] = 'C='
                            if( update_chr ): ptclObj_i.charge =  0.0
                        elif( ELCNT[6] == 2 and ELCNT[8] == 1 ):          # Benzophenone
                            ptclObj_i.tagsDict["fftype"] = 'C'
                            if( update_chr ): ptclObj_i.charge =  0.7
                        elif( ELCNT[6] == 1 and ELCNT[8] == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'C' # pmma
                            if( update_chr ): ptclObj_i.charge =  0.7

                if int(NNAB) == 2 :
                    ptclObj_i.tagsDict["fftype"] = 'C:'    # Allene
                    if( ELCNT[6] == 1 and ELCNT[7] == 1 ):   # "Benzonitrile -CN"  
                        ptclObj_i.tagsDict["fftype"] = 'CZ'

                if int(NNAB) == 1 :
                    ptclObj_i.tagsDict["fftype"] = '' # Aromatic C
                    error_line =  " WARNING!!! carbon index ",pid_i," bonded to single atom "
                    sys.exit(error_line)
            #
            # label oxygens
            #
            if( ptclObj_i.tagsDict["number"] == 8 ):
                if int(NNAB) == 1 :
                    ptclObj_i.tagsDict["fftype"] = 'O' # double bonded
                    if( update_chr ): ptclObj_i.charge =  -0.5
                if int(NNAB) == 2 :
                    ptclObj_i.tagsDict["fftype"] = 'OS' # ether
                    if( update_chr ): ptclObj_i.charge =  -0.5
                    if( ELCNT[1] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'OH' # Alcohol
                        if( update_chr ): ptclObj_i.charge =   -0.6830
                    if( ELCNT[8] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'O2' # Carboxylate
                        if( update_chr ): ptclObj_i.charge =   -0.800
                    if( ELCNT[16] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'OY' # Sulfoxide
                        if( update_chr ): ptclObj_i.charge =   -0.4200
                if int(NNAB) == 3 :
                    if( ELCNT[7] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'ON'
                        if( update_chr ): ptclObj_i.charge =  -0.118

            #
            # label nitrogens 
            #
            if ptclObj_i.tagsDict["number"] == 7 :
                if int(NNAB) == 3 :      # amide
                    ptclObj_i.tagsDict["fftype"] = 'N' 
                    if( ELCNT[1] == 3 ): # Ammonia NH3"   
                        ptclObj_i.tagsDict["fftype"] = 'NT'
                    if( ELCNT[1] == 2 ): # -NH2
                        ptclObj_i.tagsDict["fftype"] = 'N2'
                    if( ELCNT[1] == 1 ): # -NH2
                        ptclObj_i.tagsDict["fftype"] = 'N2'

                r_numb = ptclObj_i.tagsDict["ring"]
                if( r_numb != 0 ):
                    nring = self.calc_nnab(r_numb,self.ring_nbindex)
                    nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                    if ( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                        if ( int(NNAB) == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'NC'      # "Imidazole N3"     
                        if  int(NNAB) == 3 :
                            ptclObj_i.tagsDict["fftype"] = 'NA'      # "Imidazole N1"      
                    else:
                        if ( int(NNAB) == 2  ):
                            ptclObj_i.tagsDict["fftype"] = 'NB'
                if ( int(NNAB) == 1 and ELCNT[6] == 1 ) :      # Nitrile 
                    ptclObj_i.tagsDict["fftype"] = 'NZ'


            #
            # label sulfurs
            #
            if( ptclObj_i.tagsDict["number"] == 16 ):
                if int(NNAB) == 2 :
                    ptclObj_i.tagsDict["fftype"] = 'S'   #  Thioether RSR (UA)
                    if( update_chr ): ptclObj_i.charge =-0.4350
                    if( ELCNT[1] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'SH'
                        if( update_chr ): ptclObj_i.charge =  -0.335
                    if( ELCNT[1] == 2  ):
                        ptclObj_i.tagsDict["fftype"] = 'SH'
                        if( update_chr ): ptclObj_i.charge =  -0.470

                if ( int(NNAB) == 4 and ELCNT[8] == 2 and ELCNT[7] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'SY'   #  "Sulfonamide -SO2N<" 
                if ( int(NNAB) == 3 and ELCNT[8] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'SZ'   #  Sulfoxide   

            # Label chlorine 
            if ( ptclObj_i.tagsDict["number"] == 17):
                ptclObj_i.tagsDict["fftype"] = 'Cl'
                if( update_chr ): ptclObj_i.charge =  -0.2

            # Label Fluorine   
            if ( ptclObj_i.tagsDict["number"] == 9):
                ptclObj_i.tagsDict["fftype"] = 'F'
                if( update_chr ): ptclObj_i.charge =  -0.2057

        #
        # label hydrogens
        #

        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            if( ptclObj_i.tagsDict["number"] == 1 ):
                if ( NNAB > 1  ):
                    sys.exit(' over coordinated H')
                pid_j = self.bonded_nblist[ self.bonded_nbindx[pid_i] ]
                ptclObj_j = self.ptclC[pid_j]
                ELCNT_j = self.calc_elcnt(pid_j,self.bonded_nblist,self.bonded_nbindx)
                if ( ptclObj_j.tagsDict["fftype"]== 'CA' ):
                    ptclObj_i.tagsDict["fftype"] = 'HA' #
                    if( update_chr ): ptclObj_i.charge =  0.115
                if ( ptclObj_j.tagsDict["fftype"]== 'CT' ):
                    ptclObj_i.tagsDict["fftype"] = 'HC' #
                    if( update_chr ): ptclObj_i.charge =  0.06
                if ( ptclObj_j.tagsDict["fftype"]== 'CW' ):
                    ptclObj_i.tagsDict["fftype"] = 'HA' #
                if ( ptclObj_j.tagsDict["fftype"]== 'CS' ):
                    ptclObj_i.tagsDict["fftype"] = 'HA' #
                if( ELCNT_j[6] == 2 and ELCNT_j[8] == 1  ): # "Ester -OCH<"
                    ptclObj_i.tagsDict["fftype"] = 'H1'
                    if( update_chr ): ptclObj_j.charge =  0.03
                if ( ptclObj_j.tagsDict["fftype"]== 'SH' ):
                    ptclObj_i.tagsDict["fftype"] = 'HS' #
                    if( update_chr ): ptclObj_i.charge = 0.1550
                if ( ptclObj_j.tagsDict["number"] == 7 ):
                    ptclObj_i.tagsDict["fftype"] = 'H' #
                if ( ptclObj_j.tagsDict["fftype"]== 'NA' ):
                    ptclObj_i.tagsDict["fftype"] = 'H' #
                if ( ptclObj_j.tagsDict["fftype"]== 'N2' ):
                    ptclObj_i.tagsDict["fftype"] = 'H' #
                if ( ptclObj_j.tagsDict["fftype"]== 'OH' ):
                    ptclObj_i.tagsDict["fftype"] = 'HO' #


        # relabel based on neighbors
        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            N_o = self.bonded_nbindx[ pid_i ]
            N_f = self.bonded_nbindx[ pid_i + 1 ] - 1
            nring_i = ptclObj_i.tagsDict["ring"]
            if ( ptclObj_i.tagsDict["fftype"] == 'CA' ):
                for indx_j in range( N_o,N_f+1):
                    pid_j = self.bonded_nblist[indx_j]
                    ptclObj_j = self.ptclC[pid_j]
                    ELCNT_j =  self.calc_elcnt(pid_j,self.bonded_nblist,self.bonded_nbindx)
                    nring_j = ptclObj_j.tagsDict["ring"]
                    if ( ptclObj_j.tagsDict["fftype"]== 'CW' ):
                        if( ELCNT[1] == 1 and ELCNT[6] == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'CS'
                    if(  ELCNT[6] == 2 and nring_i == nring_j and ptclObj_j.tagsDict["number"] != 6 ):          # fussed 
                        if( debug):
                            print ' CB ',pid_i+1, ELCNT[6] ,  nring_i , nring_j 
                        ptclObj_i.tagsDict["fftype"] = 'CB'
                if ( ptclObj_i.tagsDict["fftype"] == 'CS' ):
                    for indx_j in range( N_o,N_f+1):
                        pid_j = self.bonded_nblist[indx_j]
                        ptclObj_j = self.ptclC[pid_j]
                        ELCNT_j = self.calc_elcnt(pid_j,self.bonded_nblist,self.bonded_nbindx)
                        if (  ptclObj_j.tagsDict["fftype"]== 'CA' ):
                            if ( ELCNT_j[6] == 3 ):
                                ptclObj_j.tagsDict["fftype"]= 'CS'
                            if ( ELCNT_j[6] == 2 and ELCNT_j[1] == 1 ):
                                ptclObj_j.tagsDict["fftype"]= 'CS'

            if ptclObj_i.tagsDict["fftype"] == 'CT' :
                for indx_j in range( N_o,N_f+1):
                    pid_j = self.bonded_nblist[indx_j]
                    ptclObj_j = self.ptclC[pid_j]
                    if ( ptclObj_j.tagsDict["fftype"]== 'O' ):
                        if( ELCNT[1] == 3 and ELCNT[8] == 1 ):
                            ptclObj_i.tagsDict["fftype"] = 'C3' # Methyloxide


        if(debug):   sys.exit('linkers')

        debug = 0
        if(debug):
            for pid_i, ptclObj_i  in self.ptclC:
                print pid_i,ptclObj_i.tagsDict["fftype"],ptclObj_i.tagsDict["number"]
            sys.exit('debug')

        # Check for unidentified atoms
        for pid_i, ptclObj_i  in self.ptclC:
            if ( ptclObj_i.tagsDict["fftype"] == '?' ):
                print ' atom ', pid_i , ptclObj_i.tagsDict["number"],' unknow '
                sys.exit(' Unknow atom ')

        return #(self)
    
    def oplsaa_atomtypes(self, update_chr = False ):
        """
        Set OPLSaa atom types 
        """
    

        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            if ptclObj_i.tagsDict["number"] == 6 :

                # simple guess based on coordination 
                if int(NNAB) == 4 :
                    ptclObj_i.tagsDict["fftype"] = 'CT' # Alkane
                    # refine guess based on nieghbors 
                    if( ELCNT[1] == 4 ):                              # Methane 
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.24
                    elif( ELCNT[6] == 1 and ELCNT[1] == 3):         # Methyl 
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.18
                    elif(   ELCNT[6] == 2 and ELCNT[1] == 2 ):      
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.12
                    elif( ELCNT[6] == 3 and ELCNT[1] == 1 ):         # Alkane 
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  -0.06
                    elif(  ELCNT[6] == 4 ):
                        ptclObj_i.tagsDict["fftype"] = 'CT'
                        if( update_chr ): ptclObj_i.charge =  0.0
                    elif(  ELCNT[1] == 1 and ELCNT[8] == 2 ):       #acetal
                        ptclObj_i.tagsDict["fftype"] = 'CO'
                        if( update_chr ): ptclObj_i.charge =  -0.4

                if  int(NNAB) == 3 :
                    r_numb = ptclObj_i.tagsDict["ring"]
                    if( r_numb != 0 ):
                        nring = self.calc_nnab(r_numb,self.ring_nbindex)
                        nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                        ptclObj_i.tagsDict["fftype"] = 'CA' # Aromatic C
                        if( ELCNT[6] == 2 and ELCNT[1] == 1 ):           # "Aromatic C"  
                            ptclObj_i.tagsDict["fftype"] = 'CA'
                            if( update_chr ): ptclObj_i.charge =  -0.1150
                        if( ELCNT[6] == 3  and nring_mod == 0 ):            # "Naphthalene Fusion C"
                            ptclObj_i.tagsDict["fftype"] = 'CA'
                            if( update_chr ): ptclObj_i.charge =  0.0
                        if( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                            if( ELCNT[6] == 3 ):            # Fused 5-6 rings 
                                ptclObj_i.tagsDict["fftype"] = 'CB'
                                if( update_chr ): ptclObj_i.charge =  0.0
                            if( ELCNT[7] == 1 or ELCNT[16] == 1  ):      # Pyrrole or thiophene 
                                if( ELCNT[6] == 1 and ELCNT[1] == 1 ):             #imidazole C2
                                    ptclObj_i.tagsDict["fftype"] = 'CW'   
                                    if( update_chr ): ptclObj_i.charge =  0.22
                                if( ELCNT[6] == 1 and ELCNT[8] == 1 ):             # Not sure about this one
                                    ptclObj_i.tagsDict["fftype"] = 'CW'   
                                    if( update_chr ): ptclObj_i.charge =  0.22
                            if(  ELCNT[7] == 2 ):             #imidazole C2
                                ptclObj_i.tagsDict["fftype"] = 'CR'   
                                if( update_chr ): ptclObj_i.charge =  0.22
                            if(  ELCNT[7] == 1 and ELCNT[16] == 1 ):             #imidazole C2
                                ptclObj_i.tagsDict["fftype"] = 'CR'
                                if( update_chr ): ptclObj_i.charge =  0.22

                        if( ELCNT[6] == 2 and  ELCNT[8] == 1  ):            # CTD strangeness ... 
                            ptclObj_i.tagsDict["fftype"] = 'CW'
                            if( update_chr ): ptclObj_i.charge =  0.0

                        debug = 0
                        if( debug  ):
                            print " r_numb ",nring,r_numb, nring_mod, ELCNT[6], ELCNT[7], ptclObj_i.tagsDict["fftype"]
                            # sys.exit(' debug ')

                    else: # not aromatic
                        if( ELCNT[6] == 2 and ELCNT[1] == 1 ):          # diene 
                            ptclObj_i.tagsDict["fftype"] = 'C='
                            if( update_chr ): ptclObj_i.charge =  0.0
                        elif( ELCNT[6] == 2 and ELCNT[8] == 1 ):          # Benzophenone
                            ptclObj_i.tagsDict["fftype"] = 'C'
                            if( update_chr ): ptclObj_i.charge =  0.7
                        elif( ELCNT[6] == 1 and ELCNT[8] == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'C' # pmma
                            if( update_chr ): ptclObj_i.charge =  0.7

                if int(NNAB) == 2 :
                    ptclObj_i.tagsDict["fftype"] = 'C:'    # Allene
                    if( ELCNT[6] == 1 and ELCNT[7] == 1 ):   # "Benzonitrile -CN"  
                        ptclObj_i.tagsDict["fftype"] = 'CZ'

                if int(NNAB) == 1 :
                    ptclObj_i.tagsDict["fftype"] = '' # Aromatic C
                    error_line =  " WARNING!!! carbon index ",pid_i," bonded to single atom "
                    sys.exit(error_line)
            #
            # label oxygens
            #
            if( ptclObj_i.tagsDict["number"] == 8 ):
                if int(NNAB) == 1 :
                    ptclObj_i.tagsDict["fftype"] = 'O' # double bonded
                    if( update_chr ): ptclObj_i.charge =  -0.5
                if int(NNAB) == 2 :
                    ptclObj_i.tagsDict["fftype"] = 'OS' # ether
                    if( update_chr ): ptclObj_i.charge =  -0.5
                    if( ELCNT[1] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'OH' # Alcohol
                        if( update_chr ): ptclObj_i.charge =   -0.6830
                    if( ELCNT[8] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'O2' # Carboxylate
                        if( update_chr ): ptclObj_i.charge =   -0.800
                    if( ELCNT[16] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'OY' # Sulfoxide
                        if( update_chr ): ptclObj_i.charge =   -0.4200
                if int(NNAB) == 3 :
                    if( ELCNT[7] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'ON'
                        if( update_chr ): ptclObj_i.charge =  -0.118

            #
            # label nitrogens 
            #
            if ptclObj_i.tagsDict["number"] == 7 :
                if int(NNAB) == 3 :      # amide
                    ptclObj_i.tagsDict["fftype"] = 'N' 
                    if( ELCNT[1] == 3 ): # Ammonia NH3"   
                        ptclObj_i.tagsDict["fftype"] = 'NT'
                    if( ELCNT[1] == 2 ): # -NH2
                        ptclObj_i.tagsDict["fftype"] = 'N2'
                    if( ELCNT[1] == 1 ): # -NH2
                        ptclObj_i.tagsDict["fftype"] = 'N2'

                r_numb = ptclObj_i.tagsDict["ring"]
                if( r_numb != 0 ):
                    nring = self.calc_nnab(r_numb,self.ring_nbindex)
                    nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                    if ( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                        if ( int(NNAB) == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'NC'      # "Imidazole N3"     
                        if  int(NNAB) == 3 :
                            ptclObj_i.tagsDict["fftype"] = 'NA'      # "Imidazole N1"      
                    else:
                        if ( int(NNAB) == 2  ):
                            ptclObj_i.tagsDict["fftype"] = 'NB'
                if ( int(NNAB) == 1 and ELCNT[6] == 1 ) :      # Nitrile 
                    ptclObj_i.tagsDict["fftype"] = 'NZ'


            #
            # label sulfurs
            #
            if( ptclObj_i.tagsDict["number"] == 16 ):
                if int(NNAB) == 2 :
                    ptclObj_i.tagsDict["fftype"] = 'S'   #  Thioether RSR (UA)
                    if( update_chr ): ptclObj_i.charge =-0.4350
                    if( ELCNT[1] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'SH'
                        if( update_chr ): ptclObj_i.charge =  -0.335
                    if( ELCNT[1] == 2  ):
                        ptclObj_i.tagsDict["fftype"] = 'SH'
                        if( update_chr ): ptclObj_i.charge =  -0.470

                if ( int(NNAB) == 4 and ELCNT[8] == 2 and ELCNT[7] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'SY'   #  "Sulfonamide -SO2N<" 
                if ( int(NNAB) == 3 and ELCNT[8] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'SZ'   #  Sulfoxide   

            # Label chlorine 
            if ( ptclObj_i.tagsDict["number"] == 17):
                ptclObj_i.tagsDict["fftype"] = 'Cl'
                if( update_chr ): ptclObj_i.charge =  -0.2

            # Label Fluorine   
            if ( ptclObj_i.tagsDict["number"] == 9):
                ptclObj_i.tagsDict["fftype"] = 'F'
                if( update_chr ): ptclObj_i.charge =  -0.2057

        #
        # label hydrogens
        #

        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            if( ptclObj_i.tagsDict["number"] == 1 ):
                if ( NNAB > 1  ):
                    sys.exit(' over coordinated H')
                pid_j = self.bonded_nblist[ self.bonded_nbindx[pid_i] ]
                ptclObj_j = self.ptclC[pid_j]
                ELCNT_j = self.calc_elcnt(pid_j,self.bonded_nblist,self.bonded_nbindx)
                if ( ptclObj_j.tagsDict["fftype"]== 'CA' ):
                    ptclObj_i.tagsDict["fftype"] = 'HA' #
                    if( update_chr ): ptclObj_i.charge =  0.115
                if ( ptclObj_j.tagsDict["fftype"]== 'CT' ):
                    ptclObj_i.tagsDict["fftype"] = 'HC' #
                    if( update_chr ): ptclObj_i.charge =  0.06
                if ( ptclObj_j.tagsDict["fftype"]== 'CW' ):
                    ptclObj_i.tagsDict["fftype"] = 'HA' #
                if ( ptclObj_j.tagsDict["fftype"]== 'CS' ):
                    ptclObj_i.tagsDict["fftype"] = 'HA' #
                if( ELCNT_j[6] == 2 and ELCNT_j[8] == 1  ): # "Ester -OCH<"
                    ptclObj_i.tagsDict["fftype"] = 'H1'
                    if( update_chr ): ptclObj_j.charge =  0.03
                if ( ptclObj_j.tagsDict["fftype"]== 'SH' ):
                    ptclObj_i.tagsDict["fftype"] = 'HS' #
                    if( update_chr ): ptclObj_i.charge = 0.1550
                if ( ptclObj_j.tagsDict["number"] == 7 ):
                    ptclObj_i.tagsDict["fftype"] = 'H' #
                if ( ptclObj_j.tagsDict["fftype"]== 'NA' ):
                    ptclObj_i.tagsDict["fftype"] = 'H' #
                if ( ptclObj_j.tagsDict["fftype"]== 'N2' ):
                    ptclObj_i.tagsDict["fftype"] = 'H' #
                if ( ptclObj_j.tagsDict["fftype"]== 'OH' ):
                    ptclObj_i.tagsDict["fftype"] = 'HO' #


        # relabel based on neighbors
        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            N_o = self.bonded_nbindx[ pid_i ]
            N_f = self.bonded_nbindx[ pid_i + 1 ] - 1
            nring_i = ptclObj_i.tagsDict["ring"]
            if ( ptclObj_i.tagsDict["fftype"] == 'CA' ):
                for indx_j in range( N_o,N_f+1):
                    pid_j = self.bonded_nblist[indx_j]
                    ptclObj_j = self.ptclC[pid_j]
                    ELCNT_j =  self.calc_elcnt(pid_j,self.bonded_nblist,self.bonded_nbindx)
                    nring_j = ptclObj_j.tagsDict["ring"]
                    if ( ptclObj_j.tagsDict["fftype"]== 'CW' ):
                        if( ELCNT[1] == 1 and ELCNT[6] == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'CS'
                    if(  ELCNT[6] == 2 and nring_i == nring_j and ptclObj_j.tagsDict["number"] != 6 ):          # fussed 
                        if( debug):
                            print ' CB ',pid_i+1, ELCNT[6] ,  nring_i , nring_j 
                        ptclObj_i.tagsDict["fftype"] = 'CB'
                if ( ptclObj_i.tagsDict["fftype"] == 'CS' ):
                    for indx_j in range( N_o,N_f+1):
                        pid_j = self.bonded_nblist[indx_j]
                        ptclObj_j = self.ptclC[pid_j]
                        ELCNT_j = self.calc_elcnt(pid_j,self.bonded_nblist,self.bonded_nbindx)
                        if (  ptclObj_j.tagsDict["fftype"]== 'CA' ):
                            if ( ELCNT_j[6] == 3 ):
                                ptclObj_j.tagsDict["fftype"]= 'CS'
                            if ( ELCNT_j[6] == 2 and ELCNT_j[1] == 1 ):
                                ptclObj_j.tagsDict["fftype"]= 'CS'

            if ptclObj_i.tagsDict["fftype"] == 'CT' :
                for indx_j in range( N_o,N_f+1):
                    pid_j = self.bonded_nblist[indx_j]
                    ptclObj_j = self.ptclC[pid_j]
                    if ( ptclObj_j.tagsDict["fftype"]== 'O' ):
                        if( ELCNT[1] == 3 and ELCNT[8] == 1 ):
                            ptclObj_i.tagsDict["fftype"] = 'C3' # Methyloxide


        # Find Ring linkers
        debug = 0
        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            N_o = self.bonded_nbindx[ pid_i ]
            N_f = self.bonded_nbindx[ pid_i + 1 ] - 1
            if ( ptclObj_i.tagsDict["fftype"] == 'CA' or  ptclObj_i.tagsDict["fftype"] == 'CB'  or  ptclObj_i.tagsDict["fftype"] == 'CW' ) :
                if( NNAB == 3 ):  # if sp2
                    for indx_j in range( N_o,N_f+1):
                        pid_j = self.bonded_nblist[indx_j]
                        ptclObj_j = self.ptclC[pid_j]
                        if ( ptclObj_j.tagsDict["fftype"]== 'CA' or  ptclObj_j.tagsDict["fftype"]== 'CB'  or  ptclObj_j.tagsDict["fftype"]== 'CW' )  :
                            NNAB_j = self.calc_nnab(pid_j,self.bonded_nbindx)
                            if( NNAB_j == 3 ):  # if sp2
                                if(  ptclObj_i.tagsDict["ring"] !=  ptclObj_j.tagsDict["ring"] ): # ring linker
                                    if(debug):
                                        print pid_i,' and ',pid_j,' link ',ptclObj_i.tagsDict["fftype"],ptclObj_i.tagsDict["fftype"]
                                        print ' ', ptclObj_i.tagsDict["ring"] , ptclObj_j.tagsDict["ring"]
                                    ptclObj_i.tagsDict["fftype"] = 'C!'
                                    ptclObj_j.tagsDict["fftype"]= 'C!'                                

        # relabel based on neighbors
        debug = 0
        for pid_i, ptclObj_i  in self.ptclC:
            if ( ptclObj_i.tagsDict["fftype"] == 'C!' ):
                r_numb = ptclObj_i.tagsDict["ring"]
                nring = self.calc_nnab(r_numb,self.ring_nbindex)
                ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
                if( debug ):
                    print ' linker found ',atom_i
                    print ' nring,ELCNT[6],ELCNT ',nring,ELCNT[6],ELCNT 
                if( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                    if( ELCNT[7] == 1 or ELCNT[16] == 1  ):      # Pyrrole or thiophene 
                        if( ELCNT[6] == 2  ):             #imidazole C2
                            N_o = self.bonded_nbindx[ pid_i ]
                            N_f = self.bonded_nbindx[ pid_i + 1 ] - 1
                            for indx_j in range( N_o,N_f+1):
                                pid_j = self.bonded_nblist[indx_j]
                                ptclObj_j =  self.ptclC[pid_j]
                                if( ptclObj_j.tagsDict["fftype"]== 'CA' ):
                                    if( debug ): print " changing ",ptclObj_j.tagsDict["fftype"], " to CS"
                                    ptclObj_j.tagsDict["fftype"]= 'CS'   
                                    if( update_chr ): ptclObj_i.charge =  0.22

        if(debug):   sys.exit('linkers')

        debug = 0
        if(debug):
            for pid_i, ptclObj_i  in self.ptclC:
                print pid_i,ptclObj_i.tagsDict["fftype"],ptclObj_i.tagsDict["number"]
            sys.exit('debug')

        # Check for unidentified atoms
        for pid_i, ptclObj_i  in self.ptclC:
            if ( ptclObj_i.tagsDict["fftype"] == '?' ):
                print ' atom ', pid_i , ptclObj_i.tagsDict["number"],' unknow '
                sys.exit(' Unknow atom ')

        return #(self)


    def biaryl_types(self, update_chr=False, debug = False):
        """
        Set biaryl types
        """
        #
        # label carbons 
        #

        if( debug): print "biaryl_types"
        n_rings = 0
        for pid_i, ptclObj_i  in self.ptclC:
            rnumb_i = ptclObj_i.tagsDict["ring"]
            if ( rnumb_i > n_rings ): n_rings = rnumb_i

        for r_numb in range( 1,n_rings+1 ):
            nring = self.calc_nnab(r_numb,self.ring_nbindex)
            ELCNT = self.calc_elcnt(r_numb,self.ring_nblist,self.ring_nbindex)

            if( debug ): print 'ring ',r_numb

            Nr_o = self.ring_nbindex[r_numb]
            Nr_f = self.ring_nbindex[r_numb+1] - 1
            for r_indx in range(Nr_o,Nr_f+1):
                pid_i = self.ring_nblist[r_indx]
                ptclObj_i =  self.ptclC[pid_i]
                NNAB_i = self.calc_nnab(pid_i,self.bonded_nbindx)
                ELCNT_i = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
                N_o = self.bonded_nbindx[ pid_i ]
                N_f = self.bonded_nbindx[ pid_i + 1 ] - 1
                NNAB_i_intra = 0
                for indx_j in range( N_o,N_f+1):
                    pid_j = self.bonded_nblist[indx_j]
                    ptclObj_j =  self.ptclC[pid_j]
                    if( ptclObj_j.tagsDict["ring"] == r_numb and ptclObj_j.tagsDict["number"] != 1 ):
                        NNAB_i_intra = NNAB_i_intra +  1

                if( debug ): print " atom ",pid_i,ptclObj_i.tagsDict["fftype"],r_numb,NNAB_i,NNAB_i_intra
                if(  NNAB_i_intra == 2 ):
                    # edge
                    if( ptclObj_i.tagsDict["number"] == 6  and NNAB_i == 3 ):
                        if( ELCNT_i[16] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                            ptclObj_i.tagsDict["fftype"] = 'CP'
                        elif( ELCNT_i[7] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                            ptclObj_i.tagsDict["fftype"] = 'CW'
                        elif( ELCNT_i[7] == 1 and ELCNT_i[6] == 1  ):
                            ptclObj_i.tagsDict["fftype"] = 'CW'
                        elif( ELCNT_i[8] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                            ptclObj_i.tagsDict["fftype"] = 'CW'
                        elif( ELCNT_i[16] == 1 and ELCNT_i[6] == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'CP'
                        elif( ELCNT_i[6] == 2 ):
                            ptclObj_i.tagsDict["fftype"] = 'CA'
                    if( ptclObj_i.tagsDict["number"] == 16  and NNAB_i == 2 and  ELCNT_i[6] == 2 ):
                        ptclObj_i.tagsDict["fftype"] = 'S'
                    if( ptclObj_i.tagsDict["number"] == 7 and  NNAB_i == 3 and  ELCNT_i[6] == 3 ):
                        ptclObj_i.tagsDict["fftype"] = 'NS'
                    if( ptclObj_i.tagsDict["number"] == 7 and  NNAB_i == 3 and  ELCNT_i[6] == 2 and  ELCNT_i[1] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'NA'

                if(  NNAB_i_intra == 3 ):
                    # fussed 
                    if( ptclObj_i.tagsDict["number"] == 6  and NNAB_i == 3 ):
                        ptclObj_i.tagsDict["fftype"] = 'CB'

                if( debug ):
                    print " atom ",atom_i,ptclObj_i.tagsDict["fftype"],NNAB_i,NNAB_i_intra

            for r_indx in range(Nr_o,Nr_f+1):
                # relable secondary atoms 
                pid_i = self.ring_nblist[r_indx]
                ptclObj_i =  self.ptclC[pid_i]
                NNAB_i = self.calc_nnab(pid_i,self.bonded_nbindx)
                ELCNT_i = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
                N_o = self.bonded_nbindx[ pid_i ]
                N_f = self.bonded_nbindx[ pid_i + 1 ] - 1
                NNAB_i_intra = 0
                CP_cnt = 0 
                CW_cnt = 0 
                for indx_j in range( N_o,N_f+1):
                    pid_j = self.bonded_nblist[indx_j]
                    ptclObj_j =  self.ptclC[pid_j]
                    if( ptclObj_j.tagsDict["ring"] == r_numb and  ptclObj_j.tagsDict["number"] != 1  ):
                        NNAB_i_intra = NNAB_i_intra +  1
                        if( ptclObj_j.tagsDict["fftype"]== "CP" ): CP_cnt = CP_cnt  + 1
                        if( ptclObj_j.tagsDict["fftype"]== "CW" ): CW_cnt = CW_cnt  + 1

                if(  NNAB_i_intra == 2 ):
                    if( ptclObj_i.tagsDict["number"] == 6  and NNAB_i == 3 ):
                        if( CP_cnt > 0  or  CW_cnt > 0 ):
                            ptclObj_i.tagsDict["fftype"] = 'CS'

                if( debug ):
                    print "    relable atom ",atom_i,ptclObj_i.tagsDict["fftype"],NNAB_i,NNAB_i_intra,CP_cnt,CW_cnt

        debug = 0         
        if( debug ):
            sys.exit('atom_types.biaryl_types')
        return ( self )


    def interring_types(self, update_chr=False ):
        """
        Set conjugated inter-ring carbons to type C! 
        """

        # Find Ring linkers
        debug = 0
        for pid_i, ptclObj_i  in self.ptclC:
            NNAB = self.calc_nnab(pid_i,self.bonded_nbindx)
            ELCNT = self.calc_elcnt(pid_i,self.bonded_nblist,self.bonded_nbindx)
            N_o = self.bonded_nbindx[ pid_i ]
            N_f = self.bonded_nbindx[ pid_i + 1 ] - 1
            if ( ptclObj_i.tagsDict["number"] == 6 ) :
                if( NNAB == 3 ):  # if sp2
                    for indx_j in range( N_o,N_f+1):
                        pid_j = self.bonded_nblist[indx_j]
                        ptclObj_j = self.ptclC[pid_j]
                        if ( ptclObj_j.tagsDict["number"] == 6 )  :
                            NNAB_j = self.calc_nnab(pid_j,self.bonded_nbindx)
                            if( NNAB_j == 3 ):  # if sp2
                                if(  ptclObj_i.tagsDict["ring"] !=  ptclObj_j.tagsDict["ring"] and ptclObj_i.tagsDict["ring"] != 0 and ptclObj_j.tagsDict["ring"] != 0): # ring linker
                                    if(debug):
                                        print atom_i,' and ',pid_j,' link ',ptclObj_i.tagsDict["fftype"],ATYPE[pid_j]
                                        print ' ', ptclObj_i.tagsDict["ring"] , ptclObj_j.tagsDict["ring"]
                                    ptclObj_i.tagsDict["fftype"] = 'C!'
                                    ptclObj_j.tagsDict["fftype"]= 'C!'

        return 


    def write_qgroup(self,debug=False):
        """
        Write charge group information

        """
        if( debug):
            print "max_qgroup_id",self.max_qgroup_id

        qgroup_qsum = np.zeros(self.max_qgroup_id+1)
        qgroup_resname = ["null"]*(self.max_qgroup_id+1)
        
        if( debug):
            print "len(qgroup_qsum)",len(qgroup_qsum)
            
        total_q = 0.0 
        for pid, ptclObj  in self.ptclC:
            if( debug ):
                print "ptclObj.tagsDict[qgroup]",ptclObj.tagsDict["qgroup"]
                 
            qgroup_qsum[ptclObj.tagsDict["qgroup"] ] += ptclObj.charge
            qgroup_resname[ptclObj.tagsDict["qgroup"] ] = ptclObj.tagsDict["resname"]
            total_q += ptclObj.charge

        if( debug):
            for qid in range(self.max_qgroup_id+1):
                print "Charge group {} with resname {} has a total charge of {} ".format(qid,qgroup_resname[qid],qgroup_qsum[qid])

            print " Total q {}".format(total_q)

        return 



    def rotate_xz(self,theta_xz,direction="counterclockwise",verbose=False):
        """

        Rotate around the y-axis in the xz plane
        
             |   cos theta_xz 0 -sin theta_xz   |
        Ry = |           0    1      0         |
             |_  sin theta_xz 0  cos theta_xz _|
      
        
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
            
        
        if( len(self.ptclC) > 0 ):
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
            

            if( len(self.ptclC[1].position ) == 3 ):
                    for pid_i, ptclObj_i  in self.ptclC:                        
                        ptclObj_i.position = Rxzdotv(ptclObj_i.position,cos_xz,sin_xz,sinprefix12,sinprefix21)

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

        
        if( len(self.ptclC) > 0 ):
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
            

            if( len(self.ptclC[1].position ) == 3 ):
                    for pid_i, ptclObj_i  in self.ptclC:    
                        ptclObj_i.position = Rxydotv(ptclObj_i.position,cos_xy,sin_xy,sinprefix12,sinprefix21)
                        
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

        
        if( len(self.ptclC) > 0 ):
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
            

            if( len(self.ptclC[1].position ) == 3 ):
                    for pid_i, ptclObj_i  in self.ptclC:    
                        ptclObj_i.position = Ryzdotv(ptclObj_i.position,cos_yz,sin_yz,sinprefix12,sinprefix21)
                        
        return


    def read_cply(self,cply_file):
        """
        Read cply file
        """
        read_lattice = False 

        # Load periodic table 
        pt = periodictable()
        
        with open(cply_file) as f:
            for line in f:
                col = line.split()
                if( read_lattice ):
                    self.latvec[lv_cnt][0] = float( col[0] )
                    self.latvec[lv_cnt][1] = float( col[1] )
                    self.latvec[lv_cnt][2] = float( col[2] )
                    if( lv_cnt == 2 ):
                        read_lattice = False 
                    lv_cnt += 1
                elif( len(col) >= 4  and col[0] != "bond" and col[0] != "#" ):
                    pt_i = Particle()
                    pt_i.type = str(col[0])
                    pt_i.position = [ float(col[1]),float(col[2]),float(col[3]) ]
                    add_dict = pt_i.tagsDict
                    el = pt.getelementWithSymbol(str(col[0]))
                    add_dict["symbol"] = str(col[0])
                    add_dict["number"] = el.number
                    add_dict["mass"] = el.mass
                    pt_i.mass = el.mass
                    add_dict["chain"] = 1
                    add_dict["ring"] = 0
                    add_dict["residue"] = 1
                    add_dict["resname"] = "RES"
                    add_dict["qgroup"] = 1
                    add_dict["fftype"] = "??"
                    add_dict["label"] =  el.symbol  
                    add_dict["cov_radii"] = el.cov_radii
                    add_dict["vdw_radii"] = el.vdw_radii
                    add_dict["lmptype"] = -1
                    if (len(col) >= 13 ):
                        add_dict["label"] = str(col[4])
                        add_dict["fftype"] = str(col[5])
                        add_dict["mass"] = float(col[6])                        
                        pt_i.mass = float(col[6])                   
                        pt_i.charge = float(col[7])
                        add_dict["qgroup"] = int(col[8])                        
                        add_dict["ring"] = int(col[9])                        
                        add_dict["residue"] = int(col[10])
                        add_dict["resname"] = str(col[11])                        
                        add_dict["chain"] = int(col[12])   
                    elif (len(col) == 7 ):
                        pt_i.charge = float(col[4])
                        add_dict["residue"] = int(col[5])
                        add_dict["resname"] = str(col[6])                                                              
                    pt_i.setTagsDict(add_dict)                    
                    self.ptclC.put(pt_i)

                elif(len(col) >= 3 ):
                    if( col[0] == "bond"):
                        b_i = int(col[1])
                        b_j = int(col[2])
                        bnd = Bond(b_i,b_j)
                        #print "process_line bond line ",col
                        self.bondC.put(bnd)
                # Key word search 
                if( len(col) > 0 ):
                    if( str(col[0]) == 'lattice' ):
                        read_lattice = True
                        lv_cnt = 0 
        return


    def write_cply(self,cply_file, write_ff=False,write_bonds=False,write_latvec=True):
        """
        Write cply file
        """

        print "writing cply",cply_file,len(self.ptclC)

        F = open(cply_file,'w')
        cply_line = "# atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]),label,fftype,ptclObj.mass,charge,qgroup,ring,residue,resname, chain, cplytag"
        F.write(cply_line)
        for pid, ptclObj  in self.ptclC:
            r_i = ptclObj.position
            atomic_symb = ptclObj.tagsDict['symbol']
            if(write_ff ):
                charge = ptclObj.charge
                residue = ptclObj.tagsDict["residue"]
                resname = ptclObj.tagsDict["resname"]
                label = ptclObj.tagsDict["label"]
                fftype = ptclObj.tagsDict["fftype"]
                qgroup = ptclObj.tagsDict["qgroup"]
                chain = ptclObj.tagsDict["chain"]
                ring = ptclObj.tagsDict["ring"]
                cply_line =  "\n %5s %16.8f %16.8f %16.8f %s %s %12.8f %12.8f  %d %d %d %s  %d "%(atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]),label,fftype,ptclObj.mass,charge,qgroup,ring,residue,resname, chain )
            else:
                cply_line =  "\n %5s %16.8f %16.8f %16.8f %s "%(atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]), cplytag )
            F.write(cply_line)

        if( write_latvec ):
            F.write("\n lattice")
            for d in range(3):
                cply_line = "\n %16.8f %16.8f %16.8f "%( self.latvec[d][0],self.latvec[d][1],self.latvec[d][2])
                F.write(cply_line)
            

        if( write_bonds ):
            for b_o, bondObj_o  in self.bondC:
                pid_i = bondObj_o.pgid1
                pid_j = bondObj_o.pgid2
                F.write("\n  bond %d %d  "%(pid_i,pid_j))
        F.close()
        
        return
