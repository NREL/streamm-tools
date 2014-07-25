"""
Class data structures for atomic data
"""

from particles import Particle
from particles import ParticleContainer
from bonds     import Bond
from bonds     import BondContainer
import pbcs

import copy
import numpy as np 
import json

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

    # def __init__(self, ptctC, bondC, angleC, dihC):
    def __init__(self, ptclC=ParticleContainer(), bondC=BondContainer()):
        """
        Constructor for a composite structure. Deepcopy of containers is used
        so this is the main container that manages memory for all sim objects

        Args:
            ptclC  (ParticleContainer)  
            bondC  (BondContainer)  
            angleC (AngleContainer)  * not implemented
            dihC   (DihedralContainer) * not implemented
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

        # Length of cartesian box size [ [xmin, xmax], [ymin, ymax]...]
        self.boxLengths = [ [0.0, 1.0], [0.0, 1.0], [0.0, 1.0] ]

        # Lattice vectors 
        self.latticevec = [  np.array([100.0,0.0,0.0]) ,  np.array( [0.0,100.0,0.0]),  np.array( [0.0,0.0,100.0]) ]


    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        print "Cleaning structureContainer"
        del self.ptclC
        del self.bondC


    def dump(self, filePrefix):
        """
        Dump a pickled version of this object

        Args:
            filePrefix (str): name of pickle file. will dump filePrefix.pkl
        """
        
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
        self.ptclC = copy.deepcopy(struc.ptclC)
        self.bondC = copy.deepcopy(struc.bondC)
        self.boxLengths = copy.deepcopy(struc.boxLengths)
        self.latticevec = copy.deepcopy(struc.latticevec)


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
        return strucStr


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator
        
        Compare global IDs of particles and reassign globalIDs for particle
        container using the max ID between the two lists. Tracks these changes
        for all other (bond, angle, dihedral) containers that reference particleIDs
        """

        idFromToDict = dict()  # Need to keep track of all ptcl ID changes at once
                               # {fromID1:toID1, fromID2:toID2...}
                               # eg {1:3, 3:5, 2:20...}
        
        bondC = BondContainer()            # Local bond container copy so ptclIDs
        bondC = copy.deepcopy(other.bondC) # inside can be changed (for adding below)

        keys1 = self.ptclC.particles.keys()    # global IDs of particles in this object
        keys2 = other.ptclC.particles.keys()   # global IDs in object being added
        self.ptclC.maxgid= max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        for ptclkey2 in other.ptclC.particles:
            self.ptclC.put(other.ptclC.particles[ptclkey2]) # Pushes ptcl to this struc's ptcl container
            fromPtclID = ptclkey2                           # Track IDs from--->to
            toPtclID   = self.ptclC.maxgid                  #  --> toID (to is the maxid of this ptclC)
            idFromToDict[fromPtclID]=toPtclID               # Store ID changes

        bondC.replacePtclIDs(idFromToDict)              # Use tracked list of ID changes
        self.bondC += bondC                             # Now add bondC with 'corrected' IDs

        # self.angleC += other.angleC                # TBI
        # angleC.replacePtclIDs(idFromToDict)        # TBI

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

        for toPtclID, ptclTuple in enumerate(self.ptclC):     # Enumerate returns (ID, obj) tuple for ptclTuple
            toPtclID +=1                                      # Sets reordering index correctly
            fromPtclID = ptclTuple[0]                         # Picks out ID from ptclTuple
            idFromToDict[fromPtclID]=toPtclID                 # Store ID changes
            ptclObj = self.ptclC.particles.pop(fromPtclID)    # Remove old ID
            self.ptclC.particles[toPtclID] = ptclObj          # reassign ptcl obj as new ID

        self.bondC.replacePtclIDs(idFromToDict)           # Use tracked list of ID changes
        # self.angleC.replacePtclIDs(idFromToDict)        # TBI

        for toBondID, bondTuple in enumerate(self.bondC):   # Enumerate returns (ID, obj) tuple for ptclTuple
            toBondID +=1                                    # Sets reordering index correctly
            fromBondID = bondTuple[0]                       # Picks out ID from ptclTuple
            bondObj = self.bondC.bonds.pop(fromBondID)      # Remove old ID
            self.bondC.bonds[toBondID] = bondObj            # reassign bond obj as new ID



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

        subAtoms = ParticleContainer(ptclIDList) # Initial ptcl container w/input IDs
        bondIDList = self.bondC.keys()           # Get keys of bond container
        subBonds = BondContainer(bondIDList)     # Intitialize subbond container

        # Grab particles from IDlist and put into sub-particle container
        for pgid in ptclIDList:
            atom = self.ptclC[pgid]
            subAtoms[pgid] = atom

        # For each bond object in container check that both
        # particles in bond are in ptcl search list
        for gid, bondObj in self.bondC:
            if ( (bondObj.pgid1 in ptclIDList) and (bondObj.pgid2 in ptclIDList) ):
                subBonds[gid] = bondObj
            else:
                # Need to remove empty key generated above
                del subBonds[gid]
                
        return StructureContainer(subAtoms, subBonds)




    def getpartnumb(self):
        """
        Return number of particles in a structure 
        """
        NP = 0
        for pid, ptclObj in self.ptclC :
            NP += 1

        return NP

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
        
        return self.latticevec

    def getVolume(self):
        """
        Calculate volume

        Method:
            Volume = ( v_i x v_j ) \dot v_k
        """

        br1 = np.cross(self.latticevec[0],self.latticevec[1])
        vol = np.dot(br1,self.latticevec[2])

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
                r_ij_sq = pbcs.sq_drij_c(r_i,r_j,self.latticevec)
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
        print "    v_i ",self.latticevec[0]
        print "    v_j ",self.latticevec[1]
        print "    v_k ",self.latticevec[2]
        print "  Bonds %d "%(len(self.bondC))

        
    def dumpLammpsInputFile(self, inputName, pairCoeffDct=dict(), bondCoeffDct=dict() ):
        """
        Write out a LAMMPS input data file from all available held
        data (particles, bonds, angles, dihedrals)

        Args:
            inputName    (str)  name of LAMMPS input file to write
            pairCoeffDct (dict) dictionary of potential parameters eg...
                                {("Si", "epsilon"):2.30, ("Si", "sigma"):1.0, ("C",  "epsilon"):0.50, ("C",  "sigma"): 0.1 }
            bondCoeffDct (dict) ""
        """

        if not isinstance(pairCoeffDct, dict):
            print "dumpLammpsInputFile: pairCoeffDct should be a python dictionary"
            sys.exit(3)

        if not isinstance(bondCoeffDct, dict):
            print "dumpLammpsInputFile: pairCoeffDct should be a python dictionary"
            sys.exit(3)

        n_atoms = len(self.ptclC)  # Obtaining particle size from container
        n_bonds = len(self.bondC)  # " "
        n_angles = 0
        n_dihedrals = 0
        n_impropers = 0

        ptclTypeInfo = self.ptclC.getTypeInfoDict()  # map of "type":[typeIndex, mass, charge]
        bondTypeInfo = self.bondC.getTypeInfoDict()  # map of "type":typeIndex

        # Returns map of type,parameter tuple and value
        # SWS: particular to this method
        n_atypes = len(ptclTypeInfo)
        n_btypes = len(bondTypeInfo)
        n_angtypes = 0 
        n_dtypes = 0 
        n_imptypes = 0

        xL = self.boxLengths[0]
        yL = self.boxLengths[1]
        zL = self.boxLengths[2]

        # Open file, write header info
        fileObj = open( inputName, 'w' )
        fileObj.write('LAMMPS Data File \n')
        fileObj.write('\n')
        fileObj.write( "%8d  atoms \n" % n_atoms )
        fileObj.write( "%8d  bonds \n" %  n_bonds )
        fileObj.write( "%8d  angles \n" % n_angles )
        fileObj.write( "%8d  dihedrals \n" %  n_dihedrals )
        fileObj.write( "%8d  impropers \n" % n_impropers  )
        fileObj.write('\n')
        fileObj.write( "%8d  atom types \n" % n_atypes  )
        fileObj.write( "%8d  bond types \n" % n_btypes )
        fileObj.write( "%8d  angle types \n" % n_angtypes )
        fileObj.write( "%8d  dihedral types \n" % n_dtypes )
        fileObj.write( "%8d  improper types \n" % n_imptypes )
        fileObj.write('\n')
        fileObj.write( "%16.8f %16.8f   xlo xhi \n" %  (xL[0] , xL[1] ) )
        fileObj.write( "%16.8f %16.8f   ylo yhi \n" %  (yL[0] , yL[1] ) )
        fileObj.write( "%16.8f %16.8f   zlo zhi \n" %  (zL[0] , zL[1] ) )
        fileObj.write('\n')

        massFormatStr = "%5d %16.8f \n"
        fileObj.write('Masses \n')
        fileObj.write('\n')
        for type, info in ptclTypeInfo.iteritems():
            tIndex = info[0]
            mass   = info[1]
            fileObj.write( massFormatStr % ( tIndex, mass ) )
        fileObj.write('\n')
        
        pairCoeffFormatStr = "%5d %12.6f %12.6f  \n"
        fileObj.write('Pair Coeffs \n')
        fileObj.write('\n')
        for typ in ptclTypeInfo.keys(): # list of types eg ["Si", "C", ..]
            info = ptclTypeInfo[typ]    # map of "type":[typeIndex, mass, charge]
            tIndex = info[0]            # type index for 'typ' (eg "Si")
            epsilon = pairCoeffDct[(typ, "epsilon")]
            sigma   = pairCoeffDct[(typ, "sigma")]
            fileObj.write( pairCoeffFormatStr % (tIndex, epsilon, sigma  ) )
        fileObj.write('\n')

        bondCoeffFormatStr = "%10d %12.6f %12.6f \n"
        if (n_bonds > 0):
            fileObj.write('Bond Coeffs \n')
            fileObj.write('\n')
            for typ in bondTypeInfo.keys(): # list of types eg ["Si", "C", ..]
                tIndex  = bondTypeInfo[typ]    # map of "type":[typeIndex, mass, charge]
                kenergy = bondCoeffDct[(typ, "Kenergy")]
                r0      = bondCoeffDct[(typ, "r0")]
                fileObj.write( bondCoeffFormatStr % (tIndex, kenergy, r0) )
            fileObj.write('\n')

        ptclFormatStr = "%5d %5d %5d %12.8f %12.6f %12.6f %12.6f \n"
        fileObj.write('Atoms \n')
        fileObj.write('\n')
        for pid, ptclObj in self.ptclC:
            pos = ptclObj.position
            mol = ptclObj.tagsDict["molnum"]
            mol = int(mol)
            typ = ptclObj.type
            typeIndex = ptclTypeInfo[typ][0]
            chg = ptclObj.charge
            fileObj.write( ptclFormatStr % (pid, mol, typeIndex, chg, pos[0], pos[1], pos[2] ) )
        fileObj.write('\n')
        
        bondFormatStr = "%9d %8d %9d %9d \n"
        if (n_bonds > 0):
            fileObj.write('Bonds \n')
            fileObj.write('\n')
        for gid, bondObj in self.bondC:
            pid1 = bondObj.pgid1
            pid2 = bondObj.pgid2
            bondType = 1
            fileObj.write( bondFormatStr % (gid, bondType, pid1, pid2) )
        if (n_bonds > 0):            
            fileObj.write('\n')
        
        # Close LAMMPS file
        fileObj.close()


        #
        # Print type mapping info
        #
        fileObj = open( inputName + ".dat", 'w' )
        fileObj.write("type   typeIndex    mass   charge \n")
        for type, info in ptclTypeInfo.iteritems():
            datFormatStr = "%s %5d %12.6f %12.6f \n"
            tIndex = info[0]
            mass   = info[1]
            charge = info[2]
            fileObj.write(datFormatStr % ( type, tIndex, mass, charge ) )         

        # Close LAMMPS mapping file
        fileObj.close()

    def putstruc_json(self, json_data ):
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
        lv_string = str( "%f %f %f %f %f %f %f %f %f " % ( self.latticevec[0][0], self.latticevec[0][1], self.latticevec[0][2], self.latticevec[1][0], self.latticevec[1][1], self.latticevec[1][2], self.latticevec[2][0], self.latticevec[2][1], self.latticevec[2][2]))

	struc_data["latticevector"] = lv_string

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


    def getsys_json(self, json_file):
        """
        Read in structure information from json file 

        Args:
            json_file (json) file with structure information 

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
            chain_i = particle_data["chain"][p_i]
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
            b_i = Bond( a_i+1, a_j+1 )            
            self.bondC.put(b_i)

        # Read in lattice vectors
        self.latticevec = []
        lv_array = struc_data["latticevector"].split()
        self.latticevec.append(  np.array( [float(lv_array[0]),float(lv_array[1]),float(lv_array[2])] ) )
        self.latticevec.append(  np.array( [float(lv_array[3]),float(lv_array[4]),float(lv_array[5])] ) )
        self.latticevec.append(  np.array( [float(lv_array[6]),float(lv_array[7]),float(lv_array[8])] ) )

    def create_top(self,ff_charges):
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
        
        for pid, ptclObj  in self.ptclC:
            ASYMB.append( ptclObj.type  )
            R.append( np.array( ptclObj.position)  )
            AMASS.append( ptclObj.mass  )
            CHARGES.append( ptclObj.charge  )
            MOLNUMB.append( ptclObj.tagsDict["chain"]  )
            RESID.append( ptclObj.tagsDict["residue"]  )
            RESN.append( ptclObj.tagsDict["resname"]  )
            
            
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
        GTYPE = top.initialize_gtype( ELN )
        CHARN = top.initialize_charn( ELN )        

        #   Build covalent nieghbor list for bonded information 
        NBLIST, NBINDEX = top.build_covnablist(ELN,R)

        NA = len(ELN)
        BONDS = top.nblist_bonds(NA,NBLIST, NBINDEX)
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
        RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)

        # Asign oplsaa atom types
        ATYPE, CHARGES = atom_types.oplsaa( ff_charges,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB)

        ATYPE , CHARGES = atom_types.biaryl_types( ff_charges, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES )
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
        for i in range( len(BONDS) ):
            #
            a_i = BONDS[i][0] 
            a_j = BONDS[i][1]
            b_i = Bond( a_i+1, a_j+1 )
            self.bondC.put(b_i)

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


    def lmp_writedata(self,data_file,norm_dihparam):
        """
        Write out lammps data file
        """

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
            
        LV[0][0] = 200.0
        LV[1][1] = 200.0
        LV[2][2] = 200.0
        
        # Find atomic number based on atomic symbol 
        ELN = elements.asymb_eln(ASYMB)
        GTYPE = top.initialize_gtype( ELN )
        NA = len(ELN)
        
        # Create neighbor list form bonds
        NBLIST,NBINDEX = groups.build_nablist_bonds(ELN,BONDS)
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
        #DIH = top.nblist_dih(NA,NBLIST, NBINDEX,options.limdih,options.limitdih_n)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)

        # Read in parameter files 
        itp_file = "oplsaa.itp"
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

        
    def calc_rdf(self, rdf_cnt_ij,bin_size,list_i,list_j,sq_r_cut):
        """
        Calculate RDF for a group of particles
        """

        #
        # Loop over list i
        #
        for p_i, ptcl_i in self.ptclC(list_i):
            r_i = np.array( ptcl_i.position )
            #
            # Loop over list j
            #
            for p_j, ptcl_j in self.ptclC(list_j):
                
                if( p_j > p_i):
                    r_j =  np.array( ptcl_j.position )
                    r_ij_sq = pbcs.sq_drij_c(r_i,r_j,self.latticevec)
                    if( r_ij_sq <= sq_r_cut ):
                        m_ij = np.sqrt(r_ij_sq)
                        bin_index = int( round( m_ij/bin_size) )
                        rdf_cnt_ij[bin_index] += 2

        return rdf_cnt_ij
    
