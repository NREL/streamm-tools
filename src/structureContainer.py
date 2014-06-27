"""
Class data structures for atomic data
"""

from particles import ParticleContainer
from bonds     import BondContainer

import copy

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
        self.boxLengths = [ [0.0, 1.0], [0.0, 1.0], [0.0,1.0] ]


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


    def replacePtclIDs(self, findPtclID, newPtclID):
        """
        Replace IDs that contain globalID of particle
        Driver for each container's replace ID method

        Args:
            findPtclID (int) globalID of particle to search for
            newPtclID  (int) globalID to replace findPtclID with
        """

        self.ptclC.replaceID(     findPtclID, newPtclID)
        self.bondC.replacePtclIDs(findPtclID, newPtclID)


    def getSubStructure(self, ptclIDList):
        """
        Return a new Structure object with partcleID's in input list
        
        Args:
            ptclIDList (list) global particles ID's for which to return structure

        Return:
            New Structure() object. IDs in new object are unique
        """

        subAtoms = ParticleContainer()
        subBonds = BondContainer()

        for pgid in ptclIDList:
            atom = self.ptclC[pgid]
            subAtoms.put(atom)

        return StructureContainer(subAtoms, subBonds)


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


    def dumpLammpsInputFile(self, inputName, coeffDict=dict()):
        """
        Write out a LAMMPS input data file from all available held
        data (particles, bonds, angles, dihedrals)

        Args:
            inputName (str)  name of LAMMPS input file to write
            coeffDict (dict) dictionary of potential parameters eg...
                      coeffDict = {("Si", "epsilon"):2.30, ("Si", "sigma"):1.0,
                                   ("C",  "epsilon"):0.50, ("C",  "sigma"): 0.1 }
        """

        if not isinstance(coeffDict, dict):
            print "dumpLammpsInputFile: coeffDict should be a python dictionary"
            sys.exit(3)

        n_atoms = len(self.ptclC)  # Obtaining particle size from container
        n_bonds = len(self.bondC)  # Need to edit
        n_angles = 0
        n_dihedrals = 0
        n_impropers = 0

        typeInfoDict = self.ptclC.getTypeInfoDict()  # map of "type":[typeIndex, mass, charge]
        typeList = typeInfoDict.keys()               # list of types eg ["Si", "C", ..]

        # Returns map of type,parameter tuple and value
        # SWS: particular to this method
        # typeParams = self.ptclC.getTypeParams()
        
        n_atypes = len(typeInfoDict)
        n_btypes = 0
        n_angtypes = 0 
        n_dtypes = 0 
        n_imptypes = 0

        xL = self.boxLengths[0]
        yL = self.boxLengths[1]
        zL = self.boxLengths[2]

        # Open file, write header info
        fileObj = open( inputName, 'w' )
        fileObj.write('Lammps data file \n')
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
        for type, info in typeInfoDict.iteritems():
            tIndex = info[0]
            mass   = info[1]
            fileObj.write( massFormatStr % ( tIndex, mass ) )
        fileObj.write('\n')
        
        nonCoeffFormatStr = "%5d %12.6f %12.6f  \n"
        fileObj.write('Nonbond Coeffs \n')
        fileObj.write('\n')

        for typ in typeList:
            info = typeInfoDict[typ]  # map of "type":[typeIndex, mass, charge]
            tIndex = info[0]          # type index for 'typ' (eg "Si")
            epsilon = coeffDict[(typ, "epsilon")]
            sigma   = coeffDict[(typ, "sigma")]
            fileObj.write( nonCoeffFormatStr % (tIndex, epsilon, sigma  ) )
        fileObj.write('\n')

        
        ptclFormatStr = "%5d %5d %5d %12.8f %12.6f %12.6f %12.6f \n"
        fileObj.write('Atoms \n')
        fileObj.write('\n')
        for pid, ptclObj in self.ptclC:
            pos = ptclObj.position
            mol = ptclObj.tagsDict["molnum"]
            mol = int(mol)
            typ = ptclObj.type
            typeIndex = typeInfoDict[typ][0]
            chg = ptclObj.charge
            fileObj.write( ptclFormatStr % (pid, mol, typeIndex, chg, pos[0], pos[1], pos[2] ) )
        fileObj.write('\n')

        # Close LAMMPS file
        fileObj.close()

        
        # Print type mapping info
        fileObj = open( inputName + ".dat", 'w' )
        fileObj.write("type   typeIndex    mass   charge \n")
        for type, info in typeInfoDict.iteritems():
            datFormatStr = "%s %5d %12.6f %12.6f \n"
            tIndex = info[0]
            mass   = info[1]
            charge = info[2]
            fileObj.write(datFormatStr % ( type, tIndex, mass, charge ) )         

        # Close LAMMPS mapping file
        fileObj.close()

    def putstruc_json(self, json_data, json_file):
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
        struc_data["structure"]["particle"] = particle_data
        struc_data["structure"]["twobody"] = twobody_data
        struc_data["structure"]["threebody"] = threebody_data
        struc_data["structure"]["fourbody"] = fourbody_data

        struc_data["number_id"] = []
        struc_data["position"] = []
        struc_data["mass"] = []
        struc_data["charge"] = []
        struc_data["chain"] = []
        struc_data["ring"] = []
        struc_data["resname"] = []
        struc_data["residue"] = []
	        
        for  pid, ptclObj in self.ptclC:
            struc_data["number_id"].append( ptclObj.position )
            struc_data["position"].append( ptclObj. )
            struc_data["mass"].append( ptclObj.mass )
            struc_data["charge"].append( ptclObj.charge )
            struc_data["chain"].append( ptclObj.chain )
            struc_data["ring"].append( ptclObj.ring )
            struc_data["resname"].append( ptclObj.resname )
            struc_data["residue"].append( ptclObj.residue )
        
            

def getstruc_json(json_file):
        """
        Return a new Structure object with partcleID's in input list
        
        Args:
            ptclIDList (list) global particles ID's for which to return structure

        Return:
            New Structure() object. IDs in new object are unique
        """

        subAtoms = ParticleContainer()
        subBonds = BondContainer()

        for pgid in ptclIDList:
            atom = self.ptclC[pgid]
            subAtoms.put(atom)

        return StructureContainer(subAtoms, subBonds)
