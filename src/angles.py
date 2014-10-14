"""
Class data structures for 2, 3, 4 point groupings of Particle objects
"""
import copy, sys

class Angle:
    """
    Data structure for describing any 3-point associatiaon of Particle-s
    """

    def __init__(self, pgid1=0, pgid2=0, pgid3=0, theta0=0.0, type="blank"):
        """
        Constructor for a general angle. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pgid1   (int)   GlobalID of Particle object in angle
            pgid2   (int)   GlobalID of Particle object in angle
            pgid3   (int)   GlobalID of Particle object in angle
            theta0  (float) Equilibrium angle (in radians)
            type   (str)   Charge value in units of [e]
        """
        
        if isinstance(pgid1, int):
            self.pgid1 = pgid1
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(pgid2, int):
            self.pgid2 = pgid2
        else:
            print "2nd arg should be int type"
            raise TypeError

        if isinstance(pgid3, int):
            self.pgid3 = pgid3
        else:
            print "3rd arg should be int type"
            raise TypeError

        if isinstance(theta0, float):
            self.theta0 = theta0
        else:
            print "4th arg should be float value"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "5th arg should be string value"
            raise TypeError

        self.lmpindx = 0
        self.g_indx = 0
        
    def __del__(self):
        """
        Destructor, clears object memory
        """
        

    def __contains__(self, pgid):
        """
        'Magic' method implementing 'in' keyword.

        Args:
            pgid (int) Particle GID to check against 'held' IDs
        """
        if ( (pgid == self.pgid1) or \
             (pgid == self.pgid2) or \
             (pgid == self.pgid3) ):
            return True
        else:
            return False



    def set_lmpindx(self,lmpindx):
        """
        Set bond type index for lammps
        """
        self.lmpindx = lmpindx
        
        
    def get_lmpindx(self):
        """
        Return bond type index for lammps
        """
        return self.lmpindx

    def set_g_indx(self,g_indx):
        """
        Set bond type index for gromacs 
        """
        self.g_indx = g_indx
        
        
    def get_g_indx(self):
        """
        Return bond type index for gromacs
        """
        return self.g_indx


class AngleContainer:
    """
    Main data structure for holding Angle objects. Map of global
    angle ID (integer) to Angle object instances
    """

    def __init__(self, idList=[], verbose=False):
        """
        Constructor: sets up a dictionary for indexing 'Angle' objects

        Args:
            idList (list): of angle IDs. If empty then ID starts at 1.
                If not empty then ID's (keys) are inititalized with Angle objects
            verbose (bool): flag for printing status/debug info        
        """
        self.verbose=verbose
        self.angles=dict()                           # Creates empty dict struc
        self.angles={key: Angle() for key in idList} # Creates empty Angle objs
                                                     #  if idList not empty

        if len(idList) == 0:         # If list not set in constructor arg
            self.maxgid=0            # default=0 if idList empty
        else:                        #
            self.maxgid=max(idList)  # take max in list for maxgid


    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        if self.verbose:
            print "Cleaning particle container"
        del self.angles


    def clear(self):
        """
        Clears angles out of AngleContainer
        """
        self.maxgid = 0
        self.angles=dict()                          # Creates empty dict struc

    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.angles)

    def __str__(self):
        """
        'Magic' method for printng contents
        """

        angleStr="\n Contains angle objects: \n"
        for gid in self.angles:
            angleStr = angleStr + str(gid) + " " + str(self.angles[gid].__dict__) + "\n"
        return angleStr


    def keys(self):
        """
        Return list of all ptcl IDs (keys) currently in container
        """
        keyList = self.angles.keys()
        return keyList


    def __setitem__(self, gid, angle):
        """
        'Magic' method implementing obj[]=value operator
        Performs deep copy of value so container is managing memory
        """
        if gid in self.angles.keys():
            self.angles[gid]=copy.deepcopy(angle)
        else:
            print "Cannot add angle object to non-existent ID"
            sys.exit(3) 


    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator
        Operations on returned elements change container
        """
        return self.angles[gid]


    def __delitem__(self, gid):
        """
        'Magic' method implementing del obj[] operator
        """
        del self.angles[gid]


    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        return self.angles.iteritems()


    def __call__(self, idSubList=None):
        """
        Callable magic method. Returns iterator to subset angles dictionary

        Args:
             idSubList (list) list of pid-s of particle objects to be returned
             
        Returns: iterator to subset of particle dictionary
        """

        subGroupDct = dict()
        
        if idSubList != None:
            for gid, angleObj in self.angles.iteritems():
                if gid in idSubList:
                    subGroupDct[gid] = angleObj
            return subGroupDct.iteritems()

        else:
            print "Callable AngleContainer requires a list of subgroup angle IDs"
            sys.exit(3)


    def __contains__(self, gid):
        """
        'Magic' method implementing in keyword (key in obj')....
        """
        return gid in self.angles


    def hasAngle(self, angleList):
        """
        Check the ptcl IDs in angleList for any angle in container that is similar
        eg angle 1-2-3 is same as angle 3-2-1

        Args: (list) ptcl IDs defining angle to search for
        
        Returns: (bool) is angle in container
        """

        for gid, angleObj in self.angles.iteritems():
            angle = [angleObj.pgid1, angleObj.pgid2, angleObj.pgid3] # Angle ID list
            angleRev = copy.deepcopy(angle)                          # Make reverse angle
            angleRev.reverse()                                       #  ID list

            if ( (angle == angleList) or (angleRev == angleList) ):
                return True
            
        return False


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator
        
        Compare global IDs of angles and reassign globalIDs for angle
        container using the max ID between the two lists

        Note: for now this reassigns ID always
        """

        keys1 = self.angles.keys()        # global IDs in this object
        keys2 = other.angles.keys()       # global IDs in object being added
        bothkeys = keys1 + keys2          # List of all keys

        if len(bothkeys) > 0:                 # If keys not empty... proceed
            self.maxgid = max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        for ptclkey2 in other.angles:
            self.put(other.angles[ptclkey2])

        return self


    def put(self, angle):
        """
        Append 'Angle' object to this container. Updates globalID for container
        by incrementing the maxgid member

        Args:
            ptcl (Particle) correctly initialized Particle object

        NOTE:
            (1) One can imagine extra conditions on angles inserted
            (2) This could check for uniqueness of all globalID's and throw error for copies
        """
        
        if isinstance(angle, Angle):
            self.maxgid += 1
            self.angles[self.maxgid] = copy.deepcopy(angle)
        else:
            print "Attempting to add non-Angle type to container"
            raise TypeError


    def replacePtclIDs(self, idFromTo):
        """
        Replace ptclIDs given a dictionary of ID changes # eg {1:3, 3:5, 2:20...}
                
        Args:
            idFromTo (dict) map of ID changes
        """

        fromIDs = idFromTo.keys()
                
        for gid in self.angles:
            
            angle = self.angles[gid]  # Angle object
            pgid1 = angle.pgid1       # ptcl1 in angle
            pgid2 = angle.pgid2       # ptcl2 in angle
            pgid3 = angle.pgid3       # ptcl3 in angle
            
            if pgid1 in fromIDs:
                toID = idFromTo[pgid1]
                angle.pgid1 = toID

            if pgid2 in fromIDs:
                toID = idFromTo[pgid2]
                angle.pgid2 = toID

            if pgid3 in fromIDs:
                toID = idFromTo[pgid3]
                angle.pgid3 = toID



    def getTypeInfoDict(self):
        """
        Return a map of type to (typeIndex, ??????)
        Method assigns a type index and checkes for consistency

        NOTE: not tested

        Returns:
            dictionary of {type:[typeIndex, ????], ....}
        """

        # Look for types and get unique list
        typeList = list()
        for gid, angleObj in self.angles.iteritems():
            angleType = angleObj.type
            typeList.append(angleType)
        typeList = list(set(typeList))
        
        # Generate list of unique type for keys to initialize dictionary
        typeIndexDict = {key:index+1 for index, key in enumerate(typeList)} # dict for typeIndex
        if self.verbose:
            print "Unique types = ", typeList
            print "typeIndexDict = ", typeIndexDict

        # Pack yypeIndex to new dictionary
        typeInfoDict = dict()
        for key in typeIndexDict:
            index  = typeIndexDict[key]
            typeInfoDict[key] = index

        return typeInfoDict
