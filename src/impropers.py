"""
Class data structures for 2, 3, 4 point groupings of Particle objects
"""
import copy, sys

class Improper:
    """
    Data structure for describing any 4-point associatiaon of Particle-s
    """

    def __init__(self, pgid1=0, pgid2=0, pgid3=0, pgid4=0, theta0=0.0, type="blank"):
        """
        Constructor for a general improper. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pgid1   (int)   GlobalID of Particle object in improper
            pgid2   (int)   GlobalID of Particle object in improper
            pgid3   (int)   GlobalID of Particle object in improper
            pgid4   (int)   GlobalID of Particle object in improper
            theta0  (float) Equilibrium improper (in radians)
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

        if isinstance(pgid4, int):
            self.pgid4 = pgid4
        else:
            print "4rd arg should be int type"
            raise TypeError

        if isinstance(theta0, float):
            self.theta0 = theta0
        else:
            print "5th arg should be float value"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "6th arg should be string value"
            raise TypeError


        self.lmpindx = 0
        self.g_indx = 0
        

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pgid1
        del self.pgid2 
        del self.pgid3 
        del self.pgid4
        del self.theta0
        del self.type
        del self.lmpindx
        del self.g_indx 

    def __contains__(self, pgid):
        """
        'Magic' method implementing 'in' keyword.

        Args:
            pgid (int) Particle GID to check against 'held' IDs
        """
        if ( (pgid == self.pgid1) or \
             (pgid == self.pgid2) or \
             (pgid == self.pgid3) or \
             (pgid == self.pgid4) ):
            return True
        else:
            return False


    def set_type(self,type):
        """
        Set bond type 
        """
        self.type = type
        
        
    def get_type(self):
        """
        Return bond type
        """
        return self.type


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



class ImproperContainer:
    """
    Main data structure for holding Improper Improper objects. Map of global
    improper ID (integer) to Improper object instances
    """

    def __init__(self, idList=[], verbose=False):
        """
        Constructor: sets up a dictionary for indexing 'Improper' objects

        Args:
            idList (list): of improper IDs. If empty then ID starts at 1.
                If not empty then ID's (keys) are inititalized with Improper objects
            verbose (bool): flag for printing status/debug info        
        """
        self.verbose=verbose
        self.impropers=dict()                              # Creates empty dict struc
        self.impropers={key: Improper() for key in idList} # Creates empty Improper objs
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
        del self.impropers
        del self.maxgid


    def clear(self):
        """
        Clears impropers out of ImproperContainer
        """
        self.maxgid = 0
        self.impropers=dict()                          # Creates empty dict struc


    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.impropers)

    def __str__(self):
        """
        'Magic' method for printng contents
        """

        improperstr="\n Contains improper objects: \n"
        for gid in self.impropers:
            improperstr = improperstr + str(gid) + " " + str(self.impropers[gid].__dict__) + "\n"
        return improperstr


    def keys(self):
        """
        Return list of all ptcl IDs (keys) currently in container
        """
        keyList = self.impropers.keys()
        return keyList


    def __setitem__(self, gid, improper):
        """
        'Magic' method implementing obj[]=value operator
        Performs deep copy of value so container is managing memory
        """
        if gid in self.impropers.keys():
            self.impropers[gid]=copy.deepcopy(improper)
        else:
            print "Cannot add improper object to non-existent ID"
            sys.exit(3) 


    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator
        Operations on returned elements change container
        """
        return self.impropers[gid]


    def __delitem__(self, gid):
        """
        'Magic' method implementing del obj[] operator
        """
        del self.impropers[gid]


    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        return self.impropers.iteritems()


    def __call__(self, idSubList=None):
        """
        Callable magic method. Returns iterator to subset impropers dictionary

        Args:
             idSubList (list) list of pid-s of particle objects to be returned
             
        Returns: iterator to subset of particle dictionary
        """

        subGroupDct = dict()
        
        if idSubList != None:
            for gid, improperObj in self.impropers.iteritems():
                if gid in idSubList:
                    subGroupDct[gid] = improperObj
            return subGroupDct.iteritems()

        else:
            print "Callable ImproperContainer requires a list of subgroup improper IDs"
            sys.exit(3)


    def __contains__(self, gid):
        """
        'Magic' method implementing in keyword (key in obj')....
        """
        return gid in self.impropers


    def hasImproper(self, improperList):
        """
        Check the ptcl IDs in improperList for any improper in container that is similar
        eg improper 1-2-3-4 is same as improper 4-3-2-1

        Args: (list) ptcl IDs defining improper to search for
        
        Returns: (bool) is improper in container
        """

        for gid, dObj in self.impropers.iteritems():
            improper = [dObj.pgid1, dObj.pgid2, dObj.pgid3, dObj.pgid4] # Improper ID list#
            improperRev = copy.deepcopy(improper)                       # Make reverse improper
            improperRev.reverse()                                       #  ID list

            if ( (improper == improperList) or (improperRev == improperList) ):
                return True
            
        return False


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator
        
        Compare global IDs of impropers and reassign globalIDs for improper
        container using the max ID between the two lists

        Note: for now this reassigns ID always
        """

        keys1 = self.impropers.keys()     # global IDs in this object
        keys2 = other.impropers.keys()    # global IDs in object being added
        bothkeys = keys1 + keys2          # List of all keys

        if len(bothkeys) > 0:                 # If keys not empty... proceed
            self.maxgid = max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        for ptclkey2 in other.impropers:
            self.put(other.impropers[ptclkey2])

        return self


    def put(self, improper):
        """
        Append 'Improper' object to this container. Updates globalID for container
        by incrementing the maxgid member

        Args:
            ptcl (Particle) correctly initialized Particle object

        NOTE:
            (1) One can imagine extra conditions on impropers inserted
            (2) This could check for uniqueness of all globalID's and throw error for copies
        """
        
        if isinstance(improper, Improper):
            self.maxgid += 1
            self.impropers[self.maxgid] = copy.deepcopy(improper)
        else:
            print "Attempting to add non-improper Improper type to container"
            raise TypeError


    def replacePtclIDs(self, idFromTo):
        """
        Replace ptclIDs given a dictionary of ID changes # eg {1:3, 3:5, 2:20...}
                
        Args:
            idFromTo (dict) map of ID changes
        """

        fromIDs = idFromTo.keys()
                
        for gid in self.impropers:
            
            improper = self.impropers[gid]  # Improper object
            pgid1 = improper.pgid1          # ptcl1 in improper
            pgid2 = improper.pgid2          # ptcl2 in improper
            pgid3 = improper.pgid3          # ptcl3 in improper
            pgid4 = improper.pgid4          # ptcl4 in improper
            
            if pgid1 in fromIDs:
                toID = idFromTo[pgid1]
                improper.pgid1 = toID

            if pgid2 in fromIDs:
                toID = idFromTo[pgid2]
                improper.pgid2 = toID

            if pgid3 in fromIDs:
                toID = idFromTo[pgid3]
                improper.pgid3 = toID

            if pgid4 in fromIDs:
                toID = idFromTo[pgid4]
                improper.pgid4 = toID


    # SWS: needs a test
    def getTypeInfoDict(self):
        """
        Return a map of type to typeIndex
        Method assigns a type index and checkes for consistency

        Returns:
            dictionary
        """

        # Look for types and get unique list
        typeList = list()
        for gid, improperObj in self.impropers.iteritems():
            improperType = improperObj.type
            typeList.append(improperType)
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
