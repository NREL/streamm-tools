"""
Class data structures for 2, 3, 4 point groupings of Particle objects
"""
import copy, sys

class Dihedral:
    """
    Data structure for describing any 4-point associatiaon of Particle-s
    """

    def __init__(self, pgid1=0, pgid2=0, pgid3=0, pgid4=0, theta0=0.0, type="blank"):
        """
        Constructor for a general dihedral. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pgid1   (int)   GlobalID of Particle object in dihedral
            pgid2   (int)   GlobalID of Particle object in dihedral
            pgid3   (int)   GlobalID of Particle object in dihedral
            pgid4   (int)   GlobalID of Particle object in dihedral
            theta0  (float) Equilibrium dihedral (in radians)
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



class DihedralContainer:
    """
    Main data structure for holding Dihedral objects. Map of global
    dihedral ID (integer) to Dihedral object instances
    """

    def __init__(self, idList=[], verbose=False):
        """
        Constructor: sets up a dictionary for indexing 'Dihedral' objects

        Args:
            idList (list): of dihedral IDs. If empty then ID starts at 1.
                If not empty then ID's (keys) are inititalized with Dihedral objects
            verbose (bool): flag for printing status/debug info        
        """
        self.verbose=verbose
        self.dihedrals=dict()                              # Creates empty dict struc
        self.dihedrals={key: Dihedral() for key in idList} # Creates empty Dihedral objs
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
        del self.dihedrals
        del self.maxgid


    def clear(self):
        """
        Clears dihedrals out of DihedralContainer
        """
        self.maxgid = 0
        self.dihedrals=dict()                          # Creates empty dict struc


    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.dihedrals)

    def __str__(self):
        """
        'Magic' method for printng contents
        """

        dihedralStr="\n Contains dihedral objects: \n"
        for gid in self.dihedrals:
            dihedralStr = dihedralStr + str(gid) + " " + str(self.dihedrals[gid].__dict__) + "\n"
            bondStr += "%d %s %s %s \n"%(gid,self.dihedrals[gid].pgid1,self.dihedrals[gid].pgid2,self.dihedrals[gid].pgid3,self.dihedrals[gid].pgid4,self.dihedrals[gid].type )
            
        return dihedralStr


    def keys(self):
        """
        Return list of all ptcl IDs (keys) currently in container
        """
        keyList = self.dihedrals.keys()
        return keyList


    def __setitem__(self, gid, dihedral):
        """
        'Magic' method implementing obj[]=value operator
        Performs deep copy of value so container is managing memory
        """
        if gid in self.dihedrals.keys():
            self.dihedrals[gid]=copy.deepcopy(dihedral)
        else:
            print "Cannot add dihedral object to non-existent ID"
            sys.exit(3) 


    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator
        Operations on returned elements change container
        """
        return self.dihedrals[gid]


    def __delitem__(self, gid):
        """
        'Magic' method implementing del obj[] operator
        """
        del self.dihedrals[gid]


    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        return self.dihedrals.iteritems()


    def __call__(self, idSubList=None):
        """
        Callable magic method. Returns iterator to subset dihedrals dictionary

        Args:
             idSubList (list) list of pid-s of particle objects to be returned
             
        Returns: iterator to subset of particle dictionary
        """

        subGroupDct = dict()
        
        if idSubList != None:
            for gid, dihedralObj in self.dihedrals.iteritems():
                if gid in idSubList:
                    subGroupDct[gid] = dihedralObj
            return subGroupDct.iteritems()

        else:
            print "Callable DihedralContainer requires a list of subgroup dihedral IDs"
            sys.exit(3)


    def __contains__(self, gid):
        """
        'Magic' method implementing in keyword (key in obj')....
        """
        return gid in self.dihedrals


    def hasDihedral(self, dihedralList):
        """
        Check the ptcl IDs in dihedralList for any dihedral in container that is similar
        eg dihedral 1-2-3-4 is same as dihedral 4-3-2-1

        Args: (list) ptcl IDs defining dihedral to search for
        
        Returns: (bool) is dihedral in container
        """

        for gid, dObj in self.dihedrals.iteritems():
            dihedral = [dObj.pgid1, dObj.pgid2, dObj.pgid3, dObj.pgid4] # Dihedral ID list#
            dihedralRev = copy.deepcopy(dihedral)                       # Make reverse dihedral
            dihedralRev.reverse()                                       #  ID list

            if ( (dihedral == dihedralList) or (dihedralRev == dihedralList) ):
                return True
            
        return False


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator
        
        Compare global IDs of dihedrals and reassign globalIDs for dihedral
        container using the max ID between the two lists

        Note: for now this reassigns ID always
        """

        keys1 = self.dihedrals.keys()     # global IDs in this object
        keys2 = other.dihedrals.keys()    # global IDs in object being added
        bothkeys = keys1 + keys2          # List of all keys

        if len(bothkeys) > 0:                 # If keys not empty... proceed
            self.maxgid = max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        for ptclkey2 in other.dihedrals:
            self.put(other.dihedrals[ptclkey2])

        return self


    def put(self, dihedral):
        """
        Append 'Dihedral' object to this container. Updates globalID for container
        by incrementing the maxgid member

        Args:
            ptcl (Particle) correctly initialized Particle object

        NOTE:
            (1) One can imagine extra conditions on dihedrals inserted
            (2) This could check for uniqueness of all globalID's and throw error for copies
        """
        
        if isinstance(dihedral, Dihedral):
            self.maxgid += 1
            self.dihedrals[self.maxgid] = copy.deepcopy(dihedral)
        else:
            print "Attempting to add non-Dihedral type to container"
            raise TypeError


    def replacePtclIDs(self, idFromTo):
        """
        Replace ptclIDs given a dictionary of ID changes # eg {1:3, 3:5, 2:20...}
                
        Args:
            idFromTo (dict) map of ID changes
        """

        fromIDs = idFromTo.keys()
                
        for gid in self.dihedrals:
            
            dihedral = self.dihedrals[gid]  # Dihedral object
            pgid1 = dihedral.pgid1          # ptcl1 in dihedral
            pgid2 = dihedral.pgid2          # ptcl2 in dihedral
            pgid3 = dihedral.pgid3          # ptcl3 in dihedral
            pgid4 = dihedral.pgid4          # ptcl4 in dihedral
            
            if pgid1 in fromIDs:
                toID = idFromTo[pgid1]
                dihedral.pgid1 = toID

            if pgid2 in fromIDs:
                toID = idFromTo[pgid2]
                dihedral.pgid2 = toID

            if pgid3 in fromIDs:
                toID = idFromTo[pgid3]
                dihedral.pgid3 = toID

            if pgid4 in fromIDs:
                toID = idFromTo[pgid4]
                dihedral.pgid4 = toID



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
        for gid, dihedralObj in self.dihedrals.iteritems():
            dihedralType = dihedralObj.type
            typeList.append(dihedralType)
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
