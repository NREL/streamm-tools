"""
Class data structures for atomic data
"""
import copy, sys

class Particle:
    """
    Data structure for describing any localized object in QM/MD simulations
    A 'Particle' has a position and type specifiers

    Note, for an object p1 of Particle type
      p1.__dict__ returns a dictionary of all class members

    Current use cases include: atoms, united-atoms and 'beads'

    Possible tags:
    ELN = []
    ASYMB = []
    CTYPE = []
    CHARGES = []
    R = []
    VEL = []
    ATYPE = []
    AMASS = []
    MOLNUMB = []
    RING_NUMB = []
    RESID = []
    RESN = []
    CHARN = []
    UNITNUMB = [] 
    UNITTYPE = []
    """

    def __init__(self, pos=[], type="blank", charge=0.0, mass=1.0):
        """
        Constructor for a general particle. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pos     (list)  Cartesian position vector
            type     (str)  String identifier (eg, Si, UA1, Bead1)
            charge (float)  Charge value in units of [e]
            mass   (float)  Mass in units of [?]
        """
        
        if isinstance(pos, list):
            self.position = pos
        else:
            print "1st arg should be list of floats"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "2nd arg should be string type"
            raise TypeError

        if isinstance(charge, float):
            self.charge = charge
        else:
            print "3rd arg should be float value"
            raise TypeError

        if isinstance(mass, float):
            self.mass = mass
        else:
            print "4th arg should be float value"
            raise TypeError

        # Tags dictionary. To be set by caller
        self.tagsDict=dict()
        self.tagsDict["type"] = type
        

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.tagsDict
        del self.position


    def setTagsDict(self, td):
        """
        Define custom dictionary for id tags (eg {"molnum":1,})
        Leaving type tag as a default

        Args:
            td (dictionary)
        """
        newdict = copy.deepcopy(td)
        type = self.tagsDict["type"]
        self.tagsDict=dict()
        self.tagsDict["type"] = type
        self.tagsDict.update(newdict)


    def isTagEqualTo(self, key, valList):
        """
        Check if tags[key]= any of the values in valList. Exits if key not in dictionary

        Args:
            key            dictionary key
            valList (list) dictionary values
            
        Returns: true if tagsDict[key]=value for any value in valueList
        """
        valueList = list()

        # If valueList not a list... convert to one for logic below
        if not isinstance(valList, list):
            valueList.append(valList)
        else:
            valueList = valList

        tagKeys = self.tagsDict.keys()
        if key not in tagKeys:
            print "Tag key in Particle dictionary does not exist"
            sys.exit(3) 

        for value in valueList:                 # Look at each search value
            if (self.tagsDict[key]==value):     # Determine if key,value match
                return True
     
        # If loop above has not returned with a found (True) state
        # then return with False (meaning not ANY value in valueList was present
        return False




class ParticleContainer:
    """
    Main data structure for holding Particle objects. Map of global
    particle ID (integer) to Particle object instances
    """

    def __init__(self, idList=[], verbose=False):
        """
        Constructor: sets up a dictionary for indexing 'Particle' objects

        Args:
            idList (list): of particle IDs. If empty then ID starts at 1.
                If not empty then ID's (keys) are inititalized with empty particle objects
            verbose (bool): flag for printing status/debug info
        """
        self.verbose=verbose
        self.particles=dict()                              # Creates empty dict struc
        self.particles={key: Particle() for key in idList} # Creates empty Particle objs
                                                           #  if idList not empty

        if len(idList) == 0:         # If list not set in constructor arg
            self.maxgid=0            # default=0 if idList empty
        else:                        #
            self.maxgid=max(idList)  # take max in list for maxgid


    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        print "Cleaning particle container"
        del self.particles


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



    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.particles)


    def __str__(self):
        """
        'Magic' method for printng contents of container
        """

        ptclStr="\n Contains particle objects: \n"
        for gid in self.particles:
            ptclStr = ptclStr + str(gid) + " " + str(self.particles[gid].__dict__) + "\n"
        return ptclStr


    def __setitem__(self, gid, ptcl):
        """
        'Magic' method implementing obj[]=value operator
        Performs deep copy of value so container is managing memory
        If gid exists in container then particle object is overwritten.
        If gid does not exist then particle object is inserted in container
        using the 'put()' method.
        """
        if gid in self.particles.keys():
            self.particles[gid]=copy.deepcopy(ptcl)
        else:
            print "Using [] operator for non-existent key"
            sys.exit(3)


    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator.
        Operations on returned elements change container
        """
        return self.particles[gid]


    def __delitem__(self, gid):
        """
        'Magic' method implementing del obj[] operator
        """
        del self.particles[gid]


    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        # return iter(self.particles)
        return self.particles.iteritems()


    def __call__(self, idSubList=None):
        """
        Callable magic method. Returns iterator to subset particle dictionary

        Args:
             idSubList (list) list of pid-s of particle objects to be returned
             
        Returns: iterator to subset of particle dictionary
        """

        subGroupDct = dict()
        
        if idSubList != None:
            for pid, ptclObj in self.particles.iteritems():
                if pid in idSubList:
                    subGroupDct[pid] = ptclObj
            return subGroupDct.iteritems()
        
        else:
            print "Callable ParticleContainer requires a list of subgroup particle IDs"
            sys.exit(3)


    def __contains__(self, gid):
        """
        'Magic' method implementing in keyword (key in obj')....
        """
        return gid in self.particles


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator
        
        Compare global IDs of particles and reassign globalIDs for particle
        container using the max ID between the two lists

        Note: for now this reassigns ID always
        """

        keys1 = self.particles.keys()     # global IDs in this object
        keys2 = other.particles.keys()    # global IDs in object being added
        self.maxgid = max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        for ptclkey2 in other.particles:
            self.put(other.particles[ptclkey2])
            
        return self


    def put(self, ptcl):
        """
        Append 'Particle' object to this container. Updates globalID for container
        by incrementing the maxgid member. Performs deep copy of ptcl so container
        is managing memory

        Args:
            ptcl (Particle) correctly initialized Particle object

        NOTE:
            (1) One can imagine extra conditions on distances between ptcl being
                inserted and current particles.
            (2) This could check for uniqueness of all globalID's and throw error
                for copies
        """
        
        if isinstance(ptcl, Particle):
            self.maxgid += 1
            self.particles[self.maxgid] = copy.deepcopy(ptcl)
        else:
            print "Attempting to add non-Particle type to container"
            raise TypeError


    def replaceID(self, findID, newID):
        """
        Replace global particle ID with another value. Attempts to make rational
        checks about new IDs

        Args:
            findGID (int) ptcl ID to search
            newGID  (int) ptcl ID to replace findGID with
        Returns:
            Errors if
              - findGID not found
              - newGID already present
        """

        keys = self.particles.keys()
        
        if findID not in keys:
            print "findID not found"
            sys.exit(3)

        if newID in keys:
            print "newID already exists"
            sys.exit(3)

        if newID < self.maxgid:
            # Pop off value for old key and reset with new key
            self.particles[newID] = self.particles.pop(findID)
            
        elif newID == self.maxgid+1:
            # If newID happens to be next ID
            self.particles[newID] = self.particles.pop(findID)
            self.maxgid += 1
        else:
            print "Error: ParticleContainer:replaceID newID >= maxgid \n"
            sys.exit(3)


    def getParticlesWithTags(self, searchDict):
        """
        Return particles whose dictionary tags match all specified key,values
        The search dictionary criteria is that ALL key,value pairs must match
        (the logic is {key1:value1} AND {key2:value2}.....

        Args:
            searchDict (dict) dictionary of search pairs {keyTag:valTag}....
            
        Returns: list of particle gid's that satisfy search
        """

        subGroupGID = list()
        tagFoundBoolList = list()
        if self.verbose:
            print "search = ", searchDict
            
        for pid, ptclObj in self.particles.iteritems():            # Iterate on particles
            for searchTag, searchVal in searchDict.iteritems():    # Iterate on dict items to search
                
                tagEqual = ptclObj.isTagEqualTo(searchTag, searchVal)  # Check if dct[searchTag]=searchVal
                tagFoundBoolList.append(tagEqual)                      # Store list of all results

            if False not in tagFoundBoolList:    # If only True-s in list then all conditions met (AND condition)
                subGroupGID.append(pid)          # so store pid of particle
            tagFoundBoolList = list()            # Reset tag found list


        return subGroupGID
        

    def getTypeInfoDict(self):
        """
        Return a map of type to (typeIndex, mass, charge)
        Method assigns a type index and checkes for consistency

        Returns:
            dictionary of {type:[typeIndex, mass, charge], ....}
        """

        # Look for types and get unique list
        typeList = list()
        for pid, ptclObj in self.particles.iteritems():
            ptclType = ptclObj.type
            typeList.append(ptclType)
        typeList = list(set(typeList))
        
        # Generate list of unique type for keys to initialize dictionary
        typeMassDict  ={key: list() for key in typeList}                    # dict for mass
        typeChargeDict={key: list() for key in typeList}                    # dict for charge
        typeIndexDict = {key:index+1 for index, key in enumerate(typeList)} # dict for typeIndex
        if self.verbose:
            print "Unique types = ", typeList
            print "typeIndexDict = ", typeIndexDict

        # Search for dataName and track with type
        for pid, ptclObj in self.particles.iteritems():

            ptclType   = ptclObj.type           # ptclType is key for dictionary
            ptclMass   = ptclObj.mass
            ptclCharge = ptclObj.charge

            dataList = typeMassDict[ptclType]   # Update pre-existing data list
            dataList.append(ptclMass)           # Append to running list
            typeMassDict[ptclType] = dataList   # Keep map of {ptclID: dataList}

            dataList = typeChargeDict[ptclType]  # Update pre-existing data list
            dataList.append(ptclCharge)          # Append to running list
            typeChargeDict[ptclType] = dataList  # Keep map of {ptclID: dataList}


        # Check for consistency mass
        for key, val in typeMassDict.iteritems():
            valSet = set(val)                  # Generate unique set of values
            if len(valSet) > 1:                # and check for length=1 (consistency)
               print "More than one mass found for ", key
               print "Check Particle() object initialization"
               sys.exit(3)
            typeMassDict[key] = val[0]         # If unique value make correct dict val

        # Check for consistency charge
        for key, val in typeChargeDict.iteritems():
            valSet = set(val)                  # Generate unique set of values
            if len(valSet) > 1:                # and check for length=1 (consistency)
               print "More than one charge found for ", key
               print "Check Particle() object initialization"
               sys.exit(3)
            typeChargeDict[key] = val[0]       # If unique value make correct dict val

        # Pack mass,charge,typeIndex to new dictionary
        typeInfoDict = dict()
        for key in typeIndexDict:
            index  = typeIndexDict[key]
            mass   = typeMassDict[key]
            charge = typeChargeDict[key]
            typeInfoDict[key] = [index, mass, charge]

        return typeInfoDict



    def scatterPlot(self):
        """
        Generate a scatter plot of particle positions
        """

        try:
            import matplotlib.pyplot as plt
        except:
            print "Error: matplotlib.pyplot module not found, check PYTHONPATH or install from source"
            sys.exit(3)

        try:
            from mpl_toolkits.mplot3d import Axes3D
        except:
            print "Error: mpl_toolkits.mplot3d module not found, check PYTHONPATH or install from source"
            sys.exit(3)

        # Unpack position data and put into list
        posDat = list()
        for pid, ptclObj in self.particles.iteritems():
            pos = ptclObj.position
            posDat.append(pos)

        # Analyze results
        col1 = [x[0] for x in posDat]
        col2 = [x[1] for x in posDat]
        col3 = [x[2] for x in posDat]

        # Make 3d scatter
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(col1, col2, col3)
        ax.set_xlabel("x [A]")
        ax.set_ylabel("y [A]")
        ax.set_zlabel("z [A]")
        
        # Name plot and output to file
        distName = "particleContainer.png"
        if self.verbose:
            print "Plotting ", distName
            
        plt.savefig(distName)
        plt.show()
        plt.close()

