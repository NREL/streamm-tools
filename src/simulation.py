"""
Class data structures for atomic data
"""
import copy, sys

class Simulation:
    """
    Data structure for describing how a simulation can be performed.
    The attributes stored can be used by a caller to construct/gather the
    appropriate components to run a simulation.

    There will be a tight coupling between objects of this class and the
    StructureContainer objects. Information from the StructureContainer objects
    will in general be used as part of the input to a simulation generated
    by using the data from this object

    NOTE:
       1. This is not the simulation code itself and does no calculations.
       2. Objects should be thought of as an instance of a simulation that
          converts a StructureContainer at step1 to a StructureContainer at step2
    """

    def __init__(self, pos=[], type="blank", charge=0.0, mass=1.0):
        """
        Constructor for a general simluation object.

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



class SimulationContainer:
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
        if self.verbose:
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
