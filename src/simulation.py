"""
Class data structures for simulation data
"""
import copy, sys

from structureContainer import StructureContainer

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

    def __init__(self, name):
        """
        Constructor for a general simulation object.
        
        Args:
            name (str): String identifier for object (eg GaussianJuly21)
        """

        if isinstance(name, str):
            self.simulationName = name   # String name of simulation code (eg GaussianJuly21)
        else:
            print "1st arg should be string name for the simulation"
            raise TypeError
            
        self.simulationExec = ""     # String name of simulation code executable (eg lmp)
        self.inputFileNames = list() # List of file name strings (SWS: full paths?)
        # others?
        # self.strucC = StructureContainer()

    def __str__(self):
        """
        Print out name of simulation object (for logging/status)
        """
        return self.simulationName


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.inputFileNames






class SimulationContainer:
    """
    Container class that associates Simulation objects with StructureContainer objects
    A Simulation is an object that is capable of 'constructing' an external simulation
    that uses a StructureContainer as an input and generates a StructureContainer as an
    output. Therefore a Simulation object should be thought of a 'function' that transforms
    a StructureContainer at state n (sC_n) to one at sC_(n+1).

    One of the primary data structures is a dictionary mapping StructureContainer-s to
    Simulation-s. An element of this dictionary will be tuple of StructureContainer-s
    mapped to a Simulation object. So:

      n      --> state index
      sC_n   --> instance of StructureContainer at state n
      simObj --> instance of Simulation

      with a dictionary element of the form  n: (sC_n, simObj, sC_n+1)
      eg.  { 1: (sC_1, simObjGaussian, sC_2), 2:(sC_2, simObjLammps, sC_3),....}
    """

    def __init__(self, verbose=False):
        """
        Constructor: sets up a dictionary for indexing 'Particle' objects

        Args:
            verbose (bool): flag for printing status/debug info
        """
        self.verbose=verbose     # verbose flag for debugging/status
        self.sc2sim = dict()     # Main object mapping structure containers to sims


    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        if self.verbose:
            print "Cleaning simulation container"

        del self.sc2sim


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
        return len(self.sc2sim)


    def __str__(self):
        """
        'Magic' method for printng contents of container
        """

        simStr="\n Contains simulations objects: \n"
        for simIndex, simTuple in self.sc2sim.iteritems():
            simObj = simTuple[1]
            simStr = simStr + " " + str(simIndex) + " " + str(simObj) + "\n"
        return simStr


    def __setitem__(self, simIndex, simTuple):
        """
        'Magic' method implementing obj[]=value operator
        """

        sC_n1  = simTuple[0]  # structureContainer at state n
        simObj = simTuple[1]  # Simulation object
        sC_n2  = simTuple[2]  # structureContainer at state n+1

        if not isinstance(sC_n1, StructureContainer):
            print "1st element of tuple should be a StructureContainer object"
            raise TypeError
            sys.exit(3)
            
        if not isinstance(simObj, Simulation):
            print "2nd element of tuple should be a Simulation object"
            raise TypeError
            sys.exit(3)
            
        if not isinstance(sC_n2, StructureContainer):
            print "3rd element of tuple should be a StructureContainer object"
            raise TypeError
            sys.exit(3)            

        if simIndex not in self.sc2sim.keys():
            e0 = copy.deepcopy(sC_n1)    
            e1 = copy.deepcopy(simObj)
            e2 = copy.deepcopy(sC_n2)
            self.sc2sim[simIndex]=(e0, e1, e2)
        else:
            print "Using [] operator for pre-existing key"
            sys.exit(3)


    def __getitem__(self, simIndex):
        """
        'Magic' method implementing obj[] operator.
        Operations on returned elements change container
        """
        return self.sc2sim[simIndex]


    def __delitem__(self, simIndex):
        """
        'Magic' method implementing del obj[] operator
        """
        del self.sc2sim[simIndex]


    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        return self.sc2sim.iteritems()


    def __contains__(self, gid):
        """
        'Magic' method implementing in keyword (key in obj')....
        """
        return gid in self.sc2sim
