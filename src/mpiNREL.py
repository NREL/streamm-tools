"""
Wrapper classes for MPI python using the Peregrine module enviroment

The derived classes implement the interface defined in the base class
'ParallelMsgr'. These derived classes currently implement communication using:
  * boost.mpi
  * mpi4py
  * serial

NOTE: continue to end of this documentation to see short tutorial methods on use

"""
import os, sys, math, random, time


class ParallelMsgr:
    """
    Base class defining the interface to parallel communication methods
    """

    def __init__(self, rk, sz, verbose=False):
        """
        Constructor

        Args:
            rk: rank of processor
            sz: size of communicator
            verbose: flag for printing output info
        """

        # Flag for debug printing
        self.verbose = verbose

        # Basic comm variables
        self.rank     = rk
        self.commSize = sz

        if (self.verbose and self.rank == 0):
            print " "
            print "Base class ParallelMsgr constructor called"


    def __del__(self):
        """
        Destructor
        """

        if (self.verbose and self.rank == 0):
            print "ParallelMsgr destructor called"
            print " "


    def getRank(self):
        """
        ID of current processor

        Returns:
            integerr rank of current processor
        """

        return self.rank


    def getCommSize(self):
        """
        Number of processors in world communicator
        """

        return self.commSize


    def barrier(self):
        """
        Must be implemented in derived 
        """

        print "Not implemented in base class"
        sys.exit(3)

    def gatherList(self, data, root=0):
        """
        Must be implemented in derived
        """

        print "gatherList not implemented in base class"
        sys.exit(3)


    def reduceSum(self, data, root=0):
        """
        Must be implemented in derived
        """

        print "reduceSum not implemented in base class"
        sys.exit(3)


    def allReduceSum(self, data):
        """
        Must be implemented in derived
        """

        print "allReduceSum not implemented in base class"
        sys.exit(3)


    def tupleOfLists2List(self, tup):
        """
        Take a tuple of lists and convert to single list
        eg ([a,b],[c,d],[e,f],...) -->  [a,b,c,d,e,f]
        """

        bigList = []
        try:
            for i in tup:
                bigList = bigList + i
        except:
            pass

        return bigList


    def splitListOnProcs(self, data):
        """
        Split up global input list equally (as possible) and return the
        'chunk' of list belonging to each processor. Note, no communication
        is performed, input data is known globally.

        Args:
            data: list of data [..]
        Returns:
            list of data 'for' this processor
        """

        parts     = self.commSize                  # Number of procesors
        numData   = float(len(data))               # Global size of data list
        dataPProc = int(math.ceil(numData/parts))  # Est. num of data on each proc

        # Make 'parts' number of chunks
        plist=[data[x:x+dataPProc] for x in xrange(0, len(data), dataPProc)]

        # Error check or return results
        if len(plist) != parts:
            if self.rank == 0:
                print " " 
                print "Partitioning failed, check data length and #-procs"
                print "   len(plist) = ", len(plist)
                print "        parts = ", parts
            sys.exit(0)
        else:
            return plist[self.rank]


class MsgrBoost(ParallelMsgr):
    """
    Derived class implementing interface in the base class using
    Boost MPI libraries
    """

    def __init__(self, commObj, verbose=False):
        """
        Constructor
        """

        # Use comm object to initialize
        self.commObj = commObj
        self.rank = commObj.rank
        self.size = commObj.size
        self.world = commObj.world

        # Ensuring base class called
        ParallelMsgr.__init__(self, self.rank, self.size, verbose)

        # Flag for debug printing
        self.verbose = verbose
        if (self.verbose and self.rank == 0):
            print "Derived class MsgrBoost constructor called"

        if (self.verbose):
            for proc in range(self.size):
                if proc == self.rank:
                    print "I am process %d of %d." % (self.rank, self.size)
                self.barrier()
            self.barrier()


    def bcast(self, data, root=0):
        """
        Broadcast from root proc to world comm
        """

        gdata = self.commObj.broadcast(self.world, data, root)
        self.barrier()
        return gdata


    def barrier(self):
        """
        Barrier for world comm
        """
        self.world.barrier()


    def gatherList(self, data, root=0):
        """
        Gather data in list from world comm on root proc
        and return as a continuous list
        """

        tupOfLists = self.commObj.gather(self.world, data, root)
        oneList= ParallelMsgr.tupleOfLists2List(self, tupOfLists)
        self.barrier()
        return oneList


    def reduceSum(self, data, root=0):
        """
        Add data elements in localData across processors and put
        result on root proc
        """

        gdata = self.commObj.all_reduce(self.world, data, lambda a,b: a+b, root)
        self.barrier()
        return gdata


    def allReduceSum(self, data):
        """
        Add data elements in localData across processors and put on
        all world procs
        """
        gdata = self.commObj.all_reduce(self.world, data, lambda a,b: a+b)
        self.barrier()
        return gdata



class MsgrMpi4py(ParallelMsgr):
    """
    Derived class implementing interface in the base class using
    Mpi4py libraries
    """

    def __init__(self, commObj, verbose=False):
        """
        Constructor
        """

        # Use comm object to initialize
        self.commObj = commObj
        self.comm = commObj.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        # Ensuring base class called
        ParallelMsgr.__init__(self, self.rank, self.size, verbose)

        # Flag for debug printing
        self.verbose = verbose
        if (self.verbose and self.rank == 0):
            print "Derived class MsgrMPI4PY constructor called"
            
        if (self.verbose):
            for proc in range(self.size):
                if proc == self.rank:
                    print "I am process %d of %d." % (self.rank, self.size)
                self.barrier()
            self.barrier()


    def bcast(self, localData, root=0):
        """
        Broadcast to world communicator
        """

        gdata = self.comm.bcast(localData,root)
        return gdata


    def barrier(self):
        """
        Barrier to world comm
        """
        self.comm.barrier()


    def gatherList(self, localData, root=0):
        """
        Gather data in list from world comm on root proc
        and return as a continuous list
        """
        listOfLists = self.comm.gather(localData, root)
        oneList = ParallelMsgr.tupleOfLists2List(self, listOfLists)
        self.barrier()
        return oneList

    def reduceSum(self, localData, root=0):
        """
        Add data elements in localData across processors and put
        result on root proc
        """
        gdata = self.comm.reduce(localData, self.commObj.SUM, root)
        self.barrier()
        return gdata


    def allReduceSum(self, localData):
        """
        Add data elements in localData across processors and put on
        all world procs
        """
        MPI = self.commObj
        gdata = self.comm.allreduce(localData, MPI.SUM)
        self.barrier()
        return gdata



class MsgrSerial(ParallelMsgr):
    """
    Derived class implementing interface in the base class using
    no MPI (serial) libraries
    """

    def __init__(self, verbose=False):
        """
        Constructor
        """
        self.rank = 0
        self.size = 1

        # Ensuring base class called
        ParallelMsgr.__init__(self, self.rank, self.size, verbose)

        # Flag for debug printing
        self.verbose = verbose
        if (self.verbose and self.rank == 0):
            print "Derived class MsgrSerial constructor called"


    def bcast(self, data, root=0):
        """
        Broadcast from root proc to world comm
        """
        return data


    def barrier(self):
        """
        Barrier for world comm
        """
        return


    def gatherList(self, data, root=0):
        """
        Gather data in list from world comm on root proc
        and return as a continuous list
        """
        return data


    def reduceSum(self, data, root=0):
        """
        Add data elements in localData across processors and put
        result on root proc
        """
        return data


    def allReduceSum(self, data):
        """
        Add data elements in localData across processors and put on
        all world procs
        """
        return data




def getMPIObject(verbose=True):
    """
    Driver for comm classes.
    Selects MPI-comm module if found and builds appropriate derived class

    Returns: an mpi object of the appropriate derived class
    """
    mpiObj=None
    try:
        from mpi4py import MPI as mpi
        mpiObj = MsgrMpi4py(mpi, verbose)
        if mpiObj.getRank() == 0:
            print "Found mpi4py module \n"
        mpiObj.barrier()
    except:
        pass

    if (mpiObj == None):
        try:
            import boost.mpi as mpi
            mpiObj = MsgrBoost(mpi, verbose)
            if mpiObj.getRank() == 0:
                print "Found boost.mpi module \n"
            mpiObj.barrier()
        except:
            pass

    if (mpiObj == None):
        mpiObj = MsgrSerial(verbose)
        print "No mpi module found, will use serial. \n"

    return mpiObj


def showExampleUsage():
    """
    Instructions and a sample code for testing the MPI routines.

    Usage is as follows:
       1. paste code below to a file --> e.g. 'commTst.py'
       2. While logged onto Peregrine load an MPI python module e.g.
          'module load boost/1.54.0.a/impi-intel'
       3. mpirun -n 2 python commTst.py

    --------------------------------------------------------------

    # commTst.py
    import mpiNREL
    p = mpiNREL.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    x = 0              # Sets x on all processors
    p.barrier()
    if rank == 0:      # Resets x on proc 0
        x = 12.3456
    p.barrier()

    # Check x setting on each processor
    for proc in range(size):
        if proc == rank:
            print "Before broadcast x = ", x, " on processor ", proc
        p.barrier()
    p.barrier()

    # Broadcast x to all processors (default root proc is '0')
    xAll = p.bcast(x)
    p.barrier()

    # Check x setting on each processor
    for proc in range(size):
        if proc == rank:
            print "After broadcast x = ", xAll, " on processor ", proc
        p.barrier()
    p.barrier()
    """

    print showExampleCode.__doc__

