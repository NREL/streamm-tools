"""
Class data structures for simulation data
"""
import copy, sys, os, shutil

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

    A template directory can be set. This is a common location for derived classes
    to access any template files needed to setup simulation

    NOTE:
    1. This is not the simulation code itself and does no calculations.
    """

    def __init__(self, name, verbose=False):
        """
        Constructor for a general simulation object.
        
        Args:
            name (str): String identifier for object (eg GaussianJuly21)
        """
        if verbose:
            print "Simulation base class constructor called"

        if isinstance(name, str):
            self.simulationName = name   # String name of simulation code (eg GaussianJuly21)
        else:
            print "1st arg should be string name for the simulation"
            raise TypeError

        # Debug/status flag
        self.verbose = verbose

        # Main attributes
        self.simulationExec = ""          # String name of simulation code executable (eg lmp)
        self.inputFileNames = list()      # List of file name strings (SWS: full paths?)
        self.simDir         = str()       # Location of simulation (where input files scripts etc should be copied)
        self.isStrucSet     = False       # Flag if object has structure container set
        self.topDir         = os.getcwd() # Current location of calling module
        self.templateDir    = "./"        # Location of directory containing any needed templates


    def __str__(self):
        """
        Print out name of simulation object (for logging/status)
        """
        return self.simulationName


    def printStruc(self):
        """
        Dump contents of held structure container
        """
        if self.isStrucSet:
            print "---------- Simulation object contains ------------"
            print self.strucC
        else:
            print "No structure container set"


    def __del__(self):
        """
        Destructor, clears object memory
        """

        if self.verbose:
            print "Cleaning simulation object"
            
        del self.inputFileNames
        if self.isStrucSet:
            del self.strucC


    def getSimName(self):
        """
        Return name used to create the simulation object
        """
        return self.simulationName


    def setTemplateDir(self, tdir):
        """
        Set template directory location
        """
        if not os.path.exists(tdir):
            print "Template directory does not exist... check full path \n"
            sys.exit(0)

        self.templateDir = tdir


    def setStructureContainer(self, strucC):
        """
        Setter for the structure container.
        """
        self.strucC = strucC
        self.isStrucSet = True


    def copyStructureContainerInto(self, strucC):
        """
        Set the structure container. Deep copy performed,
        so that external changes to structure container are not
        reflected here.
        """
        self.strucC = copy.deepcopy(strucC)        
        self.isStrucSet = True


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
        print "No Simulation:readOutput method defined for pure base class"
        sys.exit(0)


    def writeInput(self, fileName):
        """
        This is the 'effective' base class interface for a method
        that writes an input file based on the internal attributes of an instance
        of the Simulation object

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
        print "No Simulation:writeInput method defined for pure base class"
        sys.exit(0)


    def createSimulation(self):
        """
        Checks for existence of a top level simulation directory and writes out
        all files needed for running a simulation.

        Files copied/output are contained in the attribute 'inputFileNames'.

        In principle many input files/scripts could be copied to this location.
        If directory not found, then directory is created. Directory is creating
        from top level of where this class is executed
        """

        # Check for run directory
        if (not os.path.exists(self.simDir)):
            print self.simDir, "does not exist... creating"
            os.mkdir(self.simDir)


        # For all simulation files, move into run directory
        for inFile in self.inputFileNames:

            fromInFile = os.path.join(self.topDir, inFile)
            mvInFile   = os.path.join(self.topDir, self.simDir, inFile)
            shutil.move(fromInFile, mvInFile)

            if self.verbose:
                print "Moved input file to ", mvInFile
