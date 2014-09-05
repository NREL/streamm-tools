"""
Derived class for structure container
"""

try:
    import simulation as sim
except:
    print "Error: structureContainer module not found"
    print "Check setup.sh path script"
    sys.exit(3)


class SimulationGaussian1(sim.Simulation):
    """
    Dervied class implementing input/output methods specific to an application
    """

    def __init__(self, name, verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """

        # Base class constructor is called
        sim.Simulation.__init__(self, name, verbose)

        if verbose:
            print "Simulation derived class 'Gaussian1' constructor called"

        self.verbose = verbose
        self.paramRunStr = str()
        self.labelRunStr = str()


    def setCoeffs(self, pStr, lStr):
        """
        Set the options/parameter string for Gaussian run
        
        Args:
            pStr (str): Parameter string at top of input file
                        eg #P opt=Cartesian iop(1/7=12000) AM1 scf=qc
            lStr (str): Label string for the run eg  'semi-empirical pre-optimiz'
        """
        
        if not isinstance(pStr, str):
            print "setCoeffs: paramRunStr should be a string"
            sys.exit(3)

        if not isinstance(lStr, str):
            print "setCoeffs: labelRunStr should be a string"
            sys.exit(3)

        self.paramRunStr = pStr
        self.labelRunStr = lStr


    def __del__(self):
        """
        Destructor, clears object memory
        """

        # Base class destructor is called ?? needed
        sim.Simulation.__del__(self)

        if self.verbose:
            print "Cleaning derived simulation object Gaussian1"


    """
    def writeInput(self,....):
        writeInputSetting()
        writeInputDataFile()
        writeInputOther()

    def writeInputDataFile(self, inputName):        
    """


    def writeInput(self, inputName):
        """
        Write out a Gaussian input data file from all available held
        data (particles, bonds-zmatrix)
        Assumptions about one 'link' step made and naming convention of referenced 'chk' files
        and input file name

        Args:
            inputName    (str)  name of Gaussian input file (.com) to write out
        """

        jobName = inputName.split('.com')[0]  # Assumes input filename struc is ***..com
        n_atoms = len(self.strucC.ptclC)      # Obtaining particle size from container

        #
        # Open file, write header info
        #
        fileObj = open( inputName, 'w' )
        fileObj.write("%chk=" + jobName + ".chk \n")
        fileObj.write(self.paramRunStr + "\n")
        fileObj.write("\n")
        fileObj.write(jobName + " " + self.labelRunStr + "\n")
        fileObj.write("\n")
        fileObj.write("0 1 \n")

        #
        # Write out atoms section
        #
        ptclFormatStr = "%s %9.4f %9.4f %9.4f \n"
        for pid, ptclObj in self.strucC.ptclC:
            pos = ptclObj.position
            typ = ptclObj.type
            fileObj.write( ptclFormatStr % (typ, pos[0], pos[1], pos[2] ) )
        fileObj.write('\n')
        
        # Close Gaussian file
        fileObj.close()
