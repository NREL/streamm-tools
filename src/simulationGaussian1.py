"""
Derived class for simulation object
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
    The Gaussian input file has two links 1.) one ground-state 2). TD DFT for excited state calc
    e.g.

    %chk=acc1_D1_eff_R3R300_A84___n1.chk 
    #P opt b3lyp/6-31g(d) geom=check guess=check  pop=full density=current

    acc1_D1_eff_R3R300_A84___n1 ground-state 

    0 1 
     <atoms>

    --Link1-- 
    %chk=acc1_D1_eff_R3R300_A84___n1.chk 
    #P b3lyp td=nstates=12/6-31g(d) geom=check guess=check  pop=full density=current
    """

    def __init__(self, name, verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """

        # Base class constructor is called
        sim.Simulation.__init__(self, name, verbose)

        # Chop off 'n1,n2..' from name
        self.oligomerNum = name.split('_')[-1]
        if ('n' not in self.oligomerNum):
            print "Run name does not contain n1, n2 or n*"
            sys.exit(0)

        if verbose:
            print "Simulation derived class 'Gaussian1' constructor called"

        self.verbose = verbose
        self.functional = str()  
        self.basisOPT = str()     # Basis for geo optim
        self.basisEXT = str()     # Basis for TD-DFT


    def getOligomerNum(self):
        """
        Return oligomer number string from run name. Assuming run name for
        contains n1,n2....
        """
        return self.oligomerNum


    def setCoeffs(self, fct, bsOPT, bsEXT=None):
        """
        Set the main model parameters of functional and basis

        Args:
            fct   (str): Gaussian functional
            bsOPT (str): Basis set to determine accuracy for geo opt.
            bsEXT (str): Basis set to determine accuracy for TD-DFT step
                         If not specified, then set to bsOPT
        """
        
        if not isinstance(fct, str):
            print "setCoeffs: fct should be a string"
            sys.exit(3)
        else:
            self.functional = fct


        if not isinstance(bsOPT, str):
            print "setCoeffs: bsOPT should be a string"
            sys.exit(3)
        else:
            self.basisOPT = bsOPT


        if bsEXT == None:
            self.basisEXT = bsOPT
        else:
            if not isinstance(bsEXT, str):
                print "setCoeffs: bsEXT should be a string if specified"
                sys.exit(3)
            else:
                self.basisEXT = bsEXT


    def __del__(self):
        """
        Destructor, clears object memory
        """

        # Base class destructor is called ?? needed
        sim.Simulation.__del__(self)

        if self.verbose:
            print "Cleaning derived simulation object Gaussian1"


    def _parseBasisString(self, bstr):
        """
        Parse functional/basis string from database
        bstr="trans//camb3lyp/6-31g"
        """
        blist = bstr.split('/')     # [ 'trans', '', 'camb3lyp', '6-31g']
        blist[:] = (value for value in blist if value != '') # Remove all spaces
        return blist


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

        # Hardwired run/link strings (associated with this derived class)
        paramRunStrRX = "#P opt " + self.functional + "/" + self.basisOPT
        paramRunStrRX += " geom=check guess=check  pop=full density=current"
        paramRunStrTD = "#P " + self.functional + " td=nstates=12/" + self.basisEXT
        paramRunStrTD += " geom=check guess=check  pop=full density=current"

        #
        # Open file, write header info
        #
        fileObj = open( inputName, 'w' )
        fileObj.write("%chk=" + jobName + ".chk \n")
        fileObj.write(paramRunStrRX + "\n")
        fileObj.write("\n")
        fileObj.write(jobName + " ground-state \n")
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

        #
        # Write out link for time-dependent calculation
        #
        fileObj.write('\n')
        fileObj.write('--Link1-- \n')
        fileObj.write("%chk=" + jobName + ".chk \n")
        fileObj.write(paramRunStrTD + "\n")
        fileObj.write("\n")        
        fileObj.write(jobName + " TD-DFT \n")
        fileObj.write("\n")
        fileObj.write("0 1 \n")

        # Close Gaussian file
        fileObj.close()
