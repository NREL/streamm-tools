"""
Derived class for simulation object
"""

import sys, os

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

    NOTE: lc-wpbe removes 'density=current' from input lines
    """

    def __init__(self, name, verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly

        Assuming that object is created with a name related to a tag name from
        donoracceptorsystems eg.  acc1_D1_R2R200_Sp2_00_A48__Sp2_00__n2
        """

        # Base class constructor is called
        sim.Simulation.__init__(self, name, verbose)

        # Local list of split name on '_'
        # eg ['acc1', 'D1', 'R2R200', 'Sp2', '00', 'A48', '', 'Sp2', '00', '', 'n2']
        sp = name.split('_')

        # Simple check for correct tag format (first part must have acc in string)
        if ('acc' not in sp[0]):
            print "Expecting object to be created with tag name containing 'acc'"
            sys.exit(0)

        # Chop off 'n1,n2..' from name
        self.oligomerNum = sp[-1]
        if ('n' not in self.oligomerNum):
            print "Run name does not contain n1, n2 or n*"
            sys.exit(0)

        # Convert oligomer string value to integer (eg. 'n1' --> 1)
        self.oligomerInt = int(self.oligomerNum.strip('n'))


        # Find 'accuracy tag' and chop off of simulation object name
        accTag  = sp[0] + "_"                                     # eg 'acc1_'
        tmpStr = name.replace(accTag,"")                          # eg 'D1_R2R200_Sp2_00_A48__Sp2_00__n2'
        self.bareTag = tmpStr.replace("_" + self.oligomerNum, "") # eg 'D1_R2R200_Sp2_00_A48__Sp2_00_'
        self.simDir = self.bareTag                                # Set (public) sim dir in base to bareTag

        self.verbose = verbose 
        if verbose:
            print "Simulation derived class 'Gaussian1' constructor called"
            print "Setting base class simDir to bareTag = ", self.bareTag

        self.repoLocation    = self._generateDefaultRepoLocation()

        self.functional      = str()  # Functional string for Gaussian
        self.basisOPT        = str()  # Basis for geo optim
        self.basisEXT        = str()  # Basis for TD-DFT
        self.dbBasisStr      = str()  # Database func/basis 'code'
        self.metaDataDict    = dict() # Metadata dictionary (json)
        self.inputNamePrefix = str()  # Name using in sim/input files (wo the .com, .json etc)
        self.nstates         = 12     # Default to number of states included in TD-DFT calculation

        # Name of PBS template file (assumed to be in location set by base class)
        self.pbsTemplateFile = "donoracceptor.pbs.template"

        # This is a special accuracy setting needed for diffuse basis functions
        # self.accuracyStr = ""
        # self.accuracyStr = " " + "Integral=(Acc2E=11)"
        self.accuracyStr = " " + "Integral=(Acc2E=12)"

        self.groundCalc  = " " + "density=current"   # Default setting calculates ground  state density
        self.excitedCalc = " " + "density=current"   # Default setting calculates excited state density


    def getBareTag(self):
        """
        Return 'bare tag' without accuracy label
        eg D1_R2R200_Sp2_00_A48__Sp2_00__n2
        """
        sim.Simulation.__init__(self, name, verbose)
        return self.bareTag


    def getDBBasisStr(self):
        """
        Return dbBasisStr
        """
        return self.dbBasisStr


    def setRepoLocation(self,rL):
        """
        Set repo location path so PBS file can be
        constructed from a template
        
        Args:
            rL (str) full path to repo location
        """

        if isinstance(rL, str):
            self.repoLocation = rL
        else:
            print "Default results repo location set --> ", self.repoLocation


    def setDBBasisStr(self,s):
        """
        Set dbBasisStr
        
        Args:
            s (str) basis string in database
        """
        self.dbBasisStr=s


    def getOligomerNum(self):
        """
        Return oligomer number string from run name. Assuming run name for
        contains n1,n2....
        """
        return self.oligomerNum


    def setMetaData(self, mdDict):
        """
        Set dictionary containing metadata (from json file info)
        """
        self.metaDataDict = mdDict


    def setCoeffs(self, fct, bsOPT, bsEXT=None):
        """
        Set the main model parameters of functional and basis
        
        Args:
            fct   (str): Gaussian functional
            bsOPT (str): Basis set to determine accuracy for geo opt.
            bsEXT (str): Basis set to determine accuracy for TD-DFT step. If not specified, then set to bsOPT
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

        # Special settings for certain functionals
        if fct == 'lc-wpbe':
            self.excitedCalc = " "
            print "Special setting for lc-wpbe... removing density=current for excited state calc"


    def __del__(self):
        """
        Destructor, clears object memory
        """

        # Base class destructor is called ?? needed
        sim.Simulation.__del__(self)

        if self.verbose:
            print "Cleaning derived simulation object Gaussian1"


    def writeInput(self):
        """
        Write out a Gaussian input data file from all available held data (particles, bonds-zmatrix)
        Assumptions about one 'link' step made and naming convention of referenced 'chk' files
        and input file name

        Args:
            inputName (str):  name of Gaussian input file (.com) to write out
                              Assumes that structure of input name is acc.functional.basis_D1_R300....
        """

        #
        # This input name prefix will be used in all private
        # 'self._writeXXX...' methods  called here
        #
        tmpName = "acc_" + str(self.functional) + "_" + str(self.basisOPT) + "_" + \
                  self.bareTag + "_" + self.oligomerNum

        # Translate problematic characters in input file names
        self.inputNamePrefix = self._convertInputName(tmpName)

        # Write input files
        self._writeGaussianComFile()    # Write Gaussian file
        self._writeGaussianComFile_r1() # Write Gaussian file restart RX stage
        self._writeGaussianComFile_r2() # Write Gaussian file restart TD stage
        self._writeXYZFile()            # Write plain config file
        self._writeMetaDataFiles()      # Write auxillary files
        self._writePBSFile()            # Write PBS script file from template



    #########################################################################
    #                                                                       #
    #                         Private class methods                         #
    #                                                                       #
    #########################################################################

    def _parseBasisString(self, bstr):
        """
        Parse functional/basis string from database
        bstr="trans//camb3lyp/6-31g"
        """
        blist = bstr.split('/')     # [ 'trans', '', 'camb3lyp', '6-31g']
        blist[:] = (value for value in blist if value != '') # Remove all spaces
        return blist


    def _convertInputName(self, fname):
        """
        Takes .com Gaussian input name and convert problematic strings for file systems
        The replacements are:
          , --> " "
          + --> "P"
          ( --> "L"
          ) --> "R"
        """
        tmpStr1 = fname.replace(",","").replace("+","P").replace("(","L").replace(")","R")
        return tmpStr1


    def _writeMetaDataFiles(self):
        """
        Using the metadata dictionary info, write out the .meta file and .json files
        Edits the metadata info with any current info on basis, accuracy, date and number of states
        in the TD-DFT steps
        
        eg {'metadata': {'date_time': '2014-10-08T14:21:39.641616', 'acceptors': ['A48'], 'donor_substituents': ['f0', 'R2'],
                         'donors': ['D1'], 'acceptor_substituents': [], 'basis': 'b3lyp', 'terminal_substituents': [],
                         'number': 1, 'spacers': ['Sp2'], 'terminals': [], 'tag': 'D1_R2R200_Sp2_00_A48__Sp2_00_',
                         'spacer_substituents': ['f0'], 'n': 1, 'nstates': 2, 'accuracy': 'XXaccuracyXX'}}
        """

        import json, datetime, time

        # Get timestamp
        dt = datetime.datetime.fromtimestamp(time.time())

        # Set json and meta file names
        jsonName = self.inputNamePrefix + ".json"
        metaName = self.inputNamePrefix + ".meta"

        # Update the metadata dictionary
        subDict = self.metaDataDict['metadata']  # Drill down to sub dictionary

        subDict['date_time'] = dt.isoformat()
        subDict['nstates']   = self.nstates
        subDict['number']    = self.oligomerInt
        subDict['n']         = self.oligomerInt
        subDict['basis']     = self.functional + "/" + self.basisEXT
        subDict['accuracy']  = self.functional + "/" + self.basisEXT

        # Write the json file
        with open(jsonName, 'w') as outfile:
            json.dump(self.metaDataDict, outfile, indent=2)

        # Write the meta file
        metaDataOut = str(self.metaDataDict['metadata'])
        with open(metaName, 'w') as outfile:
            outfile.write("metadata = " + metaDataOut)

        # Assign sim files to base class list
        self.inputFileNames.append(jsonName)
        self.inputFileNames.append(metaName)


    def _writeGaussianComFile(self):

        # Set base class list of input files
        inputName = self.inputNamePrefix + ".com"
        self.inputFileNames.append(inputName)
        jobName = self.inputNamePrefix        # Assumes input filename struc is ***..com
 
        # Hardwired run/link strings (associated with this derived class)
        paramRunStrRX = "#P opt " + self.functional + "/" + self.basisOPT
        paramRunStrRX += " pop=full" + self.groundCalc + self.accuracyStr
        paramRunStrTD = "#P " + self.functional + " td=nstates=" + str(self.nstates) + "/" + self.basisEXT
        paramRunStrTD += " geom=check guess=check pop=full" + self.excitedCalc + self.accuracyStr

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


    def _writeGaussianComFile_r1(self):

        # Set base class list of input files
        inputName = self.inputNamePrefix + ".r1.com"
        self.inputFileNames.append(inputName)
        jobName = self.inputNamePrefix        # Assumes input filename struc is ***..com
 
        # Hardwired run/link strings (associated with this derived class)
        paramRunStrRX = "#P opt " + self.functional + "/" + self.basisOPT
        paramRunStrRX += " geom=check guess=check pop=full" + self.groundCalc + self.accuracyStr
        paramRunStrTD = "#P " + self.functional + " td=nstates=" + str(self.nstates) + "/" + self.basisEXT
        paramRunStrTD += " geom=check guess=check pop=full" + self.excitedCalc + self.accuracyStr

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


    def _writeGaussianComFile_r2(self):

        # Set base class list of input files
        inputName = self.inputNamePrefix + ".r2.com"
        self.inputFileNames.append(inputName)
        jobName = self.inputNamePrefix        # Assumes input filename struc is ***..com

        # Hardwired run/link strings (associated with this derived class)
        paramRunStrTD = "#P " + self.functional + " td=nstates=" + str(self.nstates) + "/" + self.basisEXT
        paramRunStrTD += " geom=check guess=check pop=full" + self.excitedCalc + self.accuracyStr

        #
        # Open file, write header info
        # write out link for time-dependent calculation
        #
        fileObj = open( inputName, 'w' )
        fileObj.write("%chk=" + jobName + ".chk \n")
        fileObj.write(paramRunStrTD + "\n")
        fileObj.write("\n")        
        fileObj.write(jobName + " TD-DFT \n")
        fileObj.write("\n")
        fileObj.write("0 1 \n")

        # Close Gaussian file
        fileObj.close()


    def _writeXYZFile(self):
        """
        Write out the XYZ files with only the coordinates
        """

        # Set base class list of input files
        inputName = self.inputNamePrefix + ".xyz"
        self.inputFileNames.append(inputName)

        n_atoms = len(self.strucC.ptclC)      # Obtaining particle size from container

        #
        # Open file, write header info
        #
        fileObj = open( inputName, 'w' )
        fileObj.write(str(n_atoms) + "\n")
        fileObj.write("A Structure \n")

        #
        # Write out atoms section
        #
        ptclFormatStr = "%s %9.4f %9.4f %9.4f \n"
        for pid, ptclObj in self.strucC.ptclC:
            pos = ptclObj.position
            typ = ptclObj.type
            fileObj.write( ptclFormatStr % (typ, pos[0], pos[1], pos[2] ) )
        fileObj.write('\n')

        # Close xyz file
        fileObj.close()


    def _writePBSFile(self):
        """
        Use PBS template to generate PBS run script for simulation
        """

        # Set base class list of input files
        pbsName = self.inputNamePrefix + ".pbs"
        self.inputFileNames.append(pbsName)

        from string import replace

        # Open template file (using base class template location)
        pbsTemplateName = os.path.join(self.templateDir, self.pbsTemplateFile)
        f = file(pbsTemplateName)
        templ = f.read()
        f.close()

        templ = replace(templ, "<job_name>",      self.inputNamePrefix)
        templ = replace(templ, "<repo_location>", self.repoLocation)

        # Write substituted PBS file
        f = file(pbsName, "w")
        f.write(templ)
        f.close()


    def _generateDefaultRepoLocation(self):
        """
        Generate default repo location from a prefix (defined here) and the user environment
        name and a timestamp (eg /projects/opv/PUBLISH-'user'-datestamp)
        """

        from time import strftime, gmtime
        userName=os.getenv('USER')
        timeStamp = strftime("%Y_%b_%d_%a_%H_%M", gmtime())
        repoLocation="/projects/opv/PUBLISH-" + userName + "-" + timeStamp

        if self.verbose:
            print "Repo location = ", repoLocation

        return repoLocation
