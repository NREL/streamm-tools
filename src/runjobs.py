#!/usr/bin/env python

"""
A collection of classes for useful functions. Needs to be split
"""

######################################################################
# Load modules
######################################################################
try:
    import os, sys
except:
    print "Error: os and or sys not found"
    print "Check if python is on this system"
    sys.exit(3)

try:
    import string
except:
    print "Error: string not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import math
except:
    print "Error: math not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import re
except:
    print "Error: re not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    from time import sleep
except:
    print "Error: time module not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    from glob import glob
except:
    print "Error: time module not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    from optparse import OptionParser
except:
    print "Error: optparse not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import platform
except:
    print "Error: platform module not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import socket
except:
    print "Error: socket module not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import subprocess, shlex
except:
    print "Error: subprocess,shlex module not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import shutil, copy
except:
    print "Error: shutil and/or copy module(s) not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import fnmatch
except:
    print "Error: fnmatch module not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)
######################################################################


class Misc:
    """
    Collection of useful methods that do not clearly fit into other classes
    """ 

    def __init__(self, debugFlag=False):
        """
        Constructor
        """

        self.debugFlag = debugFlag
        if (self.debugFlag):
            print " "
            print "Misc class created"
            print " debugFlag set to --> ", debugFlag
            print " "


    def __del__(self):
        """
        Destructor
        """

        if (self.debugFlag):
            print " "
            print "Misc class closed"


    def recursive_glob(self, treeroot, pattern):
        """
        Recursive glob using the os.walk command to traverse
        directory structure. Starting from 'treeroot' directory,
        return list of full path names for files containing 'pattern'
        
        Args:
            treeroot (str): Directory name to start top of search
            pattern  (str): String pattern (can use * wildcard) for search
        Returns:
            List of path name string
        """

        results = []
        for base, dirs, files in os.walk(treeroot):
            goodfiles = fnmatch.filter(files, pattern)
            results.extend(os.path.join(base, f) for f in goodfiles)
        return results



    def getList(self, x0, dx, xmax):
        """
        Generate list from x0, x0+dx,..... xmax
        """

        xList=[]    
        xList.append(x0)
        xnext = x0+dx
        while xnext<xmax:
            xList.append(xnext)
            xnext = xnext+dx
        return xList



    def getTagForSorting(self, num):
        """
        Return a string with leading 0's so file systems etc
        sort is a normal way
        """

        if (num > 9999):
            print "getTagForSorting: tag failed. Job num > 9999"
            sys.exit(1)
        tagTemplate = 100000
        tagstr = str(tagTemplate+num)
        return tagstr.lstrip('1')



    def partitionList(self, data, partLen):
        """
        Partition list into chunks w the last chunk
        is the remainder of elements:
        
        data    --> eg [1.2, 3.3, 55.6, 0.34.......3.1,38.2] (odd)
        partLen --> eg 2
        chunks = [ [1.2, 3.3], [55.6, 0.34].... [38.2] ]
        """

        # Make 'parts' number of chunks
        chunks=[data[x:x+partLen] for x in xrange(0, len(data), partLen)]
        return chunks



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



class IO:
    """
    Class for managing output to files/screen
    """

    def __init__(self, fileName, debugFlag=False):
        """
        Constructor
        
        Args:
            fileName (str): string for file name
        """

        self.fileName = fileName
        self.debugFlag = debugFlag
        self.fileObj = open(fileName, 'w')
        if (self.debugFlag):
            print " "
            print "IO file object ", fileName, " opened"
            print " debugFlag set to --> ", debugFlag
            print " Flags can be reset in scripts by developers"
            print " "


    def __del__(self):
        """
        Destructor (closes file)
        """

        self.fileObj.close()
        if (self.debugFlag):
            print " "
            print "IO file object ", self.fileName, " closed"
            print " "

    def flush(self):
        """
        Flush the print buffer to the file object
        """
        self.fileObj.flush()


    def prtdb(self, prtString):
        """
        Print debug info if set
        """

        if (self.debugFlag):
            print prtString

    def prtScreen(self, prtString):
        print prtString

    def prtFile(self, prtString):
        """
        Only print to file
        """
        self.fileObj.write(prtString)
        self.fileObj.write(" \n")


    def prt(self, prtString):
        """
        Write to screen and file
        """
        print prtString

        self.fileObj.write(prtString)
        self.fileObj.write(" \n")


    def floatToStr(self, x, num=5):
        """
        Convert/format a float value to a string
        Args:
            x (float): value
            num (int): Size of format field
        """

        formatStr = "%." + str(num) + "f"
        xStr = formatStr % x
        return xStr


#class FileEditorWithSomething:
#    def __init__(self):
#        pass


class Job:
    """
    Job class representing data/members needed
    to represent an instance of a code that will run
    with a single set of input parameters
    """

    def __init__(self, repopath, exe,
                 inputfile, script, verbose=False):
        """
        Constructor
        """

        self.repopath = repopath     # Path to repo for common files
        self.exe = exe               # Name of executable binary
        self.inputfile = inputfile   # Name of input file
        self.script = script         # Name of script for PBS
        self.suffix = script.split('.')[-1]  # Suffix on PBS script
        self.verbose = verbose       # Flag for debug printing
        self.inputEditorSet  = False # Flag if input  file editor set
        self.scriptEditorSet = False # Flag if script file editor set

        if (self.verbose):
            print " "
            print "Job object template created with verbosity on"
        else:
            print " "
            print "Job object template created with verbosity off"

        # Check for repo path
        if (not os.path.exists(repopath)):
            print " "
            print "Repo path ", repopath, " not found"
            print "Edit createJobsDirs.py:repoPath by hand for your path"
            print " "
            sys.exit(1)

        # Check for full paths of needed files
        exePath    = os.path.join(self.repopath, self.exe)
        inputPath  = os.path.join(self.repopath, self.inputfile)
        scriptPath = os.path.join(self.repopath, self.script)
        if (not os.path.exists(exePath)):
            print "Executable ", exePath, " not found"
            sys.exit(1)
        if (not os.path.exists(inputPath)):
            print "Input file template ", inputPath, " not found"
            sys.exit(1)
        if (not os.path.exists(scriptPath)):
            print "Queue script ", scriptPath, " not found"
            sys.exit(1)
            
    # Destructor
    def __del__(self):
        if (self.verbose):
            print "Destructor called"

    # Setters (can be used to reset a member for a specific run)
    def setExecutable(self, exe):
        self.exe = exe
    def setRepoPath(self, repopath):
        self.repopath = repopath
    def setInputFile(self, inputfile):
        self.inputfile = inputfile
    def setRunDir(self, rundir):
        self.rundir = rundir
        self.rundirName = os.path.basename(rundir)
    def setScript(self, script):
        self.script = script
        self.suffix = script.split('.')[-1]  # script.mb --> mbn1

    def setInputEditor(self, editor):
        """
        Expects an editor object that implements the
        appropriate 'editFile(..)' method
        """

        self.inputEditorSet = True
        self.inputEditor = editor


    def setScriptEditor(self, editor):
        """
        Expects an editor object that implements the
        appropriate 'editFile(..)' method
        """

        self.scriptEditorSet = True
        self.scriptEditor = editor

        
    def info(self):
        """
        Dump state of object
        """

        print " "
        print "Data members set to: "
        print "   repopath  = ", self.repopath
        print "   rundir    = ", self.rundir
        print "   exe       = ", self.exe
        print "   inputfile = ", self.inputfile
        print "   script    = ", self.script
        print " "


    def createJobRun(self):
        """
        Create run directory and populate with appropriate files
        """

        # Local full paths of needed files
        exePath    = os.path.join(self.repopath, self.exe)
        inputPath  = os.path.join(self.repopath, self.inputfile)
        scriptPath = os.path.join(self.repopath, self.script)

        # If run directory doesnt exist, create it
        if (os.path.exists(self.rundir)):
            print self.rundir, " already exists, skipping to next \n"
        else:
            os.makedirs(self.rundir)
        
        if (self.verbose):
            print " "
            print "Creating/setting up job directory ", self.rundir
            print "   exePath    = ", exePath
            print "   inputPath  = ", inputPath
            print "   scriptPath = ", scriptPath
            print " "

        # Copy repo files to run directory
        shutil.copy2(exePath,    self.rundir)
        shutil.copy2(inputPath,  self.rundir)
        shutil.copy2(scriptPath, self.rundir)

        # Copy repo script name to run directory with run dir name
        # Note: this CHANGES class data
        oldScriptFullPath = os.path.join(self.rundir, self.script)
        self.script = self.rundirName + "." + self.suffix
        newScriptFullPath = os.path.join(self.rundir, self.script)
        shutil.move(oldScriptFullPath, newScriptFullPath)


    def changeInputFile(self):
        """
        Use imported edit methods to change input file for run
        The imported method must implement an 'editfile(file)' method
        """

        if (self.inputEditorSet):
            currFile=os.path.join(self.rundir, self.inputfile)
            self.inputEditor.editFile(currFile)
        else:
            print "Input file editor not set. ", self.inputfile, " is unchanged"
            sys.exit(1)


    def changeScriptFile(self):
        """
        Use imported edit methods to change script file for run
        The imported method must implement an 'editfile(file)' method
        """

        if (self.scriptEditorSet):
            currFile=os.path.join(self.rundir, self.script)
            self.scriptEditor.editFile(currFile)
        else:
            print "Scriptfile editor not set. ", self.script, " is unchanged"
            sys.exit(1)



class JobStatus:
    """
    Class containing members needed to inquire about job status once
    already created and populated by necessary files/scripts
    """

    def __init__(self, verbose=False):
        """
        Constructor
        """
        self.verbose = verbose       # Flag for debug printing


    def __del__(self):
        """
        Destructor
        """
        if (self.verbose):
            print "Destructor called"
    

    def isJobInQueue(self, jobdir, jobscript):
        """
        This method will be inherited by derived classes as all
        jobs are assumed to be in queueing system and running jobs
        will have same signature
        
        NOTE: method logic assumes that .pbs file name is the job
        name to look for
        """

        # jobdir=/home/ssides/opv/opv_generator/mols/D1_R16R100_A1_R3R30000_
        # jobscript = acc1_D1_R16R100_A1_R3R30000__n1.pbs 
        # jobdirName = os.path.basename(jobscript) # D1_R16R100_A1_R3R30000_ 
        jobname = jobscript.split('.')[0]  # acc1_D1_R16R100_A1_R3R30000__n1

        #
        # Set submit command
        # Takes qstat -f output and returns 8 lines of output after the 'jobname' is found
        # Must grep this because location of job_state line changes
        #
        qcommand = " qstat -f | awk '/" + jobname + \
                   "/ {for(i=1; i<=8; i++) {getline; print}}' |grep job_state"

        # Runs submission command and grabs output
        proc = subprocess.Popen(qcommand, shell=True,
                                stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE)

        # Success code from qcommand
        return_code = proc.wait()

        # If command fails, then nothing found in queue
        if ( return_code != 0):
            jobInQueue = False
        else:
            # Read from pipes, line should be eg 'job_state = R'
            for line in proc.stdout:
                splitStr = line.split('=')
                jobstate = splitStr[-1].lstrip().rstrip()
                if ( 'R' in jobstate ):
                    jobInQueue = True
                elif ( 'Q' in jobstate ):
                    jobInQueue = True
                elif ( 'C' in jobstate ):
                    jobInQueue = True
                else:
                    jobInQueue = False
        
        return jobInQueue



    def isJobDone(self, jobdir, jobscript):
        """
        Finds 'code specific' finish condition for a particular run.
        Assumes that call is made from correct run directory
        
        This method can be redefined for other 'finish' criteria
        It must be defined as
        
        def isJobDone(jobdir, jobscript):
        ...
        ...
        return isDone (True/False)
    
        Note: This call is defaulting to a job NOT being finished if
        a specific isJobDone is not defined for a given applic.
        """

        print "No JobStatus:isJobDone method defined... set job to NOT done"
        return False



class FileEditorWithTags:
    """
    Creates objects capable of editing files by replacing string tags
    with values. The correspondence between these tags and values are
    specified by the dictionary objects:
    
     'tagValDict'    --> actual dictionary used to do string replacements
     'tagValDictTpl' --> dictionary template used by class method to automatically
    create new file edit objects with one parameter "sweeping" over a series of specified values
    """
    

    def __init__(self, tagValDictTpl):
        """
        Constructor (set base level input file values)
        These are used as template on which to do sweeps
        """
        self.tagValDictTpl = tagValDictTpl
        self.tagValDict    = tagValDictTpl


    def __del__(self):
        """
        Destructor
        """
        pass


    def setDictValue(self, str):
        """
        Setter for dictionary values from assignment string
        eg obj.setDictValue('EDIT_x=0.30')
        """
        
        splitStr=str.split('=')
        tag = splitStr[0].strip(' ')
        val = splitStr[1].strip(' ')

        if (tag not in self.tagValDict):
            print "tag ", tag, "not found in setDictValue(...)"
            sys.exit(1)
            
        try:
            # Try to convert val to number
            self.tagValDict[tag] = float(val)
        except:
            # If fails assume that val is a string
            self.tagValDict[tag] = val


    def setDictTplValue(self, str):
        """
        Setter for TEMPLATE dictionary values from assignment string
        eg obj.setDictTplValue('EDIT_x=0.30')
        
        NOTE: when one wants a sweep over a different parameter
        """
        
        splitStr=str.split('=')
        tag = splitStr[0].strip(' ')
        val = splitStr[1].strip(' ')

        if (tag not in self.tagValDictTpl):
            print "tag ", tag, "not found in setDictValue(...)"
            sys.exit(1)
            
        try:
            # Try to convert val to number
            self.tagValDictTpl[tag] = float(val)
        except:
            # If fails assume that val is a string
            self.tagValDictTpl[tag] = val


    def setDict(self, tagValDict):
        """
        Setter for new dictionary
        """
        self.tagValDict = tagValDict


    def setDictTpl(self, tagValDictTpl):
        """
        Setter for new dictionary template
        """
        self.tagValDictTpl = tagValDictTpl


    def show(self):
        """
        Print the current dictionary to be used
        to edit a file
        """
        print "tag edit info: ", self.tagValDict


    def fileLines(self, fname):
        """
        Count lines in file
        """

        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1


    def getParamSweeps(self, tag, params):
        """
        This is a factory which creates copies of this object
        Return a list of copied file edit objects generated by sweeping over
        'param' values for 'tag' in the dictionary template.
        """

        # Blank list of copied file objects
        fEditObjList=list()

        for v in params:
            fEditObj = FileEditorWithTags(self.tagValDictTpl) # Create new file edit object
            tagValDict = copy.copy(self.tagValDictTpl)        # copy dictionary from template
            tagValDict[tag] = v                               # edit the dictionary for new value
            fEditObj.setDict(tagValDict)                      # set new dictionary in copied file edit object
            fEditObjList.append(fEditObj)                     # Store fedit obj in list

        return fEditObjList


    def editFile(self, inFile):
        """
        Take each line in file and replace tagged strings
        with values to be set. tagValDict is a dictionary eg
        
        tagValDict['EDIT_x1' ] = 0.232
        tagValDict['EDIT_x2']  = 134
        tagValDict['EDIT_x3']  = 0.456
        ....
        """

        # Local copy for replace dictionary
        tagValDict = self.tagValDict
        
        outFile='tmp.in'
        f=open(inFile,  'r+')
        g=open(outFile, 'w')

        # Loop over all lines in file
        # SWS: can simplify the 'in' syntax
        for x in range(self.fileLines(inFile)):
            strline = f.readline()

            # Loop over tags and replace each
            for tag, val in tagValDict.iteritems():
                strline=strline.replace( tag, str(val) )
            g.write(strline)

        # Close local files
        f.close()
        g.close()

        # Copy over template file w/new file
        shutil.copy2(outFile, inFile)
        os.remove(outFile)




class Geometry:
    """
    Collection of geometric/math methods
    """

    def __init__(self, debugFlag=False):
        """
        Constructor
        """

        self.debugFlag = debugFlag
        if (self.debugFlag):
            print " "
            print "Geometry class created"
            print " debugFlag set to --> ", debugFlag
            print " "


    def __del__(self):
        """ Destructor """

        if (self.debugFlag):
            print " "
            print "Geometry class closed"


    def rdist(self, pt1, pt2):
        """
        Cartestian dist between pt1-pt2 (in list [] format)
        
        Args:
            pt1 : python list [x,y,z]
            pt2 : python list [x,y,z]
        """
        r2 = self.rdistSqr(pt1, pt2)
        return math.sqrt(r2)


    def rdistSqr(self, pt1, pt2):
        """
        Squared Cartestian dist between pt1-pt2 (in list [] format)
        
        Args:
            pt1 : python list [x,y,z]
            pt2 : python list [x,y,z]
        Returns:
            float value for cartestian distance
        """

        pt1 = map(float,pt1)   # Ensure args are type float
        pt2 = map(float,pt2)

        delR= map(float.__sub__, pt1, pt2)   # Subtracts coord lists elementwise
        dR2 = map(pow, delR, [2]*len(delR))  # Square each element
        return sum(dR2)                             # Sum elements


    def rdistPBCbox(self, pt1, pt2, bsize):
        """
        Cartestian dist between pt1-pt2 (in list [] format)
        Finds minimun distance assuming points are located inside a rectangular
        box with sizes in bsize = [Lx, Ly, Lz]
        
        Args:
            pt1 : python list [x,y,z]
            pt2 : python list [x,y,z]
            bsize : python list [Lx, Ly, Lz]
        Returns:
            float value for cartestian distance in periodic box
        """

        r2 = self.rdistPBCboxSqr(pt1, pt2, bsize)
        return math.sqrt(r2)


    def rdistPBCboxSqr(self, pt1, pt2, bsize):
        """
        Squrared Cartestian dist between pt1-pt2 (in list [] format)
        Finds minimun distance assuming points are located inside a rectangular
        box with sizes in bsize = [Lx, Ly, Lz]
        
        Args:
            pt1 : python list [x,y,z]
            pt2 : python list [x,y,z]
            bsize : python list [Lx, Ly, Lz]
        Returns:
           float value for cartestian distance squared periodic box
        """

        pt1 = map(float,pt1)   # Ensure args are type float
        pt2 = map(float,pt2)
        bsize=map(float,bsize)

        diff    = map(float.__sub__, pt1, pt2)     # Subtracts coord lists elementwise
        adelR   = map(math.fabs, diff)             # Take abs val elem-wise
        delRBox = map(float.__sub__, bsize, adelR) # Subtract coord diff from box

        for i,val in enumerate(delRBox):    # Select smallest coord diff elem-wise
            minDel = min(adelR[i], val)     # Store min results in delR list
            adelR[i] = minDel               #

        dR2 = map(pow, adelR, [2]*len(adelR)) # Square each element
        return sum(dR2)                       # Sum elements and take root


    def rdistInBounds(self, pt1, pt2, rcutoff, bsize):
        """
        Only checks x-coord with cutoff to determine
        bound for Cartesian distance.
        """
        idir = 0
        dr  = math.fabs(pt1[idir] - pt2[idir])
        if (dr >= rcutoff):
            return False
        else:
            return True


    def rdistInBoundsPBC(self, pt1, pt2, rcutoff, bsize):
        """
        Only checks idir-coord with cutoff to determine
        bound for Cartesian distance.
        NOTE: TEsting
        """
        idir = 0
        dr  = math.fabs(pt1[idir] - pt2[idir])
        drL = bsize[idir] - dr
        if (dr >= rcutoff) and (drL >= rcutoff):
            return False

        idir = 1
        dr  = math.fabs(pt1[idir] - pt2[idir])
        drL = bsize[idir] - dr
        if (dr >= rcutoff) and (drL >= rcutoff):
            return False
        else:
            return True



class Histogram1D:
    """
    Collection of binning/histogram methods in 1D
    """

    def __init__(self, xmin, xmax, numbins, debugFlag=False):
        """
        Constructor
        """

        self.xmin = xmin
        self.xmax = xmax
        self.numbins = numbins
        self.weights = None
        
        self.debugFlag = debugFlag
        if (self.debugFlag):
            print " "
            print "1D histogram class created"
            print " debugFlag set to --> ", debugFlag
            print " "

        # Associated histogram parameters
        self.binsize  = (xmax-xmin)/float(numbins)

        # Positions of midpoints of each bin based on xmin/max
        self.binvals = []
        binpos = xmin + (self.binsize/2.0)
        for ibin in range(numbins):
            self.binvals.append(binpos)
            binpos = binpos + self.binsize


    def __del__(self):
        """
        Destructor
        """

        if (self.debugFlag):
            print " "
            print "Histogram1D class closed"


    def setWeights(self, wts):
        """        
        'weights' is the factor that bin counting is weighted by. If no
        weights are given then the default for all data is --> '1'. If
        weights are given is should be a list with length equal to len(data)
        """
        self.weights = wts


    def getHistogramNorm(self, data):
        """
        Returns tuple with bin positions and normalized histogram for all
        values in 'data', this is a driver for getBinning().

        NOTE: !!!!!!! ----- The exact meaning of this norm factor needs clarification

        Args:
            data (list) positions that will be binned given params used to construct this object
        Returns:
            tuple with (bin-values, histogram)
        """

        binCounts, binWeights = self._getBinning(data)
        totalCounts = sum(binCounts)
        norm = self.binsize*float(totalCounts)
        if (norm <= 0.0):
            norm = 1.0

        normHist = [x/norm for x in binCounts]
        return (self.binvals, normHist)


    def getHistogram(self, data):
        """
        Returns tuple with bin positions and raw binned histogram (raw counts) for all
        values in 'data', this is a driver for getBinning().

        Args:
            data (list) positions that will be binned given params used to construct this object
        Returns:
            tuple with (bin-values, raw histogram)
        """

        binCounts, binWeights = self._getBinning(data)
        return (self.binvals, binCounts)



    def getHistogramWeighted(self, data):
        """
        Returns tuple with bin positions and raw binned histogram (raw counts) for all
        values in 'data', this is a driver for getBinning().

        Args:
            data (list) positions that will be binned given params used to construct this object
        Returns:
            tuple with (bin-values, raw histogram)
        """

        binCounts, binWeights = self._getBinning(data)
        return (self.binvals, binWeights)



    def getHistogramWeightedAvePerBin(self, data):
        """
        Returns tuple with bin positions and binned histogram for all
        values in 'data', this is a driver for getBinning() using the total weight
        per bin / (number of counts) for each bin

        Args:
            data (list) positions that will be binned given params used to construct this object
        Returns:
            tuple with (bin-values, raw histogram)
        """

        binCounts, binWeights = self._getBinning(data)

        for ibin, ibinWeight in enumerate(binWeights):

            if binCounts[ibin] == 0:
               binWeights[ibin] = binWeights[ibin]
            else:
               binWeights[ibin] = binWeights[ibin]/binCounts[ibin]

        return (self.binvals, binWeights)



    def getBinSize(self):
        """
        Return calculated binsize
        """
        return self.binsize


    def _getBinning(self, data):
        """
        Uses the output from bin1DList to get the binning info
        (used to generate the histogram data). Calculates number of values
        in each bin. Calculates total weighted counts for each bin.

        Note: If weights not set then binCounts == binWeights

        Args:
            data (list) raw 1D list of data values
            
        Return: tuple of binCounts, binWeights. 
        """

        binList = self._bin1DList(data)
        binCounts = [0 for x in range(self.numbins)]
        binWeights = [0 for x in range(self.numbins)]

        # print "binList = ", binList, len(binList)
        # print "binCounts = ", binCounts, len(binCounts)

        if None in binList:
            print "_bin1DList call failed... values outside of limits"
            sys.exit(3)

        if (self.weights != None):
            if (len(data) != len(self.weights)):
                print "getBinning: len(weights) != to len(data)"
                sys.exit(3)

        # If weights have been set...
        if (self.weights == None):
            for ibin in binList:
                binCounts[ibin] = binCounts[ibin] + 1
                binWeights[ibin] = binWeights[ibin] + 1
        else:
            for i,ibin in enumerate(binList):
                count = self.weights[i]
                binCounts[ibin]  = binCounts[ibin] + 1
                binWeights[ibin] = binWeights[ibin] + count

        return binCounts, binWeights



    def _getBinVals(self):
        """
        Return position list of midpoints of bins
        """
        return self.binvals


    def _bin1D(self, val):
        """
        Bins returned run from 0-->(numbins-1)
        If value outside of range set in constructor then
        'None' is returned
        """

        # if ( (val <= self.xmin) or (val >= self.xmax) ):
        if ( (val < self.xmin) or (val > self.xmax) ):
            return None
        else:
            ibin = int( (val-self.xmin) / self.binsize )
            return ibin


    def _bin1DList(self, valList):
        """
        Driver for bin1D(val) for a list of values
        """
        return map(self._bin1D, valList)
