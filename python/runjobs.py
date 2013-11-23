#!/usr/bin/env python

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



######################################################################
#
# Collection of useful methods
# 
######################################################################
class Misc:

    def __init__(self, debugFlag=False):
        self.debugFlag = debugFlag
        if (self.debugFlag):
            print " "
            print "Misc class created"
            print " debugFlag set to --> ", debugFlag
            print " "

    # Destructor (closes file)
    def __del__(self):
        if (self.debugFlag):
            print " "
            print "Misc class closed"

    ####################################################################
    #
    # Recursive glob using the os.walk command to traverse
    # directory structure. Starting from 'treeroot' directory,
    # return list of full path names for files containing 'pattern'
    #
    #  treeroot - Directory name to start top of search
    #  pattern  - String pattern (can use * wildcard) for search
    #
    ####################################################################
    def recursive_glob(self, treeroot, pattern):
        results = []
        for base, dirs, files in os.walk(treeroot):
            goodfiles = fnmatch.filter(files, pattern)
            results.extend(os.path.join(base, f) for f in goodfiles)
        return results


    ###############################################
    #
    # Generate list from x0, x0+dx,..... xmax
    #
    ###############################################
    def getList(self, x0, dx, xmax):
        xList=[]    
        xList.append(x0)
        xnext = x0+dx
        while xnext<xmax:
            xList.append(xnext)
            xnext = xnext+dx
        return xList

    #########################################################
    #
    # Return a string with leading 0's so file systems etc
    # sort is a normal way
    #
    #########################################################
    def getTagForSorting(self, num):
        if (num > 9999):
            print "getTagForSorting: tag failed. Job num > 9999"
            sys.exit(1)
        tagTemplate = 100000
        tagstr = str(tagTemplate+num)
        return tagstr.lstrip('1')

######################################################################


######################################################################
#
# Class for managing output
#
######################################################################
class IO:

    # fileName -- string for file name
    def __init__(self, fileName, debugFlag=False):
        self.fileName = fileName
        self.debugFlag = debugFlag
        self.fileObj = open(fileName, 'w')
        if (self.debugFlag):
            print " "
            print "IO file object ", fileName, " opened"
            print " debugFlag set to --> ", debugFlag
            print " Flags can be reset in scripts by developers"
            print " "

    # Destructor (closes file)
    def __del__(self):
        self.fileObj.close()
        if (self.debugFlag):
            print " "
            print "IO file object ", self.fileName, " closed"
            print " "

    def flush(self):
        self.fileObj.flush()

    # Print debug info if set
    def prtdb(self, prtString):
        if (self.debugFlag):
            print prtString

    def prtScreen(self, prtString):
        print prtString

    def prtFile(self, prtString):
        self.fileObj.write(prtString)
        self.fileObj.write(" \n")

    # Write to screen and file
    def prt(self, prtString):
        print prtString
        self.fileObj.write(prtString)
        self.fileObj.write(" \n")
######################################################################


#class FileEditorWithSomething:
#    def __init__(self):
#        pass

#######################################################################
#
# Job class representing data/members needed
# to represent an instance of a code that will run
# with a single set of input parameters
#
#######################################################################
class Job:

    # Constructor
    def __init__(self, repopath, exe,
                 inputfile, script, verbose=False):

        self.repopath = repopath     # Path to repo for common files
        self.exe = exe               # Name of executable binary
        self.inputfile = inputfile   # Name of input file
        self.script = script         # Name of script for PBS
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
    def setScript(self, script):
        self.script = script

    # Expects an editor object that implements the
    # appropriate 'editFile(..)' method
    def setInputEditor(self, editor):
        self.inputEditorSet = True
        self.inputEditor = editor

    # Expects an editor object that implements the
    # appropriate 'editFile(..)' method
    def setScriptEditor(self, editor):
        self.scriptEditorSet = True
        self.scriptEditor = editor
        
    # Dump state of object
    def info(self):
        print " "
        print "Data members set to: "
        print "   repopath  = ", self.repopath
        print "   rundir    = ", self.rundir
        print "   exe       = ", self.exe
        print "   inputfile = ", self.inputfile
        print "   script    = ", self.script
        print " "

    # Create run directory and populate with appropriate files
    def createJobRun(self):

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


    # Use imported edit methods to change input file for run
    # The imported method must implement an 'editfile(file)' method
    def changeInputFile(self):

        if (self.inputEditorSet):
            currFile=os.path.join(self.rundir, self.inputfile)
            self.inputEditor.editFile(currFile)
        else:
            print "Input file editor not set. ", self.inputfile, " is unchanged"
            sys.exit(1)

    # Use imported edit methods to change script file for run
    # The imported method must implement an 'editfile(file)' method
    def changeScriptFile(self):

        if (self.scriptEditorSet):
            currFile=os.path.join(self.rundir, self.script)
            self.scriptEditor.editFile(currFile)
        else:
            print "Scriptfile editor not set. ", self.script, " is unchanged"
            sys.exit(1)
#######################################################################################


#####################################################################################
#
# Creates objects capable of editing files by replacing string tags
# with values. The correspondence between these tags and values are
# specified by the dictionary objects:
#            
#     'tagValDict' --> actual dictionary used to do string replacements
#  'tagValDictTpl' --> dictionary template used by class method to automatically
#                      create new file edit objects with one parameter
#                      "sweeping" over a series of specified values            
#
#####################################################################################            
class FileEditorWithTags:
    
    # Constructor (set base level input file values)
    # These are used as template on which to do sweeps
    def __init__(self, tagValDictTpl):
        self.tagValDictTpl = tagValDictTpl
        self.tagValDict    = tagValDictTpl

    # Destructor
    def __del__(self):
        pass

    # Setter for dictionary values from assignment string
    # eg obj.setDictValue('EDIT_x=0.30')
    def setDictValue(self, str):
        
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

    # Setter for TEMPLATE dictionary values from assignment string
    # eg obj.setDictTplValue('EDIT_x=0.30')
    #
    # NOTE: when one wants a sweep over a different parameter
    def setDictTplValue(self, str):
        
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


            
    # Setter for new dictionary
    def setDict(self, tagValDict):
        self.tagValDict = tagValDict

    # Setter for new dictionary template
    def setDictTpl(self, tagValDictTpl):
        self.tagValDictTpl = tagValDictTpl

    # Print the current dictionary to be used
    # to edit a file
    def show(self):
        print "tag edit info: ", self.tagValDict

    # Count lines in file
    def fileLines(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    #########################################################################
    #
    # This is a factory which creates copies of this object
    # Return a list of copied file edit objects generated by sweeping over
    # 'param' values for 'tag' in the dictionary template.
    #
    #########################################################################
    def getParamSweeps(self, tag, params):

        # Blank list of copied file objects
        fEditObjList=list()

        for v in params:
            fEditObj = FileEditorWithTags(self.tagValDictTpl) # Create new file edit object
            tagValDict = copy.copy(self.tagValDictTpl)        # copy dictionary from template
            tagValDict[tag] = v                               # edit the dictionary for new value
            fEditObj.setDict(tagValDict)                      # set new dictionary in copied file edit object
            fEditObjList.append(fEditObj)                     # Store fedit obj in list

        return fEditObjList

    ###########################################################
    #
    # Take each line in file and replace tagged strings
    # with values to be set. tagValDict is a dictionary eg
    #
    #    tagValDict['EDIT_x1' ] = 0.232
    #    tagValDict['EDIT_x2']  = 134
    #    tagValDict['EDIT_x3']  = 0.456
    #    ....
    #
    ###########################################################
    def editFile(self, inFile):

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
#######################################################################################
