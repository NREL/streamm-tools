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

    # Count lines in file
    def fileLines(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def getIDs(self, inFile):

        outFile='tmp.in'
        f=open(inFile,  'r+')
        g=open(outFile, 'w')

        # Loop over all lines in file
        for x in range(self.fileLines(inFile)):
            strline = f.readline()
            cmdStr="qdel " + strline
            print "Executing ", cmdStr
            os.system(cmdStr)
            sleep(0.1)

        # Close local files
        f.close()
        g.close()


u=Misc()
u.getIDs("jobs.tst")
