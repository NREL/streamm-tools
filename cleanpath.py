#!/usr/bin/env python

import os, sys, math, random, time

def rmFromEnvironPath(pathVar, rmSearchList):
    """
    Search path variable and remove paths if any part contains the search string
    Writes an external script with correct shell commands and executes. Order of
    path variable is preserved for remaining paths

    Args:
        pathVar      (str) : name of environment variable (eg PYTHONPATH)
        rmSearchList (list): list of search strings for this is any part of path
    """

    newPath=str()                        # New path variable
    scriptName="clean-"+ pathVar + ".sh" # Unique script name
    path = os.getenv(pathVar)            # Get path from shell
    if (not path==None):
        pathList = path.split(':')       # Split path into list elements
    else:
        pathList = []


    # Construct list with search elements removed
    newPathList=[]
    for p in pathList:
        if any(x in p for x in rmSearchList):
            print "Removing ", p, " from ", pathVar
        else:
            newPathList.append(p)

    for p in newPathList:           # Construct path setting from new path list
        newPath += p + ":"          # Build up path
    newPath=newPath[:-1]            # Remove trailing ':'

    #
    # Construct shell script with new path variable set
    # and the appropriate EXPORT command
    #
    fObj=open(scriptName, 'w')             # If file doesnt exists create
    newPath=pathVar + "=" + newPath        # Construct new path wo/appending old var
    fObj.write(newPath + '\n')
    fObj.write('echo "Sourcing script ' + scriptName + '" \n')
    fObj.write('export ' + pathVar + '\n')
    fObj.close()

    # print "newPath---> ", newPath, "\n"
    return


#
# This instance will append environment commands to same cleanPath.sh file
# 
rmFromEnvironPath('PYTHONPATH', ['/tools/src', '/tools-', './'] )
rmFromEnvironPath('PATH',       ['/tools/src', 'tools-', './', '/scripts'] )
rmFromEnvironPath('TOOLS_PATH', ['/tools'] )
