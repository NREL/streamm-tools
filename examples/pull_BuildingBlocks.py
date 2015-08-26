#!/usr/bin/env python
# -*- coding: utf-8 -*-
# encoding: utf-8
"""externals

Created by kurt on 2014-07-23
Edited by Scott Sides between 2014-08-01 and 2015-08-31



About
========
Add the relative target path and the git url to the `externals` list,
one for each item you want checked-out.

"""

import sys, os, traceback, optparse
import time
import re
import shutil
import subprocess
from optparse import OptionParser


#############################################################################
# Command line option parse

parser = OptionParser()
usage = """

%prog [options]

     Manages pulling external git repos into project
     Default is to pull from the github.com/NREL external site

     Place in your git repo root, edit the external list, then run the command: externals

     Version flags are (default is external_git)
       * internal_git
       * external_git
"""
parser.set_usage(usage)

parser.add_option("--version",
                  dest="versionStr",
                  default='external_git',
                  type="string",
                  help="Chooses internal/external github versions of repo. Release versions should be external \n")

# Acessing options
(options,args) = parser.parse_args()
versionStr     = options.versionStr
#############################################################################



NRELinternal_externals = [
('./BuildingBlocks-release','git@github.nrel.gov:streamm/BuildingBlocks-release.git'),
]

NRELexternal_externals = [
('./BuildingBlocks-release','http://github.com/NREL/streamm-BuildingBlocks-release.git'),
]

#
# Choosing version repo flag
#
if ( versionStr == 'internal_git' ):
    externals = NRELinternal_externals
elif ( versionStr == 'external_git' ):
    externals = NRELexternal_externals
else:
    print "Repo version string not recognized"
    sys.exit(0)


def pull_externals():
    this_cwd = os.getcwd()
    for i in externals:
        loc, repo = i[0], i[1]
        target_path = '%s/%s' % (this_cwd, loc)
        if not os.path.exists(target_path):
            cmd = "git clone %s %s" % (repo, target_path)
            ex = subprocess.Popen(cmd, cwd=this_cwd, shell=True)
            output = ex.communicate()[0]
            print(output)
        else:
            cmd = "git pull"
            print "Repo ", repo
            ex = subprocess.Popen(cmd, cwd=target_path, shell=True)
            output = ex.communicate()[0]
            print(output)
            print " "

def twiddle_ignore():
    for i in externals:
        loc = i[0]
        with open('.gitignore', 'a+') as f:
            if not any(loc == x.rstrip('\r\n') for x in f):
                f.write(loc + '\n')
                print("\n\n + added %s to .gitignore" % loc)

def main (restore=False):
    try:
        pull_externals()
    except:
        print "pull_externals failed in tools-opv/pull_externals"
        print "Check ssh keys and re-generate if necessary"
    twiddle_ignore()

if __name__ == '__main__':
    main()
