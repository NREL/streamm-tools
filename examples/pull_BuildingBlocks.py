#!/usr/bin/env python
# -*- coding: utf-8 -*-
# encoding: utf-8
"""externals

Created by kurt on 2014-07-23

USAGE
========
Place in your git repo root
Edit the external list
Then run the command:
  externals


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


externals = [
('./BuildingBlocks-release','git@github.nrel.gov:streamm/BuildingBlocks-release.git'),
]


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
