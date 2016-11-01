#! /usr/bin/env python
"""
This will clean up files created during testing
the test directory
"""

import os,shutil

keepexts = ['.py','itp']
for subdir, dirs, files in os.walk(os.path.dirname(__file__)):
    for file in files:
        if( file[-3:] not in  keepexts):
            print "File to remove ",file
            os.remove(os.path.join(os.path.dirname(__file__), file))
    for dir in dirs:
        print dir
        shutil.rmtree(os.path.join(os.path.dirname(__file__), dir))
