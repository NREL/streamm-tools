#! /usr/bin/env python
"""
Subroutines for file operations
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# Email travis.kemper@nrel.gov
# Version 1.00 


def file_exists(filename):
    try:
        with open(filename) as f:
            return True
    except IOError:
        return False

def getlines(F_IN):
    F = open( F_IN , 'r' )
    Lines = F.readlines()
    F.close()
    return Lines
