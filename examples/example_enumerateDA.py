#!/usr/bin/env python

import os, sys

def main():
    """
    Example of running donor-acceptor systems through the enumerate.py script.
    Needs a 'subset.test' file and access to the BuildingBlocks repos.
    """

    print " "
    print "Removing old mols-* directory"
    os.system("rm -rf example_enumerateDA-mols-dir")

    print "Running enumerate.py (donoracceptorsystems) test"
    print " "
    run_py = "python ../da_builder/enumerate.py -d -c DA -b ../../BuildingBlocks-private -s example_enumerateDA-subset.txt -f example_enumerateDA-mols-dir --opvPath=../.. --repoPath=./test-repo > example_enumerateDA.dat"
    os.system(run_py)

    # Partial output for diff-ing
    # os.system("tail -4 example_enumerateDA.dat")

    os.system("mv example_enumerateDA.dat example_enumerateDA-mols-dir")
    print "Results for builder in example_enumerateDA-mols-dir"
    print " "


if __name__ == '__main__':
    main()
