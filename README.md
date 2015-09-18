tools
=====

Common Tools for all STREAMM projects
Note:
  Jan  23 2014 -- Moved job creation/submission/monitoring python scripts to tools-public
  July 08 2015 -- Big reorg/cleaning in preparation for release


*****************
File structure:
*****************

- README.md  -- Repo notes
- analysis   -- 3rd party codes used for post-processing
- config.sh  -- Configuration script that sets PYTHONPATH and other needed environment variables
- da_builder -- Donor-Acceptor builder source code
- doc        -- Sphinx documentation and scripts for creating Python API from docstrings
- examples   -- High-level specific examples using the tools/scripts and tools/src code (documented in /doc)
- scripts    -- High-level drivers using tools/src code
- src        -- Main classes implementing the STREAMM tools functionality



Source code
===========

Source code is hosted in a git repository at

https://github.com/NREL/streamm-tools

and is available under the Apache License, Version 2.0 (see
the file LICENSE.txt).


Configure 
===========

The config.sh script is provided to set various environment variables
and the PYTHONPATH required for the source code modules to run
correctly. Before using STREAMM tools or running the tests execute 

source config.sh

********************************
Development/Release workflow
********************************

1. Develop in a branch

2. Merge necessary branches into master

3. Test master and update docs

4. Once ready for a new release, tag the master branch
   *  git tag -a v1.1.4 -m "Give release a descriptive name"
   *  git push origin --tags

5. Push master and tags for repo to external NREL github
   *  git remote add streamm-tools  https://github.com/NREL/streamm-tools.git
   *  git push -u streamm-tools  master
   *  git push -u streamm-tools  v1.1.4
   *  git remote rm -u streamm-tools
