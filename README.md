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

********************************
Development/Release workflow
********************************

1. Develop in a branch

2. Merge necessary branches into master

3. Test master and update docs

4. Once ready for a new release, tag the master branch (version number is actual example)
   *  git tag -a v1.1.4 -m "Give release a descriptive name"
   *  git push origin --tags

5. Push master and tags for repo to external NREL github
   *  git remote add 'name-of-repo' https://github.com/NREL/'name-of-repo'.git
   *  git push -u 'name-of-repo'  master
   *  git push -u 'name-of-repo'  v1.1.4
   *  git remote rm 'name-of-repo'

6. Making documentation available
   * Develop docs in master (as always)
   * Assuming step 5 complete and streamm-tools repo has been updated with tagged version
   * Email Toby Walhers to git pull from the tagged version on the external site
   * Toby will build and publish docs
