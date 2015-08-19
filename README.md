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

Test line
