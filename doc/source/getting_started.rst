.. _getting_started:

.. index:: regression tests, unit tests, test scripts, github

*********************************
Getting started
*********************************

.. _download-from-github:

Clone Github repository 
============================

The STREAMM *tools* code is maintained in a Github repository at https://github.com/NREL/streamm-tools hosted by
the National Renewable Energy Lab (NREL). The main STREAMM repo contains many projects that use the tools code as 
well as the separate suite of unit/regression tests. The web based interface allows one to download the tools repo
using the GUI. The tools repo can be downloaded using the linux command line ::

     git clone https://github.com/NREL/streamm-tools.git

The tools repo contains the following directories and files

- README.md  -- Repo notes
- analysis   -- 3rd party codes used for post-processing
- config.sh  -- Configuration script that sets PYTHONPATH and other needed environment variables
- da_builder -- Donor-Acceptor builder source code
- doc        -- Sphinx documentation and scripts for creating Python API from docstrings
- examples   -- High-level specific examples using the tools/scripts and tools/src code (documented in /doc)
- scripts    -- High-level drivers using tools/src code
- src        -- Main classes implementing the STREAMM tools functionality




..  
   _install-tools:

    Installation
   ============================

    The curent release of STREAMM *tools* is a pure python based code with
    minimal dependencies. Therefore, no compolation is
    necessary. 

..  _configure-tools:


Configure tools
============================

The config.sh script is provided to set various environment variables
and the PYTHONPATH required for the source code modules to run
correctly. Before using STREAMM tools or running the tests execute ::

    source config.sh

This will set the correct environment for the current terminal
session. To set permanently, note the output of the script and then
set in your local environment. 


Running tests
========================

The tests for the STREAMM tools is in a separate Github repository https://github.com/NREL/streamm-tools-tests hosted by the National
Renewable Energy Lab (NREL). The tools repo can be downloaded using the linux command line ::

     git clone https://github.com/NREL/streamm-tools-tests.git

If the tools repo has been configured correctly the check.sh script
can be used to run the units tests. The help info for check.sh is
output from ::

    check.sh

To run a particular test ::

    check.sh 'test-name.py'

and to run all tests ::

    check.sh all



Test descriptions
========================

.. automodule:: test_angleContainer
   :members: main

.. automodule:: test_bondContainer
   :members: main

.. automodule:: test_checkTypes
   :members: main

.. automodule:: test_dihedralContainer
   :members: main

.. automodule:: test_dumpLammps
   :members: main

.. automodule:: test_periodictable
   :members: main

.. automodule:: test_ffparameters
   :members: main

.. automodule:: test_n2_mpiBase
   :members: main

.. automodule:: test_nX_mpiBase_splitData
   :members: main

.. automodule:: test_particleConstructors
   :members: main

.. automodule:: test_particleContainer
   :members: main

.. automodule:: test_particleSetInfo
   :members: main

.. automodule:: test_ptclContainerConstructor
   :members: main

.. automodule:: test_searchTags
   :members: main

.. automodule:: test_simulation
   :members: main

.. automodule:: test_strucAdd
   :members: main

.. automodule:: test_strucAdd2
   :members: main

.. automodule:: test_strucAddBig
   :members: main

.. automodule:: test_strucCompressID
   :members: main

.. automodule:: test_strucCompressIDWithAngles
   :members: main

.. automodule:: test_strucDumpSave
   :members: main

.. automodule:: test_strucEmpty
   :members: main

.. automodule:: test_strucEmpty2
   :members: main

.. automodule:: test_strucSetPtclPos
   :members: main

.. automodule:: test_strucWithAngles
   :members: main

.. automodule:: test_strucWithDihedrals
   :members: main

.. automodule:: test_subStructure
   :members: main

.. automodule:: test_subStructureWithAngles
   :members: main

.. automodule:: test_subStructureWithDihedrals
   :members: main
