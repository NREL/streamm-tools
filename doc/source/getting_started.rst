.. _getting_started:

.. index:: regression tests, unit tests, test scripts, github

*********************************
Installation / Getting starting
*********************************

.. _download-from-github:

Download from Github
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


..  _configure-tools:

Configure tools
============================

The config.sh script is provided to set various environment variables
and the PYTHONPATH required for the source code modules to run
correctly. Before using STREAMM tools or running the tests execute ::

    source config.sh

This will set the correct environment for the current terminal
session. To set permanently, note the output of the script and then
set in your local environment. A number of packages are required for
the complete functionality of the STREAMM tools.
A current list of python module dependencies on the Peregrine cluster
at NREL (our default test platform) is as follows

Cython==0.20.1
GridDataFormats==0.2.4
MDAnalysis==0.8.1
MDAnalysisTests==0.8.1
backports.ssl-match-hostname==3.4.0.2
biopython==1.64
decorator==3.4.0
h5py==2.2.1
ipython==2.0.0
matplotlib==1.3.1
memory-profiler==0.31
mpi4py==1.3.1
networkx==1.9
nose==1.3.1
numpy==1.8.1
pandas==0.13.1
psutil==2.1.1
psycopg2==2.5.3
pyparsing==2.0.1
python-dateutil==2.2
pytz==2014.2
pyzmq==14.1.1
requests==2.3.0
scikit-learn==0.14.1
scipy==0.13.3
six==1.6.1
tornado==3.2
virtualenv==1.11.4
wsgiref==0.1.2

Python 2.7.3 or greater should be installed to run STREAMM. Development and testing of the STREAMM code has been performed with
Python 2.7.3 on the Peregrine cluster at NREL and Python 2.7.7 on Mac OSX. We recommend not using Python 3.x and greater as no
testing has been done for those versions. 

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
