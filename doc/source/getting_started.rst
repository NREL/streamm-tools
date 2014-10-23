.. _getting_started:

*******************************
Getting starting
*******************************



.. _download-from-github:

Download from Github:
============================

The STREAMM *tools* code is maintained in a Github repository at https://github.nrel.gov/streamm/tools hosted by the National Renewable Energy Lab (NREL). The main STREAMM repo contains many projects that use the tools code as well as the separate suite of unit/regression tests. The web based interface allows one to download the tools repo using the GUI. The tools repo can be downloaded using the linux command line ::

     git clone https://github.nrel.gov:streamm/tools.git

The tools repo contains the following directories and files

- AtomicPy      -- legacy python files (converting to src)
- README.md -- Repo notes
- config.sh      -- Configuration script that sets PYTHONPATH and other needed environment variables
- analysis        -- Code for analysis, visualization
- doc               -- Sphinx documentation and scripts for creating Python API from docstrings
- examples     -- Examples of using STREAMM tools for specific applications
- simmod        -- ?
- src                -- Main classes implementing the STREAMM tools functionality



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



