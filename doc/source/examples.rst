
.. _examples:

*************************
Examples
*************************

Selected tools-tests
========================================================
These examples are closely related to regression tests in the tools-tests repo. Each of the selections below is available as
an IPython notebook, so you may explore these examples further. After clicking on the example, there is a download link. Click on the
download and drag the desired file to a location on your computer. The IPython web interface for a notebook can be started from
the command line by e.g. ::

   ipython notebook example_test.ipynb

.. toctree::
   :maxdepth: 1

   example_particleContainer.rst
   example_subStructure.rst
   example_searchTags.rst
   example_strucAdd.rst

.. Removing this code since OPV-database description is out for release
.. More complicated example with simulation objects
.. ========================================================
.. Functional/basis 'spamming' description (SCOTT)


Molecular generation
========================================================

Examples of running the donoracceptorsystems.py script to build complex molecules.
Bash scripts of these examples are in ::

   /tools/examples/mol_gen_ex1.sh
   /tools/examples/mol_gen_ex2.sh
   /tools/examples/mol_gen_ex3.sh

To run these examples, various templates from the OPV project must be available as well as files describing
the structure of molecular building blocks (see below). To set these paths and check/copy the appropriate
files into the examples directory, a user can run ::

   mol_gen_setup.sh

which will indicate whether the necessary paths/files can be set.



Files and directories
-------------------------------------------------------

- tools/da_builder/donoracceptorsystems.py: Concatenates molecular cply files based on connectivity tags 
- BuildingBlocks repos: must contain the directories

   - acceptors   
   - donors            
   - functional_groups spacers           
   - terminals

The backbone of the molecule can be composed of  acceptors, spacers or
donors, and can be  decorated with molecules from the
functional_groups directory and terminated with molecules from the
terminals directory. Run ::

	python ../da_builder/donoracceptorsystems.py -h

to see a list of the options  


Donor-Acceptor Builder tutorials
-------------------------------------------------------

.. toctree::
   :maxdepth: 1

   mol_gen_ex1.rst
   mol_gen_ex2.rst
   mol_gen_ex3.rst



Molecular dynamics
========================================================

Create input files for molecular dynamics runs based the Molecular
generation example::
   
   /tools/examples/mol_md.sh

This script acts on example 3 of the Donor-Acceptor Builder
tutorials, by setting some local variables:: 

    MOLDIR=mols/thiophene_R_hexane_/
    MOLID=acc1_thiophene_R_hexane__n5

which is the P3HT example, but any structure can be used. 


Molecular dynamics file generation examples 
-------------------------------------------------------

.. toctree::
   :maxdepth: 1

   xyz2gromacs.rst
   xyz2data.rst
   replicate_data.rst

Coming soon
========================================================

- analyize MD output 
