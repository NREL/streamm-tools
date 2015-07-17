.. _examples:

*************************
Examples
*************************

Selected tools-tests
========================================================
These examples are closely related to regression tests in the
tools-tests repo. 

.. toctree::
   :maxdepth: 1

   example_particleContainer.rst
   example_subStructure.rst


IPython notebook downloads
========================================================
Each of the selections above is available as an IPython notebook, so you may explore these examples further.
After clicking on the download link, drag the desired file to a location on your computer.

:download:`example_subStructure.ipynb      <example_subStructure.ipynb>`

:download:`example_particleContainer.ipynb <example_particleContainer.ipynb>`


The IPython web interface for a notebook can be started from the command line by e.g. ::

   ipython notebook example_test.ipynb

.. Removing this code since OPV-database description is out for release
.. More complicated example with simulation objects
.. ========================================================
.. Functional/basis 'spamming' description (SCOTT)


Molecular generation
========================================================

Examples of running the donoracceptorsystems.py script to build complex molecules 

A bash file of these examples in ::
   
   /tools/examples/mol_gen.sh 


Files and directories
-------------------------------------------------------

- donoracceptorsystems.py: Concatenates molecular cply files based on connectivity tags 
- BuildingBlocks: must contain the directories 

   - acceptors   
   - donors            
   - functional_groups spacers           
   - terminals

The backbone of the molecule can be composed of  acceptors, spacers or
donors, and can be  decorated with molecules from the
functional_groups directory  and terminated with molecules from the
terminals directory. Run ::

	python ../da_builder/donoracceptorsystems.py -h

to see a list of the options  


Donor-Acceptor Builder tutorials
-------------------------------------------------------

.. toctree::
   :maxdepth: 1

   example_da_1.rst
   example_da_2.rst
   example_da_3.rst



Coming soon
========================================================

- running MD sim (LAMMPS, Gromacs)
- calculation of RDFsd

- setup of system inputs
- running MD sim (LAMMPS, Gromacs)
- calculation of RDFs
