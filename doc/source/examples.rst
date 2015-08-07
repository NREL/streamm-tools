
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


Molecular generation and replication 
========================================================


To run these examples the reference coordinate files in the
BuildingBlocks-release repo are needed, your PATH needs to include
`tools/da_builder/`  and your PYTHONPATH needs to include `tools/src/`
needs to set correctly to include . By running::

   source config.sh

in the tools directory and::

   examples_setup.sh

in the directory where you are running the examples. The functionality
of :ref:`donoracceptorsystems.py<donoracceptorsystems>` will be tested
and the BuildingBlocks-release repo will be cloned in your current directory. 

.. toctree::
   :maxdepth: 1

   benzene.rst
   C60.rst
   thiophene.rst
   oligothiophene.rst
   oligomethylthiophene.rst
   P3HT.rst
   replicateP3HT.rst
   replicateP3HTC60benzene.rst
