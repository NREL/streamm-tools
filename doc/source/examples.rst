
.. _examples:

*************************
Examples
*************************

To run these examples your PATH needs to include `tools/da_builder/`
and your PYTHONPATH needs to include `tools/src/`.  See :ref:`configure-tools <configure-tools>`

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

These examples create molecular structures files, and input files for
`Gaussian <http://www.gaussian.com/>`_ , `GROMACS
<http://www.gromacs.org/>`_  and `LAMMPS <http://lammps.sandia.gov/>`_
by reading in reference structures from  the `BuildingBlocks`
repository and reference force-field parameters from `conj.itp`.
The `BuildingBlocks` repository clonned manually::

    git clone https://github.com/NREL/streamm-BuildingBlocks

or using ::
     
     ./examples_setup.sh

which will also check that the `PATH` and `PYTHONPATH` enviromental variables were set correctly
by the :ref:`configure-tools <configure-tools>` script. 

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
