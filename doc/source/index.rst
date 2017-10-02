.. highlight:: rst


.. ::
    
    # with overline, for parts
    * with overline, for chapters
    =, for sections
    -, for subsections
    ^, for subsubsections
    ", for paragraphs


STREAMM
#######

 The Simulation Toolkit for Renewable Energy and Advanced Materials Modeling (STREAMM) is a python package that generates structures, input files for quantum chemical and molecular dynamics codes.
 STREAMM does not directly conduct simulations rather it is meant to drive quantum chemical and molecular dynamics codes to allow for high-throughput computational analysis of materials.
 The streamm package is written for `python2.7 <https://www.python.org/download/releases/2.7/>`_ and was written using the core modules of the `pymatgen <http://pymatgen.org/>`_ code as the initial core modules.


:STREAMM: Copyright (C) 2015, Dr. Scott W. Sides, Dr. Travis W. Kemper, Dr. Ross E. Larsen and Dr. Peter Graf.
:Website: `<http://github.com/NREL/streamm-tools>`_

Quick install
*************

Install the python package with pip::

    $ pip install streamm

Access the streamm modules with

.. code:: python

    import streamm
    
This will allow access to the streamm modules.

Contents
********

.. toctree::
    :maxdepth: 2
     
    installation_instructions.rst 
    getting_started.rst
    functionality.rst
    how_to.rst
    examples.rst
    api_documentation.rst
    publication_highlights.rst
    
Release Notes
*************

v0.3.1 -- September 2017
========================

* Add pymatgen (https://github.com/materialsproject/pymatgen) dependency 
* Create tests for each module in the directory tests/
* Split up modules into directories
* Move functions dependent on mpi to util directory

v0.3.0 -- August 2017
======================

* Update the structure of the code to allow `setup.py` installation 


v0.2.0 -- August 28 2015 
========================

* Initial release

`NREL <http://www.nrel.gov/>`_ is a National Laboratory of the U.S. Department of Energy,
Office of Energy Efficiency and Renewable Energy, operated by the Alliance for Sustainable Energy, LLC.

Licenses
========

Licensed under the Apache License, Version 2.0

.. toctree::
    :maxdepth: 2
     
    license.rst 

Referencing STREAMM
*******************

When referencing the STREAMM toolkit in publications, this website can be cited as::

  Dr. Scott W. Sides, Dr. Travis W. Kemper, Dr. Ross E. Larsen and Dr. Peter Graf. "STREAMM (Simulation Toolkit for
  Renewable Energy and Advanced Materials Modeling)," National Renewable Energy Lab, 21 Sept. 2015. <http://github.com/NREL/streamm-tools>.
  
  
Also reference the Materials genome project code pymatgen::
    
    Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier, Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A. Persson, Gerbrand Ceder. Python Materials Genomics (pymatgen) : A Robust, Open-Source Python Library for Materials Analysis. Computational Materials Science, 2013, 68, 314–319. doi:10.1016/j.commatsci.2012.10.028
