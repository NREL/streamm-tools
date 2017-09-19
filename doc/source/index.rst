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

The Simulation Toolkit for Renewable Energy and Advanced Materials Modeling (STREAMM) is a python module that generates structures,
input files for quantum chemical and molcular dynamics codes.
STREAMM does not directly conduct simulations, rather it is ment to drive and connect
quantum chemical and molecular dynamics codes to allow for high-throughput compulational analysis of materials. 

:STREAMM: Copyright (C) 2015, Dr. Scott W. Sides, Dr. Travis W. Kemper, Dr. Ross E. Larsen and Dr. Peter Graf.
:Website: `<http://github.com/NREL/streamm-tools>`_

Quick install
*************

Install the python package with pip::

    $ pip install streamm

Access the streamm modules with::

    $ import streamm
    
This will allow access to the streamm modules.

Contents
********

.. toctree::
    :maxdepth: 2
     
    installation_instructions.rst 
    getting_started.rst
    api_documentation.rst
    publication_highlights.rst
    
Release Notes
*************

- v0.3.1 -- September 2017

    -- Split up modules into directories
    -- Create tests for each module in the directory tests/
    -- Add pymatgen (https://github.com/materialsproject/pymatgen) dependency 
    -- Move functions dependent on mpi to util directory

- v0.3.0 -- August 2017

    -- Major update to code structure allowing for `python setup.py` installation 


- v0.2.0 -- August 28 2015 -- Initial release

- `NREL <http://www.nrel.gov/>`_ is a National Laboratory of the U.S. Department of Energy,
  Office of Energy Efficiency and Renewable Energy, operated by the Alliance for Sustainable Energy, LLC.

- Licensed under the Apache License, Version 2.0


Referencing STREAMM
*******************

When referencing the STREAMM toolkit in publications, this website can be cited as::

  Dr. Scott W. Sides, Dr. Travis W. Kemper, Dr. Ross E. Larsen and Dr. Peter Graf. "STREAMM (Simulation Toolkit for
  Renewable Energy and Advanced Materials Modeling)," National Renewable Energy Lab, 21 Sept. 2015. <http://github.com/NREL/streamm-tools>.

