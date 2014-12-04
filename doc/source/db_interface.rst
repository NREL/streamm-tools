.. _db_interface:

.. index:: query, SQL, database, OPV

*********************************************
Interface to NREL database
*********************************************

One of the projects at NREL that uses the STREAMM toolkit is the
Organic Photovoltaics (OPV) modeling effort. High-throughput density
functional theory (DFT) calculations are used to screen candidate
molecules so as to let synthetic chemists focus on only the most
promising materials for energy applications. The simulation results
from ~10,000's of calculations are stored in the NREL OPV database
The STREAMM toolkit is used not only to generate the simulation files
needed to perform the DFT calculations, but also to access and
effectively mine the massive amounts of data generated on the
Peregrine compute cluster. The python modules and scripts that enable
efficient access to the OPV database is described below.


Description of OPV data
=======================================

The SQL database is organized into two major sections. The raw data
from DFT calculations is in the *vw_oligomer* section. The
extrapolated results that use the oligomer data is in the
*vw_structure* section. The various keywords are shown in the tables
below:

.. toctree::
   :maxdepth: 1

   table_oligomer.rst
   table_structure.rst


Desciption of opvSQL python module
======================================================

.. automodule:: opvSQL
   :members:
   :show-inheritance:


Examples of using STREAMM classes with database
======================================================

simulationGaussian example
