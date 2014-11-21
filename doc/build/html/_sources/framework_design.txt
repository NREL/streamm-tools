.. _framework_design:

.. index:: Particle, Bond, Angle, Dihedral, Simulation, Structure, Class


*********************************************
STREAMM design
*********************************************

Object-oriented design is used throughout  STREAMM to enable ease of
use and the ability to update and extend its functionality. 


StructureContainer Class
============================

The StructureContainer class is a data structure for describing a
collection of discrete particles. These particles have unique integer
IDs and spatial positions and can have user specified attributes such
as mass, charge and type. These particles can have 2-body (bond),
3-body (angle) and 4-body (dihedral) interactions, which are
implemented by the python classes Particle, Bond, Angle and Dihedral.

Each of these classes has an associated container. For example,
multiple Particle objects are held in an instance of a
ParticleContainer class. These containers are essentially modified
dictionaries that map a unique integer ID to an object.
All IDs referred to by the subcontainers are consistent.
The embedded docstrings for the modules containing the classes above
and their related containers are:

- :mod:`particles`
- :mod:`bonds`
- :mod:`angles`
- :mod:`dihedrals`

Finally, a StructureContainer object is a container for the Particle, Bond,
Angle and Dihedral containers. This class is in the module
:mod:`structureContainer`

.. figure:: figures/strucC.png
   :align: center
   :scale: 40%



Simulation Class
============================

The embedded docstrings for the simulation modules are:

- :mod:`simulation`
- :mod:`simulationLAMMPS1`
- :mod:`simulationGaussian1`



Algorithms
============================
This may or may not be part of STREAMM design

.. toctree::
   :hidden:

   docstr_particles.rst
   docstr_bonds.rst
   docstr_angles.rst
   docstr_dihedrals.rst
   docstr_structureContainer.rst
   docstr_simulation.rst
   docstr_simulationLAMMPS1.rst
   docstr_simulationGaussian1.rst
