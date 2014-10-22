.. _framework_design:

.. index:: Particle, Bond, Angle, Dihedral, Simulation, Structure, Class


*********************************************
STREAMM design
*********************************************


StructureContainer Class
============================

The StructureContainer class is a data structure for describing a
collection of discrete particles. These particles have unique integer
IDs and spatial positions and can have user specified attributes such
as mass, charge and type. These particles can have 2-body (bond),
3-body (angle) and 4-body (dihedral) interactions, which are
implemented by the python classes

- Particle
- Bond
- Angle
- Dihedral

Each of these classes has an associated container. For example,
multiple Particle objects are held in an instance of a
ParticleContainer class. These containers are essentially modified
dictionaries that map a unique integer ID to an object. Finally, a
StructureContainer object is a container for the Particle, Bond,
Angle and Dihedral containers

.. figure:: figures/strucC.png
   :align: center
   :scale: 40%

The contained ParticleContainer object has a special role... this object is a
dictionary that tracks the globalIDs for all particles in this Structure.
Global IDs for a structure are unique. All IDs referred to by the
subcontainers are consistent.

..  .. inheritance-diagram:: structureContainer.StructureContainer


Simulation Class
============================


Algorithms
============================
