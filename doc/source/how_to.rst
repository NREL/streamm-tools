.. _how_to:

How To
******

.. _read_xyz:

Read `.xyz`
===========

The easiest way to create an organic structure to read into streamm is to use a molecular viewer, such as `Avogadro <https://avogadro.cc/>`_.
You can generate a :class:`streamm.Buildingblock <streamm.structures.buildingblock.Buildingblock>` structure by inserting and connecting atoms or importing fragments.

.. Note::

    If you use Avogadro there are fragments available in `Build>Insert>Fragment`


Once you have created your structure export it in the `.xyz <https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html>`_  format,
either place the `.xyz` file in the working directory of your project or navigate to its location using the `os.chdir()` command.
The structure can be read in with the :func:`read_xyz() <streamm.structures.structure.Structure.read_xyz>` function
by creating an empty object and either setting the `tag` of the object to the same name as prefix of the `.xyz` file

.. code:: python

    mol101 = streamm.Buildingblock('mol101')
    mol101.read_xyz()

or using the keyword `xyz_file` can be used to set the file name.


.. code:: python

    mol101.read_xyz(xyz_file='mol101.xyz')


.. _change_units:

Change units
============

You can change the units of object properties with the :func:`update_units() <pymatgen_core.core.units.ObjectUnits.update_units>` function.
This function takes a dictionary with the unit type as the key and the new unit to be used as the value.
For example, the units of a Buildingblock can be changed to `nm` by running

.. code:: python

    mol101.update_units({'length':'nm'})
    
    
    
    
    

