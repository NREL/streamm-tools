.. _lattice_example:
  
lattice_example
========================
 

.. code:: python

    from pprint import pprint

.. code:: python

    import logging
    logging.basicConfig(filename='lattice_example.log',level=logging.DEBUG)

The ``lattice`` object keeps track of the lattice used in ``Structure``
object

.. code:: python

    import pymatgen_core.core.lattice as lattice

.. code:: python

    lat = lattice.Lattice()

.. code:: python

    print lat
    print lat.unit_conf['length']

The default size of the lattice is 100.0 angstroms

The lattice also has lattice constants

.. code:: python

    print lat.constants
    print lat.unit_conf['length']

Which are returned as [a,b,c,alpha,beta,gamma]

We can calculate the distance between two points in the lattice

Let’s turn on periodic boundary conditions

.. code:: python

    lat.pbcs = [True,True,True]

.. code:: python

    pos_i = [25.0,25.0,25.0]
    pos_j = [-50.0,25.0,25.0]

.. code:: python

    dr_ij = lat.d_pos(pos_i,pos_j)
    print dr_ij

If we want a tuple of the vector and the magnitude we can use

.. code:: python

    dr_ij,mag_dr_ij =  lat.delta_pos(pos_i,pos_j)
    print dr_ij,mag_dr_ij

We can also turn pbcs off and calculate the distance

.. code:: python

    lat.pbcs = [False,False,False]

.. code:: python

    print lat.delta_pos(pos_i,pos_j)

The size of the lattice can be changed using the ``matrix`` or the
``constants`` ``setter``

.. code:: python

    lat.matrix = [ 12,0,0,0,12,0,0,0,12 ]

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length']

To set to a triclinic lattice

.. code:: python

    lat.constants = [ 12,8,15,60.0,120.0,80.0 ]

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length']

Let’s turn pbcs’s back on and calculate the distance

.. code:: python

    lat.pbcs = [True,True,True]

.. code:: python

    print pos_i,pos_j

.. code:: python

    dr_ij,mag_dr_ij =  lat.delta_pos(pos_i,pos_j)
    print dr_ij,mag_dr_ij

Change the units to ``nm``

.. code:: python

    lat.update_units({'length':'nm'})

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length']

If you need your angles in radians

.. code:: python

    lat.update_units({'angle':'radian'})

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length'],lat.unit_conf['angle']

We can export the lattice object as json object and dump it into a file

.. code:: python

    lat_json = lat.export_json('lat_ex',write_file=True)

Delete the lattice object

.. code:: python

    del lat

Create a new blank object

.. code:: python

    lat = lattice.Lattice()

And read in the file to get the properties of the lattice back

.. code:: python

    lat.import_json('lat_ex',read_file=True)

Handy for saving or exporting to javascript

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length'],lat.unit_conf['angle']

Cool, aye!
