
.. code:: python

    %load_ext autoreload
    %autoreload 2

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


.. parsed-literal::

    100.000000 0.000000 0.000000
    0.000000 100.000000 0.000000
    0.000000 0.000000 100.000000
    ang


The default size of the lattice is 100.0 angstroms

The lattice also has lattice constants

.. code:: python

    print lat.constants
    print lat.unit_conf['length']


.. parsed-literal::

    [ 100.  100.  100.   90.   90.   90.]
    ang


Which are returned as [a,b,c,alpha,beta,gamma]

We can calculate the distance between to points in the lattice

Let's turn Periodic boundary conditions

.. code:: python

    lat.pbcs = [True,True,True]

.. code:: python

    pos_i = [25.0,25.0,25.0]
    pos_j = [-50.0,25.0,25.0]

.. code:: python

    dr_ij = lat.d_pos(pos_i,pos_j)
    print dr_ij


.. parsed-literal::

    [ 25.   0.   0.]


If we want a tuple of the vector and the magnitude we can use

.. code:: python

    dr_ij,mag_dr_ij =  lat.delta_pos(pos_i,pos_j)
    print dr_ij,mag_dr_ij


.. parsed-literal::

    [ 25.   0.   0.] 25.0


We can also turn pbcs off and calculated the distance

.. code:: python

    lat.pbcs = [False,False,False]

.. code:: python

    print lat.delta_pos(pos_i,pos_j)


.. parsed-literal::

    (array([-75.,   0.,   0.]), 75.0)


The size of the lattice can be changed using the ``matrix`` or the
``constants`` ``setter``

.. code:: python

    lat.matrix = [ 12,0,0,0,12,0,0,0,12 ]

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length']


.. parsed-literal::

    [[ 12.   0.   0.]
     [  0.  12.   0.]
     [  0.   0.  12.]]
    [ 12.  12.  12.  90.  90.  90.]
    ang


To set to a triclinic lattice

.. code:: python

    lat.constants = [ 12,8,15,60.0,120.0,80.0 ]

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length']

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

Cool, aye!

.. code:: python

    lat = 
