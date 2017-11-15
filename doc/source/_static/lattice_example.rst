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


.. parsed-literal::

    [ 25.   0.   0.]


If we want a tuple of the vector and the magnitude we can use

.. code:: python

    dr_ij,mag_dr_ij =  lat.delta_pos(pos_i,pos_j)
    print dr_ij,mag_dr_ij


.. parsed-literal::

    [ 25.   0.   0.] 25.0


We can also turn pbcs off and calculate the distance

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


.. parsed-literal::

    [[ 10.39230485   0.          -6.        ]
     [  3.91349423   5.71704145   4.        ]
     [  0.           0.          15.        ]]
    [12.0, 8.0, 15.0, 60.0, 120.0, 80.0]
    ang


Let’s turn pbcs’s back on and calculate the distance

.. code:: python

    lat.pbcs = [True,True,True]

.. code:: python

    print pos_i,pos_j


.. parsed-literal::

    [25.0, 25.0, 25.0] [-50.0, 25.0, 25.0]


.. code:: python

    dr_ij,mag_dr_ij =  lat.delta_pos(pos_i,pos_j)
    print dr_ij,mag_dr_ij


.. parsed-literal::

    [-2.25386608  0.          3.        ] 3.75232092392


Change the units to ``nm``

.. code:: python

    lat.update_units({'length':'nm'})

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length']


.. parsed-literal::

    [[ 1.03923048  0.         -0.6       ]
     [ 0.39134942  0.57170414  0.4       ]
     [ 0.          0.          1.5       ]]
    [1.2, 0.79999999999999993, 1.4999999999999998, 60.0, 120.0, 80.0]
    nm


If you need your angles in radians

.. code:: python

    lat.update_units({'angle':'radian'})

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length'],lat.unit_conf['angle']


.. parsed-literal::

    [[ 1.03923048  0.         -0.6       ]
     [ 0.39134942  0.57170414  0.4       ]
     [ 0.          0.          1.5       ]]
    [1.2, 0.79999999999999993, 1.4999999999999998, 1.0471975511965976, 2.0943951023931953, 1.3962634015954636]
    nm radian


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


.. parsed-literal::

    Reading lat_ex_lat.json


Handy for saving or exporting to javascript

.. code:: python

    print lat.matrix
    print lat.constants
    print lat.unit_conf['length'],lat.unit_conf['angle']


.. parsed-literal::

    [[ 1.03923048  0.         -0.6       ]
     [ 0.39134942  0.57170414  0.4       ]
     [ 0.          0.          1.5       ]]
    [ 1.2         0.8         1.5         1.04719755  2.0943951   1.3962634 ]
    nm radian


Cool, aye!
