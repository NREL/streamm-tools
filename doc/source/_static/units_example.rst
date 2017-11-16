.. _units_example:
  
units_example
========================
 

.. code:: python

    from pprint import pprint

The units model is taken from the pymatgen code and modified to contain
more unit conversions for the streamm code.

.. code:: python

    import pymatgen_core.core.units as units 

Each object in streamm has a property ``_unit_conf`` that specifies the
units for the object.

The default ``unit_conf`` is in units

.. code:: python

    pprint(units.unit_conf)


.. parsed-literal::

    {u'amount': u'atom',
     u'angle': u'degree',
     u'capacitance': u'F',
     u'charge': u'e',
     u'conductance': u'S',
     u'current': u'A',
     u'density': u'amu_nm^3',
     u'electric_dipole_moment': u'D',
     u'emf': u'V',
     u'energy': u'Ha',
     u'force': u'GN',
     u'frequency': u'Hz',
     u'harm_bond_coeff': u'kCalmolsqang',
     u'intensity': u'cd',
     u'length': u'ang',
     u'magnetic_flux': u'Wb',
     u'mass': u'amu',
     u'memory': u'Kb',
     u'power': u'GW',
     u'pressure': u'KPa',
     u'resistance': u'ohm',
     u'temperature': u'K',
     u'time': u'ns',
     u'volume': u'nm^3'}


The units of most objects in streamm can be changed using an
``update_units`` function associated with that object

The units are then converted using a ``units_instance``

For example to convert units of length

.. code:: python

    Unit_instance = units.partial(units.FloatWithUnit, unit_type='length')

.. code:: python

    value = 10.0 
    old_unit = 'ang'
    new_unit = 'nm'
    print "Conversion of {} to {} is {}".format(old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))


.. parsed-literal::

    Conversion of ang to nm is 1.0 nm


.. code:: python

    value = 10.0 
    old_unit = 'mile'
    new_unit = 'ang'
    print "Conversion of {} to {} is {}".format(old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))


.. parsed-literal::

    Conversion of mile to ang is 1.609344e+14 ang


To convert units of energy

.. code:: python

    Unit_instance = units.partial(units.FloatWithUnit, unit_type='energy')

.. code:: python

    value = 37.500000
    old_unit = 'kCalmol'
    new_unit = 'kJmol'
    print "Conversion of {} {} to {} is {}".format(value,old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))


.. parsed-literal::

    Conversion of 37.5 kCalmol to kJmol is 156.9 kJmol


.. code:: python

    value = 1.0 
    old_unit = 'eV'
    for new_unit in ['Ha','J','wavenumber','kJmol','kCalmol']:
        print "Conversion of {} {} to {} is {}".format(value,old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))


.. parsed-literal::

    Conversion of 1.0 eV to Ha is 0.03674932248 Ha
    Conversion of 1.0 eV to J is 1.6021766208e-19 J
    Conversion of 1.0 eV to wavenumber is 8065.5440048 wavenumber
    Conversion of 1.0 eV to kJmol is 96.4853328825 kJmol
    Conversion of 1.0 eV to kCalmol is 23.0605480121 kCalmol


We have special units for bond stretching of energy/length^2

.. code:: python

    Unit_instance = units.partial(units.FloatWithUnit, unit_type='harm_bond_coeff')

.. code:: python

    value = 367.000000
    old_unit = 'kCalmolsqang'
    new_unit = 'kJmolsqnm'
    print "Conversion of {} {} to {} is {}".format(value,old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))


.. parsed-literal::

    Conversion of 367.0 kCalmolsqang to kJmolsqnm is 153552.8 kJmolsqnm


Well, that’s handy!
