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

.. code:: python

    value = 10.0 
    old_unit = 'mile'
    new_unit = 'ang'
    print "Conversion of {} to {} is {}".format(old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))

To convert units of energy

.. code:: python

    Unit_instance = units.partial(units.FloatWithUnit, unit_type='energy')

.. code:: python

    value = 37.500000
    old_unit = 'kCalmol'
    new_unit = 'kJmol'
    print "Conversion of {} {} to {} is {}".format(value,old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))

.. code:: python

    value = 1.0 
    old_unit = 'eV'
    for new_unit in ['Ha','J','wavenumber','kJmol','kCalmol']:
        print "Conversion of {} {} to {} is {}".format(value,old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))

We have special units for bond stretching of energy/length^2

.. code:: python

    Unit_instance = units.partial(units.FloatWithUnit, unit_type='harm_bond_coeff')

.. code:: python

    value = 367.000000
    old_unit = 'kCalmolsqang'
    new_unit = 'kJmolsqnm'
    print "Conversion of {} {} to {} is {}".format(value,old_unit,new_unit,Unit_instance(value,old_unit).to(new_unit))

Well, that’s handy!
