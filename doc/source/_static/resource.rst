.. _resource:

resource
========


.. code:: python

    from pprint import pprint

In this example, we will use the resource object to keep track of file
locations

.. code:: python

    import streamm

.. code:: python

    res_local = streamm.Resource('local')

This reads in the current directory as the root directory for a project

.. code:: python

    pprint(res_local.dir)


.. parsed-literal::

    {u'home': '/Users/tkemper/Development/STREAMM/streamm-tools/examples',
     u'launch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/scratch',
     u'materials': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/materials',
     u'scratch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/scratch',
     u'scripts': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/scripts',
     u'storage': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/storage',
     u'templates': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/templates'}


These directories can be created using the make\_dir() function

.. code:: python

    res_local.make_dir()

Then the ``Calculation`` object can use the location of these
directories and files within them to copy files to the correct locations

Meh, whatever.
