.. _resource_example:
  
resource_example
========================
 

In this example, we will use the resource object to keep track of file
locations

.. code:: python

    from pprint import pprint

.. code:: python

    import streamm
    import os 

.. code:: python

    local = streamm.Resource('local')

Set the ``templates`` directory

.. code:: python

    local.dir['templates'] =  os.path.join(local.dir['home'],'..','templates','')

This reads in the current directory as the root directory for a project

.. code:: python

    pprint(local.dir)


.. parsed-literal::

    {u'home': '/Users/my_username/Development/streamm-tools/examples',
     u'launch': u'/Users/my_username/Development/streamm-tools/examples/scratch',
     u'materials': u'/Users/my_username/Development/streamm-tools/examples/materials',
     u'scratch': u'/Users/my_username/Development/streamm-tools/examples/scratch',
     u'scripts': u'/Users/my_username/Development/streamm-tools/examples/scripts',
     u'storage': u'/Users/my_username/Development/streamm-tools/examples/storage',
     u'templates': '/Users/my_username/Development/streamm-tools/examples/../templates/'}


These directories can be created using the make\_dir() function

.. code:: python

    local.make_dir()

Then the ``Calculation`` object can use the location of these
directories and files within them to copy files to the correct locations

.. code:: python

    local_json = local.export_json()

We can also setup a resource we can access using ssh calls

.. code:: python

    remote = streamm.Resource('remote')

Set the type to ``ssh``, this will trigger some if statements in the
calculation object to scp calculation files to the external resource.

.. code:: python

    remote.meta['type'] = "ssh"

Enter your username and the address of the resource

.. code:: python

    remote.ssh['username'] = 'my_username'
    remote.ssh['address'] = 'system_address'

Then add the direcotry structure to the ``dir`` dictionary of the
resource

.. code:: python

    remote.dir['storage'] = '/storage/%s'%(remote.ssh['username'])
    remote.dir['scratch'] = '/scratch/%s'%(remote.ssh['username'])
    remote.dir['home'] = local.dir['home']
    remote.dir['launch'] = local.dir['launch']
    remote.dir['templates'] = local.dir['templates']

If you are running on a remote resource you can decide whether to set
the properties['exe\_command'] to ``qsub`` to submit the calculation to
a queuing system or ``./`` to run the calculation on the same node the
script is running on.

.. code:: python

    remote.properties['exe_command'] = 'qsub '

.. code:: python

    ssh_json = remote.export_json()

If you are running a streamm script on a remote resource you will want
to set the type to ``local``

.. code:: python

    remote.meta['type'] = "local"

Also, you will want to set the launch directory to the scratch location

.. code:: python

    remote.dir['launch'] = remote.dir['scratch']

This for the other examples like ``P3HT_ET`` were input files are copied
to the ``launch`` directory rather than the directly to the scratch
directory in case the example is accessing a remote resource via ``ssh``

Meh, whatever.
