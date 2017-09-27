
.. code:: python

    from __future__ import division, unicode_literals

.. code:: python

    import os 
    from pprint import pprint

The aim of the streamm package is to design a project with numerous
calculations run them on local and remote resources and collect the
output for analysis, thus facilitating high-throughput computational
material design.

To accomplish, the directory structure is contained within a resource as
a dictionary. Resources, structures, and forcefields are contained
within a calculation object. Sets of calculations are contained within a
project

So let's first create a resource object that we will use to set the
directory locations of all the subsequent calculation objects

.. code:: python

    import streamm

.. code:: python

    import streamm.calculations.resource as resource  

.. code:: python

    res_i = resource.Resource('local')

This sets the current working directory as the root/home directory by
default

.. code:: python

    pprint(res_i.dir)


.. parsed-literal::

    {u'home': '/Users/tkemper/Development/STREAMM/streamm-tools/examples',
     u'launch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/scratch',
     u'materials': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/materials',
     u'scratch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/scratch',
     u'scripts': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/scripts',
     u'storage': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/storage',
     u'templates': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/templates'}


Let's create a new home directory called 'example\_proj'

.. code:: python

    EXAMPLE_DIR = res_i.dir['home']

.. code:: python

    new_root = os.path.join(EXAMPLE_DIR,'example_proj')

.. code:: python

    res_i.set_home(new_root)

.. code:: python

    pprint(res_i.dir)


.. parsed-literal::

    {u'home': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj',
     u'launch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch',
     u'materials': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/materials',
     u'scratch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch',
     u'scripts': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scripts',
     u'storage': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/storage',
     u'templates': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/templates'}


However, we want to use structures from our previous structures and
forcefields examples, so let's set the materials directory to examples/

.. code:: python

    res_i.dir['materials'] = EXAMPLE_DIR

To write out input files we will use the templates provided in the
streamm package

Set the template dir dictionary entry to the location of templates
directory

.. code:: python

    res_i.dir['templates'] =  os.path.join(EXAMPLE_DIR,'..','templates','')

.. code:: python

    print res_i.dir['templates']


.. parsed-literal::

    /Users/tkemper/Development/STREAMM/streamm-tools/examples/../templates/


This also contains the properties dictionary, which can be used to write
.pbs scripts on clusters

.. code:: python

    pprint(res_i.properties)


.. parsed-literal::

    {u'allocation': u'',
     u'exe_command': u'./',
     u'feature': u'24core',
     u'nodes': 1,
     u'nproc': 1,
     u'pmem': 1500,
     u'ppn': 1,
     u'queue': u'batch',
     u'walltime': 24}


By default the resource type is 'local'; however, setting type to 'ssh'
will invoke an scp command when copying files

Okay create the directories we need for our calculation

.. code:: python

    res_i.make_dir()

Now we should have a directory 'example\_proj/' with materials, scratch,
scripts, storage and templates directories

We can create a gaussian calculation

.. code:: python

    import streamm.calculations.gaussian as gaussian  

.. code:: python

    calc_i = gaussian.Gaussian('methane_HF')

Set the resource and all the directories

.. code:: python

    calc_i.set_resource(res_i)

.. code:: python

    pprint(calc_i.dir)


.. parsed-literal::

    {u'home': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj',
     u'launch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_HF/',
     u'materials': '/Users/tkemper/Development/STREAMM/streamm-tools/examples',
     u'scratch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_HF/',
     u'scripts': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scripts',
     u'storage': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/storage/methane_HF/',
     u'templates': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/../templates/'}


Make the calculation directories

.. code:: python

    calc_i.make_dir()

Let's assign a structure to this calculation

First copy the .xyz file from the materials directory to our scratch
directory using the cp\_file() function.

This takes an type and key to set the calc\_i.files[type][key]
dictionary

.. code:: python

    file_type = 'input'
    file_key = 'xyz'
    file_name = "methane.xyz"
    from_dirkey = 'materials'
    to_dirkey = 'scratch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Generally, the materials directory is thought to contain a repository of
material files, and local versions in the scratch directory should be
made in case modifications are necessary

Change to the scratch directory

.. code:: python

    pprint(calc_i.dir['scratch'])


.. parsed-literal::

    u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_HF/'


.. code:: python

    os.chdir(calc_i.dir['scratch'])

Read in methane .xyz file from the structures example

.. code:: python

    calc_i.strucC.read_xyz('methane.xyz')

.. code:: python

    print calc_i.strucC.n_particles


.. parsed-literal::

    5


Now that we have a structure and parameters for each interaction we can
create an input file for a simulation

Get the bash run script for Gaussian. By setting the file\_key to run,
this will be the script that executed when the run() function is called

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "gaussian.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)


Get the .com template

.. code:: python

    
    file_type = 'templates'
    file_key = 'com'
    file_name = "gaussian.com"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)


Make sure we are in the scratch directory

.. code:: python

    pprint(os.getcwd())


.. parsed-literal::

    '/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_HF'


Load the template files into memory

.. code:: python

    calc_i.load_str('templates','com')        
    calc_i.load_str('templates','run')

Set the properties strings in the template files to have the values we
want

.. code:: python

    calc_i.properties['commands'] = 'HF/3-21G SP'
    calc_i.properties['charge'] = 0
    calc_i.properties['spin_mult'] = 1
    calc_i.properties['coord'] = calc_i.strucC.write_coord()

Replace the strings in the files['input']['com']

.. code:: python

    calc_i.replacewrite_prop('com','input','com','%s.com'%(calc_i.tag))

Add the name of the .com file to the properties, and replace the strings
in the files['input']['run']

.. code:: python

    calc_i.properties['input_com'] = calc_i.files['input']['com']
    calc_i.replacewrite_prop('run','scripts','run','%s.sh'%(calc_i.tag))

Save a .json file in the home directory

.. code:: python

    os.chdir(calc_i.dir['home'])
    calc_i.dump_json()

Go to scratch directory and see if there is a completed output file for
the calculation

.. code:: python

    os.chdir(calc_i.dir['scratch'])
    calc_i.check()

Check the status

.. code:: python

    pprint("Calculation:{} has status:{}".format(calc_i.tag,calc_i.meta['status']))


.. parsed-literal::

    u'Calculation:methane_HF has status:written'


If you have gaussian installed on your machine and g09 in your PATH you
can run the bash script

.. code:: python

    calc_i.run()

You can read in the data from the log file

.. code:: python

    calc_i.add_file('output','log','{}.log'.format(calc_i.strucC.tag))

.. code:: python

    calc_i.check()
    if(calc_i.meta['status'] == 'finished' ):
        calc_i.analysis()

Then compress the results and copy them to storage

.. code:: python

    calc_i.store()

Next we can follow a similar procedure to run a LAMMPS MD simulation

.. code:: python

    import streamm.calculations.lammps as lammps  

.. code:: python

    calc_j = lammps.LAMMPS('methane_lmp')

Set the resource

.. code:: python

    calc_j.set_resource(res_i)

Make directories

.. code:: python

    calc_j.make_dir()

.. code:: python

    pprint(calc_j.dir)


.. parsed-literal::

    {u'home': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj',
     u'launch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_lmp/',
     u'materials': '/Users/tkemper/Development/STREAMM/streamm-tools/examples',
     u'scratch': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_lmp/',
     u'scripts': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scripts',
     u'storage': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/storage/methane_lmp/',
     u'templates': u'/Users/tkemper/Development/STREAMM/streamm-tools/examples/../templates/'}


This takes an type and key to set the calc\_i.files[type][key]
dictionary

.. code:: python

    file_type = 'input'
    file_key = 'xyz'
    file_name = "methane.xyz"
    from_dirkey = 'materials'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    os.chdir(calc_j.dir['scratch'])

Read in methane .xyz file from the structures example

.. code:: python

    calc_j.strucC.read_xyz('methane.xyz')

.. code:: python

    print calc_j.strucC.n_particles


.. parsed-literal::

    5


Set the forcefield particletypes

.. code:: python

    for pkey,p in calc_j.strucC.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC'

Set neighbor list

.. code:: python

    calc_j.strucC.bonded_nblist = calc_j.strucC.guess_nblist(0,radii_buffer=1.25)

Find bonds and bond angles based on neighbor list

.. code:: python

    calc_j.strucC.bonded_bonds()
    calc_j.strucC.bonded_angles()

Copy the pickled forcefield parameter file to scratch and read it in

.. code:: python

    file_type = 'input'
    file_key = 'param'
    file_name = "oplsaa.pkl"
    from_dirkey = 'materials'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    import streamm.forcefields.parameters as parameters 

.. code:: python

    calc_j.paramC = parameters.read_pickle('oplsaa')

.. code:: python

    print calc_j.paramC


.. parsed-literal::

    
        Parameters 
          LJ parameters 2 
          Bond parameters 2 
          Angle parameters 2 
          Dihedral parameters 1 
          Imporper Dihedral parameters 0 
    


.. code:: python

    for ptkey,pt in calc_j.paramC.particletypes.iteritems():
        print ptkey,pt,pt.unit_conf['energy'],pt.unit_conf['length']


.. parsed-literal::

    0  CT epsilon:0.066 sigma:3.5 kCalmol ang
    1  HC epsilon:0.03 sigma:2.5 kCalmol ang


.. code:: python

    for btkey,bt in calc_j.paramC.bondtypes.iteritems():
        print btkey,bt,bt.unit_conf['harm_bond_coeff'],pt.unit_conf['length']


.. parsed-literal::

    0  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   kCalmolsqang ang
    1  bond  CT - CT type harmonic 
      harmonic r_0 = 1.530000 K = 268.000000 lammps index 0  gromacs index 0   kCalmolsqang ang


.. code:: python

    for atkey,at in calc_j.paramC.angletypes.iteritems():
        print atkey,at,at.unit_conf['energy'],at.unit_conf['length']


.. parsed-literal::

    0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   kCalmol ang
    1  angle  HC - CT - CT type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   kCalmol ang


Use the set\_ffparam() function to iterate through the structure
container and set parameters based on ffkeys

.. code:: python

    calc_j.set_ffparam()

Now we have a structure that has forcefield parameters for each
particle, bond and bond angle

Let's get the input file template

.. code:: python

    file_type = 'templates'
    file_key = 'in'
    file_name = "lammps_sp.in"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Bash run file

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "lammps.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Got to scratch dir

.. code:: python

    os.chdir(calc_j.dir['scratch'])

Read in template files

.. code:: python

    calc_j.load_str('templates','in')
    calc_j.load_str('templates','run')

Write LAMMPS data file

.. code:: python

    calc_j.write_data()

Replace properties strings in template and write template

.. code:: python

    calc_j.replacewrite_prop('in','input','in','%s.in'%(calc_j.tag))

Set .in file in properties and write run script

.. code:: python

    calc_j.properties['input_in'] = calc_j.files['input']['in']
    calc_j.replacewrite_prop('run','scripts','run','%s.sh'%(calc_j.tag))

Save a .json file in the home directory

.. code:: python

    os.chdir(calc_j.dir['home'])
    calc_j.dump_json()

Go to scratch directory and see if there is a completed output file for
the calculation

.. code:: python

    os.chdir(calc_j.dir['scratch'])
    calc_j.check()

.. code:: python

    pprint("Calculation:{} has status:{}".format(calc_j.tag,calc_j.meta['status']))


.. parsed-literal::

    u'Calculation:methane_lmp has status:written'


So now we have two calculations, let's put them in a project so we can
operate on them both at the same time

.. code:: python

    import streamm.calculations.project as project  

.. code:: python

    import copy

.. code:: python

    proj_i = streamm.Project('example_proj')

.. code:: python

    proj_i.calculations[calc_i.tag] = copy.deepcopy(calc_i)
    proj_i.calculations[calc_j.tag] = copy.deepcopy(calc_j)

Now we can check the status of each calculation with a single command

.. code:: python

    proj_i.check()


.. parsed-literal::

    Calculation methane_lmp has status written
    Calculation methane_HF has status written


We can run each simulation

.. code:: python

    proj_i.run()


.. parsed-literal::

    /Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_lmp
    /Users/tkemper/Development/STREAMM/streamm-tools/examples/example_proj/scratch/methane_HF


We can tar up the results and copy the tar files to a storage location

.. code:: python

    proj_i.store()

Neat-O!
