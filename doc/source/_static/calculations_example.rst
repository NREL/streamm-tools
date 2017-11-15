.. _calculations_example:
  
calculations_example
===============
 

.. code:: ipython2

    import os 
    from pprint import pprint
    from pathlib2 import Path

The aim of the streamm package is to design a project with numerous
calculations, run them on local and remote resources and collect the
output for analysis, thus facilitating high-throughput computational
material design.

To accomplish this, the directory structure is contained within a
resource as a dictionary. Resources, structures, and forcefields are
contained within a calculation object. Sets of calculations are
contained within a project

So let’s first create a resource object that we will use to set the
directory locations of all the subsequent calculation objects

.. code:: ipython2

    import streamm

We need some output from another example so check that it has been
created first

.. code:: ipython2

    need_files = ['ethane_struc.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run buildingblocks_example.ipynb")
            os.system("jupyter nbconvert --to python  buildingblocks_example.ipynb")
            os.system("python buildingblocks_example.py")

.. code:: ipython2

    need_files = ['oplsaa_param.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run forcefields_example.ipynb")
            os.system("jupyter nbconvert --to python  forcefields_example.ipynb")
            os.system("python forcefields_example.py")

.. code:: ipython2

    import streamm.calculations.resource as resource  

.. code:: ipython2

    res_i = resource.Resource('local')

This sets the current working directory as the root/home directory by
default

.. code:: ipython2

    pprint(res_i.dir)

.. code:: ipython2

    EXAMPLE_DIR = res_i.dir['home']

However, we want to use structures from our previous structures and
forcefields examples, so let’s set the materials directory to examples/

.. code:: ipython2

    res_i.dir['materials'] = res_i.dir['home']

To write out input files we will use the templates provided in the
streamm package

Set the template dir dictionary entry to the location of templates
directory

.. code:: ipython2

    res_i.dir['templates'] =  os.path.join(EXAMPLE_DIR,'..','templates','')

.. code:: ipython2

    print res_i.dir['templates']

This also contains the properties dictionary, which can be used to write
.pbs scripts on clusters

.. code:: ipython2

    pprint(res_i.properties)

By default the resource type is ‘local’; however, setting type to ‘ssh’
will invoke an scp command when copying files

Okay create the directories we need for our calculation

.. code:: ipython2

    res_i.make_dir()

Now we should have materials, scratch, scripts, storage and templates
directories

We can create a gaussian calculation

.. code:: ipython2

    import streamm.calculations.gaussian as gaussian  

.. code:: ipython2

    calc_i = gaussian.Gaussian('ethane_HF')

Set the resource and all the directories

.. code:: ipython2

    calc_i.set_resource(res_i)

.. code:: ipython2

    pprint(calc_i.dir)

Make the calculation directories

.. code:: ipython2

    calc_i.make_dir()

Let’s assign a structure to this calculation

First copy the .xyz file from the materials directory to our scratch
directory using the cp_file() function.

.. code:: ipython2

    os.chdir(calc_i.dir['home'])

This takes an type and key to set the calc_i.files[type][key] dictionary

.. code:: ipython2

    file_type = 'input'
    file_key = 'xyz'
    file_name = "ethane_struc.json"
    from_dirkey = 'materials'
    to_dirkey = 'scratch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Generally, the materials directory is thought to contain a repository of
material files, and local versions in the scratch directory should be
made in case modifications are necessary

Change to the scratch directory

.. code:: ipython2

    pprint(calc_i.dir['scratch'])

.. code:: ipython2

    os.chdir(calc_i.dir['scratch'])

Read in methane ``.json`` file from the structures example

.. code:: ipython2

    calc_i.strucC.tag = 'ethane'
    calc_i.strucC.import_json(read_file=True)

.. code:: ipython2

    print(calc_i.strucC.print_properties())

Now that we have a structure and parameters for each interaction we can
create an input file for a simulation

Get the bash run script for Gaussian. By setting the file_key to run,
this will be the script that executed when the run() function is called

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "gaussian.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)


Get the .com template

.. code:: ipython2

    
    file_type = 'templates'
    file_key = 'com'
    file_name = "gaussian.com"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)


Make sure we are in the scratch directory

.. code:: ipython2

    pprint(os.getcwd())

Load the template files into memory

.. code:: ipython2

    calc_i.load_str('templates','com')        
    calc_i.load_str('templates','run')

Set the properties strings in the template files to have the values we
want

.. code:: ipython2

    calc_i.properties['commands'] = 'HF/3-21G SP'
    calc_i.properties['charge'] = 0
    calc_i.properties['spin_mult'] = 1
    calc_i.properties['coord'] = calc_i.strucC.write_coord()

Replace the strings in the files[‘input’][‘com’]

.. code:: ipython2

    calc_i.replacewrite_prop('com','input','com','%s.com'%(calc_i.tag))

Add the name of the .com file to the properties, and replace the strings
in the files[‘input’][‘run’]

.. code:: ipython2

    calc_i.properties['input_com'] = calc_i.files['input']['com']
    calc_i.replacewrite_prop('run','scripts','run','%s.sh'%(calc_i.tag))

Save a .json file in the home directory

.. code:: ipython2

    os.chdir(calc_i.dir['home'])
    calc_json = calc_i.export_json()

Go to scratch directory and see if there is a completed output file for
the calculation

.. code:: ipython2

    os.chdir(calc_i.dir['scratch'])
    calc_i.check()

Check the status

.. code:: ipython2

    pprint("Calculation:{} has status:{}".format(calc_i.tag,calc_i.meta['status']))

If you have gaussian installed on your machine and g09 in your PATH you
can run the bash script

.. code:: ipython2

    calc_i.run()

You can read in the data from the log file

.. code:: ipython2

    calc_i.add_file('output','log','{}.log'.format(calc_i.strucC.tag))

.. code:: ipython2

    calc_i.check()
    if(calc_i.meta['status'] == 'finished' ):
        calc_i.analysis()

Then compress the results and copy them to storage

.. code:: ipython2

    calc_i.store()

Next we can follow a similar procedure to run a LAMMPS MD simulation

.. code:: ipython2

    import streamm.calculations.lammps as lammps  

.. code:: ipython2

    calc_j = lammps.LAMMPS('ethane_lmp')

Set the resource

.. code:: ipython2

    calc_j.set_resource(res_i)

Make directories

.. code:: ipython2

    calc_j.make_dir()

.. code:: ipython2

    pprint(calc_j.dir)

This takes an type and key to set the calc_i.files[type][key] dictionary

.. code:: ipython2

    file_type = 'input'
    file_key = 'xyz'
    file_name = "ethane_struc.json"
    from_dirkey = 'materials'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    os.chdir(calc_j.dir['scratch'])

Read in the ethane .json file from the structures example

.. code:: ipython2

    calc_j.strucC.tag = 'ethane'
    calc_j.strucC.import_json(read_file=True)

.. code:: ipython2

    print(calc_j.strucC.print_properties())

Set the forcefield particletypes

.. code:: ipython2

    for pkey,p in calc_j.strucC.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC'

Copy the forcefield parameter .json file to scratch and read it in

.. code:: ipython2

    file_type = 'input'
    file_key = 'param'
    file_name = "oplsaa_param.json"
    from_dirkey = 'materials'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    import streamm.forcefields.parameters as parameters 

.. code:: ipython2

    calc_j.paramC = parameters.Parameters('oplsaa')

.. code:: ipython2

    calc_j.paramC.import_json()

.. code:: ipython2

    print calc_j.paramC

.. code:: ipython2

    for ptkey,pt in calc_j.paramC.particletypes.iteritems():
        print ptkey,pt,pt.unit_conf['energy'],pt.unit_conf['length']

.. code:: ipython2

    for btkey,bt in calc_j.paramC.bondtypes.iteritems():
        print btkey,bt,bt.unit_conf['harm_bond_coeff'],pt.unit_conf['length']

.. code:: ipython2

    for atkey,at in calc_j.paramC.angletypes.iteritems():
        print atkey,at,at.unit_conf['energy'],at.unit_conf['length']

Use the set_ffparam() function to iterate through the structure
container and set parameters based on ``paramkeys``

.. code:: ipython2

    calc_j.set_ffparam()

Now we have a structure that has forcefield parameters for each
particle, bond and bond angle

Let’s get the input file template

.. code:: ipython2

    file_type = 'templates'
    file_key = 'in'
    file_name = "lammps_sp.in"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Bash run file

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "lammps.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    calc_j.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Got to scratch dir

.. code:: ipython2

    os.chdir(calc_j.dir['scratch'])

Read in template files

.. code:: ipython2

    calc_j.load_str('templates','in')
    calc_j.load_str('templates','run')

Write LAMMPS data file

.. code:: ipython2

    calc_j.write_data()

Replace properties strings in template and write template

.. code:: ipython2

    calc_j.replacewrite_prop('in','input','in','%s.in'%(calc_j.tag))

Set .in file in properties and write run script

.. code:: ipython2

    calc_j.properties['input_in'] = calc_j.files['input']['in']
    calc_j.replacewrite_prop('run','scripts','run','%s.sh'%(calc_j.tag))

Save a .json file in the home directory

.. code:: ipython2

    os.chdir(calc_j.dir['home'])
    calc_json = calc_j.export_json()

Go to scratch directory and see if there is a completed output file for
the calculation

.. code:: ipython2

    os.chdir(calc_j.dir['scratch'])
    calc_j.check()

.. code:: ipython2

    pprint("Calculation:{} has status:{}".format(calc_j.tag,calc_j.meta['status']))

So now we have two calculations, let’s put them in a project so we can
operate on them both at the same time

.. code:: ipython2

    import streamm.calculations.project as project  

.. code:: ipython2

    import copy

.. code:: ipython2

    proj_i = streamm.Project('example_proj')

.. code:: ipython2

    proj_i.add_calc(calc_i,deepcopy=True)
    proj_i.add_calc(calc_j,deepcopy=True)

Now we can check the status of each calculation with a single command

.. code:: ipython2

    proj_i.check()

We can run each simulation

.. code:: ipython2

    proj_i.run()

We can tar up the results and copy the tar files to a storage location

.. code:: ipython2

    proj_i.store()

And dump the details of the project to a json file

.. code:: ipython2

    os.chdir(calc_i.dir['home'])
    proj_i.export_json()

.. code:: ipython2

    del proj_i

.. code:: ipython2

    proj_i = streamm.Project('example_proj')

.. code:: ipython2

    proj_i.import_json()

.. code:: ipython2

    proj_i.check()

Neat-O!
