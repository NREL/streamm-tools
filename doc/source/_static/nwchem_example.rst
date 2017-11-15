.. _nwchem_example:
  
nwchem_example
===============
 

.. code:: ipython2

    import os 
    from pprint import pprint

Check that output from other examples has been generated

.. code:: ipython2

    from pathlib2 import Path

.. code:: ipython2

    need_files = ['ethane_struc.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run buildingblocks_example.ipynb")
            os.system("jupyter nbconvert --to python  buildingblocks_example.ipynb")
            os.system("python buildingblocks_example.py")

.. code:: ipython2

    need_files = ['local_res.json','remote_res.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run resource_example.ipynb")
            os.system("jupyter nbconvert --to python  resource_example.ipynb")
            os.system("python resource_example.py")

.. code:: ipython2

    import streamm

Now let’s create project and resource to keep track of our work

.. code:: ipython2

    nwchem_example = streamm.Project('nwchem_example')
    res_local = streamm.Resource('local')

.. code:: ipython2

    res_local.make_dir()

Update relative location of templates directory

.. code:: ipython2

    res_local.dir['templates'] =  os.path.join(res_local.dir['home'],'..','templates','')

Make sure this is the location of the templates directory that comes
with the streamm git repository https://github.com/NREL/streamm-tools

.. code:: ipython2

    print res_local.dir['templates']

Create the local directories that will store our files

.. code:: ipython2

    nwchem_example.make_dir()

Tell the project about our directories

.. code:: ipython2

    nwchem_example.set_resource(res_local)

Read in the ethane structure created in the buildingblocks_example.ipynb
example

.. code:: ipython2

    ethane = streamm.Buildingblock('ethane')

.. code:: ipython2

    ethane.import_json()

.. code:: ipython2

    print ethane.print_properties()

Set the paramkeys so we can identify force field parameters later on

.. code:: ipython2

    for pkey,p in ethane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC'

.. code:: ipython2

    for pk,p in ethane.particles.iteritems():
        p.residue = 1
        p.resname = 'ETH'

Set ``rsite``\ ’s to hydrogens to be replaced during join

.. code:: ipython2

    ethane.particles[1].rsite = 'RH'
    ethane.particles[5].rsite = 'RH'

Run ``find_rsites()`` to populate ``func`` list

.. code:: ipython2

    ethane.find_rsites()

.. code:: ipython2

    print ethane.show_rsites()

.. code:: ipython2

    import copy

Create octane from ethane

Copy ethane to a new Buildingblock octane

.. code:: ipython2

    octane = copy.deepcopy(ethane)

.. code:: ipython2

    from streamm.structures.buildingblock import attach

Then attach 3 more ethanes to make an octane

.. code:: ipython2

    for i in range(3):
        octane = attach(octane,ethane,'RH',1,'RH',0)

.. code:: ipython2

    octane.tag = 'octane'

.. code:: ipython2

    print octane.n_particles

.. code:: ipython2

    octane.write_xyz()

Update the tag

.. code:: ipython2

    octane.tag = 'octane'

Rename the residue and resname for octane

.. code:: ipython2

    for pk,p in octane.particles.iteritems():
        p.residue = 2
        p.resname = "OCT"

.. code:: ipython2

    octane.write_xyz()

Create NWChem Calculation object

.. code:: ipython2

    nwchem_octane = streamm.NWChem('nw_octane_OPT')

Add calculation to project

.. code:: ipython2

    nwchem_example.add_calc(nwchem_octane)

Set the structure of the calculation to octane

.. code:: ipython2

    nwchem_octane.strucC = octane

Set the resource to be local

.. code:: ipython2

    nwchem_octane.set_resource(res_local)

Make the local directories

.. code:: ipython2

    nwchem_octane.make_dir()

Change to the ``scratch`` directory

.. code:: ipython2

    os.chdir(nwchem_octane.dir['scratch'])

Copy the template files to the scratch directory

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    nwchem_octane.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    nwchem_octane.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Read in the template files and add them to the ``str`` dictionary

.. code:: ipython2

    nwchem_octane.load_str('templates','nw')        
    nwchem_octane.load_str('templates','run')

Set the properties dictionary to desired calculation details

.. code:: ipython2

    nwchem_octane.properties['basis'] = '6-31g'
    nwchem_octane.properties['method'] = 'UHF'
    nwchem_octane.properties['charge'] = 0
    nwchem_octane.properties['spin_mult'] = 1
    nwchem_octane.properties['task'] = 'SCF optimize'
    nwchem_octane.properties['coord'] = nwchem_octane.strucC.write_coord()

.. code:: ipython2

    pprint(nwchem_octane.properties)

Replace the keys in the template strings and write the input files

.. code:: ipython2

    nwchem_octane.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_octane.tag))

Add the input file to the properties to be written into the run file

.. code:: ipython2

    nwchem_octane.properties['input_nw'] = nwchem_octane.files['input']['nw']
    nwchem_octane.replacewrite_prop('run','scripts','run','%s.sh'%(nwchem_octane.tag))

Add the log file to the files dictionary

.. code:: ipython2

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_octane.tag)
    nwchem_octane.add_file(file_type,file_key,file_name)

Change back to the root directory and write a json file

.. code:: ipython2

    os.chdir(nwchem_example.dir['home'])
    nwchem_example.export_json()

Change back to scratch

.. code:: ipython2

    print nwchem_example.dir['scratch']

.. code:: ipython2

    os.chdir(nwchem_example.dir['scratch'])

Run the bash script for the calculation or submit the job to the cluster

.. code:: ipython2

    print nwchem_octane.tag

.. code:: ipython2

    nwchem_octane.run()

Check the status of all the calculations in the project

.. code:: ipython2

    nwchem_example.check()

Run the analysis

.. code:: ipython2

    nwchem_octane.analysis()

Tar and zip the results and copy them to a storage location

.. code:: ipython2

    nwchem_example.store()

Save json in home directory

.. code:: ipython2

    os.chdir(nwchem_example.dir['home'])
    nwchem_example.export_json()
