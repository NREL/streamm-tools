.. _nwchem_example:
  
nwchem_example
========================
 

.. code:: python

    import os 
    from pprint import pprint

Check that output from other examples has been generated

.. code:: python

    from pathlib2 import Path

.. code:: python

    need_files = ['ethane_struc.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run buildingblocks_example.ipynb")
            os.system("jupyter nbconvert --to python  buildingblocks_example.ipynb")
            os.system("python buildingblocks_example.py")

.. code:: python

    need_files = ['local_res.json','remote_res.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run resource_example.ipynb")
            os.system("jupyter nbconvert --to python  resource_example.ipynb")
            os.system("python resource_example.py")

.. code:: python

    import streamm

Now let's create project and resource to keep track of our work

.. code:: python

    nwchem_example = streamm.Project('nwchem_example')
    res_local = streamm.Resource('local')

.. code:: python

    res_local.make_dir()

Update relative location of templates directory

.. code:: python

    res_local.dir['templates'] =  os.path.join(res_local.dir['home'],'..','templates','')

Make sure this is the location of the templates directory that comes
with the streamm git repository https://github.com/NREL/streamm-tools

.. code:: python

    print res_local.dir['templates']


.. parsed-literal::

    /Users/tkemper/Development/streamm-tools/examples/../templates/


Create the local directories that will store our files

.. code:: python

    nwchem_example.make_dir()

Tell the project about our directories

.. code:: python

    nwchem_example.set_resource(res_local)

Read in the ethane structure created in the
buildingblocks\_example.ipynb example

.. code:: python

    ethane = streamm.Buildingblock('ethane')

.. code:: python

    ethane.import_json()

.. code:: python

    print ethane.print_properties()


.. parsed-literal::

     n_particles:8 
     n_bonds:7
     n_angles:12
     n_dihedrals:9
     n_impropers:0


Set the paramkeys so we can identify force field parameters later on

.. code:: python

    for pkey,p in ethane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC'

.. code:: python

    for pk,p in ethane.particles.iteritems():
        p.residue = 1
        p.resname = 'ETH'

Set ``rsite``'s to hydrogens to be replaced during join

.. code:: python

    ethane.particles[1].rsite = 'RH'
    ethane.particles[5].rsite = 'RH'

Run ``find_rsites()`` to populate ``func`` list

.. code:: python

    ethane.find_rsites()

.. code:: python

    print ethane.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom H (H) index:5 n_bonds:1] 
    


.. code:: python

    import copy

Create octane from ethane

Copy ethane to a new Buildingblock octane

.. code:: python

    octane = copy.deepcopy(ethane)

.. code:: python

    from streamm.structures.buildingblock import attach

Then attach 3 more ethanes to make an octane

.. code:: python

    for i in range(3):
        octane = attach(octane,ethane,'RH',1,'RH',0)

.. code:: python

    octane.tag = 'octane'

.. code:: python

    print octane.n_particles


.. parsed-literal::

    26


.. code:: python

    octane.write_xyz()

Update the tag

.. code:: python

    octane.tag = 'octane'

Rename the residue and resname for octane

.. code:: python

    for pk,p in octane.particles.iteritems():
        p.residue = 2
        p.resname = "OCT"

.. code:: python

    octane.write_xyz()

Create NWChem Calculation object

.. code:: python

    nwchem_octane = streamm.NWChem('nw_octane_OPT')

Add calculation to project

.. code:: python

    nwchem_example.add_calc(nwchem_octane)

Set the structure of the calculation to octane

.. code:: python

    nwchem_octane.strucC = octane

Set the resource to be local

.. code:: python

    nwchem_octane.set_resource(res_local)

Make the local directories

.. code:: python

    nwchem_octane.make_dir()

Change to the ``scratch`` directory

.. code:: python

    os.chdir(nwchem_octane.dir['scratch'])

Copy the template files to the scratch directory

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    nwchem_octane.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    nwchem_octane.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Read in the template files and add them to the ``str`` dictionary

.. code:: python

    nwchem_octane.load_str('templates','nw')        
    nwchem_octane.load_str('templates','run')

Set the properties dictionary to desired calculation details

.. code:: python

    nwchem_octane.properties['basis'] = '6-31g'
    nwchem_octane.properties['method'] = 'UHF'
    nwchem_octane.properties['charge'] = 0
    nwchem_octane.properties['spin_mult'] = 1
    nwchem_octane.properties['task'] = 'SCF optimize'
    nwchem_octane.properties['coord'] = nwchem_octane.strucC.write_coord()

.. code:: python

    pprint(nwchem_octane.properties)


.. parsed-literal::

    {u'allocation': u'',
     u'basis': '6-31g',
     u'charge': 0,
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     'coord': u'     C       4.02000000       0.00000000       3.79009235 \n     H       5.21999999      -0.00000000       3.79009235 \n     H       3.62000000       0.97979589       4.35577777 \n     H       3.62000000      -0.97979589       4.35577777 \n     C       3.57333333       0.00000000       2.52672823 \n     H       3.97333333      -0.97979589       1.96104281 \n     H       3.97333333       0.97979589       1.96104281 \n     C       2.23333333       0.00000000       2.52672823 \n     H       1.83333334       0.97979589       3.09241365 \n     H       1.83333334      -0.97979589       3.09241365 \n     C       1.78666667       0.00000000       1.26336412 \n     H       2.18666666      -0.97979589       0.69767869 \n     H       2.18666666       0.97979589       0.69767869 \n     C       0.44666667       0.00000000       1.26336412 \n     H       0.04666667       0.97979589       1.82904954 \n     H       0.04666667      -0.97979589       1.82904954 \n     C       0.00000000       0.00000000       0.00000000 \n     H       0.40000000      -0.97979589      -0.56568542 \n     H       0.40000000       0.97979589      -0.56568542 \n     C      -1.34000000       0.00000000       0.00000000 \n     H      -1.74000000       0.97979589       0.56568542 \n     H      -1.74000000      -0.97979589       0.56568542 \n     C      -1.78666667       0.00000000      -1.26336412 \n     H      -2.98666666       0.00000000      -1.26336412 \n     H      -1.38666667      -0.97979589      -1.82904954 \n     H      -1.38666667       0.97979589      -1.82904954 \n',
     u'exe_command': u'./',
     u'feature': u'24core',
     u'finish_str': u'Total times  cpu:',
     u'maxiter': 100,
     u'method': 'UHF',
     u'nodes': 1,
     u'nproc': 1,
     u'pmem': 1500,
     u'ppn': 1,
     u'queue': u'batch',
     'scratch': u'/Users/tkemper/Development/streamm-tools/examples/scratch/nw_octane_OPT/',
     u'spin_mult': 1,
     u'task': 'SCF optimize',
     'uncompress': 'tar -xzf ',
     u'walltime': 24}


Replace the keys in the template strings and write the input files

.. code:: python

    nwchem_octane.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_octane.tag))

Add the input file to the properties to be written into the run file

.. code:: python

    nwchem_octane.properties['input_nw'] = nwchem_octane.files['input']['nw']
    nwchem_octane.replacewrite_prop('run','scripts','run','%s.sh'%(nwchem_octane.tag))

Add the log file to the files dictionary

.. code:: python

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_octane.tag)
    nwchem_octane.add_file(file_type,file_key,file_name)

Change back to the root directory and write a json file

.. code:: python

    os.chdir(nwchem_example.dir['home'])
    nwchem_example.export_json()




.. parsed-literal::

    {u'calculations': {'nw_octane_OPT': u'nwchem'},
     u'meta': {'date': '2017-11-15T16:56:54.736823',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['local']}



Change back to scratch

.. code:: python

    print nwchem_example.dir['scratch']


.. parsed-literal::

    /Users/tkemper/Development/streamm-tools/examples


.. code:: python

    os.chdir(nwchem_example.dir['scratch'])

Run the bash script for the calculation or submit the job to the cluster

.. code:: python

    print nwchem_octane.tag


.. parsed-literal::

    nw_octane_OPT


.. code:: python

    nwchem_octane.run()

Check the status of all the calculations in the project

.. code:: python

    nwchem_example.check()


.. parsed-literal::

    Calculation nw_octane_OPT has status written


Run the analysis

.. code:: python

    nwchem_octane.analysis()


.. parsed-literal::

    File nw_octane_OPT.log not found 


Tar and zip the results and copy them to a storage location

.. code:: python

    nwchem_example.store()

Save json in home directory

.. code:: python

    os.chdir(nwchem_example.dir['home'])
    nwchem_example.export_json()




.. parsed-literal::

    {u'calculations': {'nw_octane_OPT': u'nwchem'},
     u'meta': {'date': '2017-11-15T16:56:54.736823',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['local']}



