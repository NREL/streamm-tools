.. _Methane_Remote:
  
Methane_Remote
========================
 

In this getting started example we will calculate the electronic
properties of methane with NWChem and the structural properties with
LAMMPS

.. code:: python

    import os 
    from pprint import pprint

.. code:: python

    import numpy as np
    import decimal

.. code:: python

    from pathlib2 import Path
    import os

.. code:: python

    import streamm

Now let’s create project and resource to keep track of our work

.. code:: python

    Methane_example = streamm.Project('Methane_example')

And a resource object to keep track of where our files are

.. code:: python

    need_files = ['local_res.json','remote_res.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run resource_example.ipynb")
            os.system("jupyter nbconvert --to python  resource_example.ipynb")
            os.system("python resource_example.py")

Load resources from resources example

.. code:: python

    res_local = streamm.Resource('local')

The calc resource can be changed to local or remote host resource

.. code:: python

    res_calc = streamm.Resource('remote')

.. code:: python

    res_local.import_json()
    res_calc.import_json()

.. code:: python

    pprint(res_calc.properties)


.. parsed-literal::

    {u'allocation': u'orgopv',
     u'e-mail': u'tkemper@nrel.gov',
     u'exe_command': u'qsub ',
     u'feature': u'24core',
     u'nodes': 1,
     u'nproc': 24,
     u'pmem': 1500,
     u'ppn': 24,
     u'queue': u'short',
     u'scratch': u'/scratch/tkemper',
     u'walltime': 4}


.. code:: python

    Methane_example.set_resource(res_calc)

.. code:: python

    res_calc.make_dir()

Create .xyz file using a molecular viewer, such as Avogadro
(https://avogadro.cc/) or explicitly as in the structure.ipynb example.

.. code:: python

    ME = streamm.Buildingblock('methane')

.. code:: python

    ME.read_xyz()

.. code:: python

    print(ME.write_xyz_str())


.. parsed-literal::

     5 
     methane 
         C       0.00000000       0.00000000       0.00000000 
         H       0.69282032       0.69282032       0.69282032 
         H      -0.69282032      -0.69282032       0.69282032 
         H      -0.69282032       0.69282032      -0.69282032 
         H       0.69282032      -0.69282032      -0.69282032 
    


Looks good let’s set up some calculations

.. code:: python

    calc_i = streamm.Gaussian('g_methane_HF')

.. code:: python

    Methane_example.add_calc(calc_i) # Add it to the project 

.. code:: python

    calc_i.strucC = ME               # set the strucC to the structure container 

.. code:: python

    print calc_i.tag


.. parsed-literal::

    g_methane_HF


Let’s use the remote resource to run this calculation

.. code:: python

    calc_i.set_resource(res_calc)

-  home : directory is the root directory for the calculation/project
-  template : directory of template input and run files to be modified
   to run the calculation
-  materials : directory to store structure files (.xyz)
-  launch : directory to temporarily store files before they are copied
   to the remote resource
-  scratch : directory to run the calculation
-  storage : directory to store completed calculation data

.. code:: python

    pprint(calc_i.dir)


.. parsed-literal::

    {u'home': u'/Users/tkemper/Development/streamm-tools/examples',
     u'launch': u'/Users/tkemper/Development/streamm-tools/examples/scratch/g_methane_HF/',
     u'materials': u'/Users/tkemper/Development/streamm-tools/examples/materials',
     u'scratch': u'/scratch/tkemper/g_methane_HF/',
     u'scripts': u'/Users/tkemper/Development/streamm-tools/examples/scripts',
     u'storage': u'/mss/users/tkemper/g_methane_HF/',
     u'templates': u'/Users/tkemper/Development/streamm-tools/examples/../templates/'}


.. code:: python

    calc_i.make_dir()

.. code:: python

    print calc_i.dir['launch']


.. parsed-literal::

    /Users/tkemper/Development/streamm-tools/examples/scratch/g_methane_HF/


.. code:: python

    os.chdir(calc_i.dir['launch'])

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "gaussian_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'com'
    file_name = "gaussian.com"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    calc_i.load_str('templates','com')        
    calc_i.load_str('templates','run')

.. code:: python

    calc_i.properties['commands'] = 'HF/3-21G SP'
    calc_i.properties['charge'] = 0
    calc_i.properties['spin_mult'] = 1
    calc_i.properties['coord'] = calc_i.strucC.write_coord()

.. code:: python

    calc_i.replacewrite_prop('com','input','com','%s.com'%(calc_i.tag))

.. code:: python

    calc_i.properties['input_com'] = calc_i.files['input']['com']
    calc_i.replacewrite_prop('run','scripts','run','%s.pbs'%(calc_i.tag))

.. code:: python

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(calc_i.tag)
    calc_i.add_file(file_type,file_key,file_name)

.. code:: python

    file_type = 'output'
    file_key = 'fchk'
    file_name = "%s.fchk"%(calc_i.tag)
    calc_i.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: python

    os.chdir(calc_i.dir['home'])
    Methane_example.export_json()




.. parsed-literal::

    {u'calculations': {'g_methane_HF': u'gaussian'},
     u'meta': {'date': '2017-11-15T17:01:08.993609',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['peregrine']}



.. code:: python

    os.chdir(calc_i.dir['launch'])

.. code:: python

    calc_i.push()

.. code:: python

    calc_i.run()

Cool. While that is in the queue let’s setup some more jobs

Let’s also run a NWChem calculation

.. code:: python

    nwchem_i = streamm.NWChem('nw_methane_HF')

.. code:: python

    Methane_example.add_calc(nwchem_i)

.. code:: python

    nwchem_i.strucC = ME

.. code:: python

    print nwchem_i.tag


.. parsed-literal::

    nw_methane_HF


.. code:: python

    nwchem_i.set_resource(res_calc)

.. code:: python

    pprint(nwchem_i.properties['scratch'])


.. parsed-literal::

    u'/scratch/tkemper/nw_methane_HF/'


.. code:: python

    nwchem_i.make_dir()

.. code:: python

    os.chdir(nwchem_i.dir['launch'])

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    nwchem_i.load_str('templates','nw')        
    nwchem_i.load_str('templates','run')

.. code:: python

    nwchem_i.properties['basis'] = '6-31g'
    nwchem_i.properties['method'] = 'UHF'
    nwchem_i.properties['charge'] = 0
    nwchem_i.properties['spin_mult'] = 1
    nwchem_i.properties['task'] = 'SCF '
    nwchem_i.properties['coord'] = nwchem_i.strucC.write_coord()

.. code:: python

    pprint(nwchem_i.properties)


.. parsed-literal::

    {u'allocation': u'orgopv',
     u'basis': '6-31g',
     u'charge': 0,
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     'coord': u'     C       0.00000000       0.00000000       0.00000000 \n     H       0.69282032       0.69282032       0.69282032 \n     H      -0.69282032      -0.69282032       0.69282032 \n     H      -0.69282032       0.69282032      -0.69282032 \n     H       0.69282032      -0.69282032      -0.69282032 \n',
     u'e-mail': u'tkemper@nrel.gov',
     u'exe_command': u'qsub ',
     u'feature': u'24core',
     u'finish_str': u'Total times  cpu:',
     u'maxiter': 100,
     u'method': 'UHF',
     u'nodes': 1,
     u'nproc': 24,
     u'pmem': 1500,
     u'ppn': 24,
     u'queue': u'short',
     u'scratch': u'/scratch/tkemper/nw_methane_HF/',
     u'spin_mult': 1,
     u'task': 'SCF ',
     'uncompress': 'tar -xzf ',
     u'walltime': 4}


.. code:: python

    nwchem_i.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_i.tag))

.. code:: python

    nwchem_i.properties['input_nw'] = nwchem_i.files['input']['nw']
    nwchem_i.replacewrite_prop('run','scripts','run','%s.pbs'%(nwchem_i.tag))

.. code:: python

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_i.tag)
    nwchem_i.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: python

    os.chdir(nwchem_i.dir['home'])
    Methane_example.export_json()




.. parsed-literal::

    {u'calculations': {'g_methane_HF': u'gaussian', 'nw_methane_HF': u'nwchem'},
     u'meta': {'date': '2017-11-15T17:01:08.993609',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['peregrine']}



.. code:: python

    os.chdir(nwchem_i.dir['launch'])

.. code:: python

    nwchem_i.push()

.. code:: python

    nwchem_i.run()

Okay we have a couple calculations now, so let’s check their status

.. code:: python

    Methane_example.check()


.. parsed-literal::

    Calculation g_methane_HF has status finished
    Calculation nw_methane_HF has status finished


Run the check() function until they show as finished

.. code:: python

    os.chdir(nwchem_i.dir['launch'])

Store the calculation in compressed files

.. code:: python

    nwchem_i.store()

Download the compressed output files

.. code:: python

    nwchem_i.pull()

.. code:: python

    nwchem_i.analysis()


.. parsed-literal::

    Running analysis on  nw_methane_HF.log


.. code:: python

    print nwchem_i.properties['alpha_energies']


.. parsed-literal::

    [-0.9041047, -0.5161086, -0.5161086, -0.5161086, 0.2264494, 0.2911283, 0.2911283, 0.2911283, 0.7887659, 0.7887659, 0.7887659, 1.085928, 1.135132, 1.135132]


.. code:: python

    print nwchem_i.properties['N_alpha_occ']


.. parsed-literal::

    4


.. code:: python

    os.chdir(nwchem_i.dir['home'])
    Methane_example.export_json()




.. parsed-literal::

    {u'calculations': {'g_methane_HF': u'gaussian', 'nw_methane_HF': u'nwchem'},
     u'meta': {'date': '2017-11-15T17:01:08.993609',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['peregrine']}



Neat!

Now let’s optimize the structure and calculate the ESP charges

.. code:: python

    nwchem_opt = streamm.NWChem('nw_methane_OPT')

.. code:: python

    Methane_example.add_calc(nwchem_opt)

.. code:: python

    nwchem_opt.strucC = ME

.. code:: python

    print nwchem_opt.tag


.. parsed-literal::

    nw_methane_OPT


.. code:: python

    nwchem_opt.set_resource(res_calc)

.. code:: python

    nwchem_opt.make_dir()

.. code:: python

    print nwchem_opt.dir['launch']


.. parsed-literal::

    /Users/tkemper/Development/streamm-tools/examples/scratch/nw_methane_OPT/


.. code:: python

    os.chdir(nwchem_opt.dir['launch'])

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_opt.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_opt.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    nwchem_opt.load_str('templates','nw')        
    nwchem_opt.load_str('templates','run')

.. code:: python

    nwchem_opt.properties['basis'] = '6-31g'
    nwchem_opt.properties['method'] = 'UHF'
    nwchem_opt.properties['charge'] = 0
    nwchem_opt.properties['spin_mult'] = 1
    nwchem_opt.properties['task'] = 'SCF optimize'
    nwchem_opt.properties['coord'] = nwchem_opt.strucC.write_coord()

.. code:: python

    pprint(nwchem_opt.properties)


.. parsed-literal::

    {u'allocation': u'orgopv',
     u'basis': '6-31g',
     u'charge': 0,
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     'coord': u'     C       0.00000000       0.00000000       0.00000000 \n     H       0.69282032       0.69282032       0.69282032 \n     H      -0.69282032      -0.69282032       0.69282032 \n     H      -0.69282032       0.69282032      -0.69282032 \n     H       0.69282032      -0.69282032      -0.69282032 \n',
     u'e-mail': u'tkemper@nrel.gov',
     u'exe_command': u'qsub ',
     u'feature': u'24core',
     u'finish_str': u'Total times  cpu:',
     u'maxiter': 100,
     u'method': 'UHF',
     u'nodes': 1,
     u'nproc': 24,
     u'pmem': 1500,
     u'ppn': 24,
     u'queue': u'short',
     u'scratch': u'/scratch/tkemper/nw_methane_OPT/',
     u'spin_mult': 1,
     u'task': 'SCF optimize',
     'uncompress': 'tar -xzf ',
     u'walltime': 4}


.. code:: python

    nwchem_opt.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_opt.tag))

.. code:: python

    nwchem_opt.properties['input_nw'] = nwchem_opt.files['input']['nw']
    nwchem_opt.replacewrite_prop('run','scripts','run','%s.pbs'%(nwchem_opt.tag))

.. code:: python

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_opt.tag)
    nwchem_opt.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: python

    os.chdir(nwchem_opt.dir['home'])
    Methane_example.export_json()




.. parsed-literal::

    {u'calculations': {'g_methane_HF': u'gaussian',
      'nw_methane_HF': u'nwchem',
      'nw_methane_OPT': u'nwchem'},
     u'meta': {'date': '2017-11-15T17:01:08.993609',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['peregrine']}



.. code:: python

    os.chdir(nwchem_opt.dir['launch'])

.. code:: python

    nwchem_opt.push()

.. code:: python

    nwchem_opt.run()

.. code:: python

    nwchem_opt.check()

.. code:: python

    print nwchem_opt.meta


.. parsed-literal::

    {'date': '2017-11-15T17:20:45.530701', 'status': 'written', 'resource': 'peregrine', 'software': u'nwchem'}


.. code:: python

    Methane_example.check()


.. parsed-literal::

    Calculation nw_methane_OPT has status finished
    Calculation g_methane_HF has status finished
    Calculation nw_methane_HF has status stored


Again wait until the calculations show finished

.. code:: python

    nwchem_opt.store()

.. code:: python

    nwchem_opt.pull()

.. code:: python

    nwchem_opt.analysis()


.. parsed-literal::

    Running analysis on  nw_methane_OPT.log


.. code:: python

    print nwchem_opt.strucC.write_xyz_str()


.. parsed-literal::

     5 
     methane 
         C       0.00000000      -0.00000001      -0.00000001 
         H       0.62474531       0.62474532       0.62474532 
         H      -0.62474532      -0.62474531       0.62474533 
         H      -0.62474533       0.62474533      -0.62474531 
         H       0.62474533      -0.62474533      -0.62474532 
    


.. code:: python

    os.chdir(nwchem_opt.dir['materials'])

.. code:: python

    nwchem_opt.strucC.tag = '{}_{}'.format(nwchem_opt.strucC.tag,nwchem_opt.tag)

.. code:: python

    nwchem_opt.strucC.write_xyz()

.. code:: python

    nwchem_esp = streamm.NWChem('nw_methane_ESP')

.. code:: python

    ME_OPT = streamm.Buildingblock('methane_nw_methane_OPT')

.. code:: python

    ME_OPT.read_xyz()

.. code:: python

    print(ME.write_xyz_str())


.. parsed-literal::

     5 
     methane_nw_methane_OPT 
         C       0.00000000      -0.00000001      -0.00000001 
         H       0.62474531       0.62474532       0.62474532 
         H      -0.62474532      -0.62474531       0.62474533 
         H      -0.62474533       0.62474533      -0.62474531 
         H       0.62474533      -0.62474533      -0.62474532 
    


.. code:: python

    Methane_example.add_calc(nwchem_esp)

.. code:: python

    nwchem_esp.strucC = ME_OPT

.. code:: python

    print nwchem_esp.tag


.. parsed-literal::

    nw_methane_ESP


.. code:: python

    nwchem_esp.set_resource(res_calc)

.. code:: python

    pprint(nwchem_esp.properties['scratch'])


.. parsed-literal::

    u'/scratch/tkemper/nw_methane_ESP/'


.. code:: python

    nwchem_esp.make_dir()

.. code:: python

    print nwchem_esp.dir['launch']


.. parsed-literal::

    /Users/tkemper/Development/streamm-tools/examples/scratch/nw_methane_ESP/


.. code:: python

    os.chdir(nwchem_esp.dir['launch'])

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_esp.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem_esp.nw"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_esp.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    nwchem_esp.load_str('templates','nw')        
    nwchem_esp.load_str('templates','run')

.. code:: python

    nwchem_esp.properties['basis'] = '6-31g'
    nwchem_esp.properties['method'] = 'UHF'
    nwchem_esp.properties['charge'] = 0
    nwchem_esp.properties['spin_mult'] = 1
    nwchem_esp.properties['task'] = 'SCF '
    nwchem_esp.properties['coord'] = nwchem_esp.strucC.write_coord()

.. code:: python

    pprint(nwchem_esp.properties)


.. parsed-literal::

    {u'allocation': u'orgopv',
     u'basis': '6-31g',
     u'charge': 0,
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     'coord': u'     C       0.00000000      -0.00000001      -0.00000001 \n     H       0.62474531       0.62474532       0.62474532 \n     H      -0.62474532      -0.62474531       0.62474533 \n     H      -0.62474533       0.62474533      -0.62474531 \n     H       0.62474533      -0.62474533      -0.62474532 \n',
     u'e-mail': u'tkemper@nrel.gov',
     u'exe_command': u'qsub ',
     u'feature': u'24core',
     u'finish_str': u'Total times  cpu:',
     'input_nw': 'nw_methane_ESP.nw',
     u'maxiter': 100,
     u'method': 'UHF',
     u'nodes': 1,
     u'nproc': 24,
     u'pmem': 1500,
     u'ppn': 24,
     u'queue': u'short',
     u'scratch': u'/scratch/tkemper/nw_methane_ESP/',
     u'spin_mult': 1,
     u'task': 'SCF ',
     'uncompress': 'tar -xzf ',
     u'walltime': 4}


.. code:: python

    nwchem_esp.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_esp.tag))

.. code:: python

    nwchem_esp.properties['input_nw'] = nwchem_esp.files['input']['nw']
    nwchem_esp.replacewrite_prop('run','scripts','run','%s.pbs'%(nwchem_esp.tag))

.. code:: python

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_esp.tag)
    nwchem_esp.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: python

    os.chdir(nwchem_esp.dir['home'])
    Methane_example.export_json()




.. parsed-literal::

    {u'calculations': {'g_methane_HF': u'gaussian',
      'nw_methane_ESP': u'nwchem',
      'nw_methane_HF': u'nwchem',
      'nw_methane_OPT': u'nwchem'},
     u'meta': {'date': '2017-11-15T17:01:08.993609',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['peregrine']}



.. code:: python

    os.chdir(nwchem_esp.dir['launch'])

.. code:: python

    nwchem_esp.push()

.. code:: python

    nwchem_esp.run()

.. code:: python

    print nwchem_esp.tag,nwchem_esp.files['output']


.. parsed-literal::

    nw_methane_ESP {'log': 'nw_methane_ESP.log'}


.. code:: python

    nwchem_esp.check()

.. code:: python

    print nwchem_esp.meta


.. parsed-literal::

    {'date': '2017-11-15T17:24:28.884272', 'status': 'written', 'resource': 'peregrine', 'software': u'nwchem'}


.. code:: python

    Methane_example.check()


.. parsed-literal::

    Calculation nw_methane_OPT has status stored
    Calculation g_methane_HF has status stored
    Calculation nw_methane_HF has status stored
    Calculation nw_methane_ESP has status finished


.. code:: python

    nwchem_esp.store()

.. code:: python

    nwchem_esp.pull()

.. code:: python

    nwchem_esp.analysis()


.. parsed-literal::

    Running analysis on  nw_methane_ESP.log


.. code:: python

    for pk,p in nwchem_esp.strucC.particles.iteritems():
        print pk,p.charge


.. parsed-literal::

    0 -0.315785
    1 0.074534
    2 0.080417
    3 0.080417
    4 0.080417


Now we have an optimized molecular geometry and ESP charges
