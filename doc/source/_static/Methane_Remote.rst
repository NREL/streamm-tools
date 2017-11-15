.. _Methane_Remote:
  
Methane_Remote
========================
 

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

In this getting started example we will calculate the electronic
properties of methane with NWChem and the structural properties with
LAMMPS

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

    res_calc = streamm.Resource('res_calc')

.. code:: python

    res_local.import_json()
    res_calc.import_json()

.. code:: python

    pprint(res_calc.properties)

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

Looks good let’s set up some calculations

.. code:: python

    calc_i = streamm.Gaussian('g_methane_HF')

.. code:: python

    Methane_example.add_calc(calc_i) # Add it to the project 

.. code:: python

    calc_i.strucC = ME               # set the strucC to the structure container 

.. code:: python

    print calc_i.tag

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

.. code:: python

    calc_i.make_dir()

.. code:: python

    print calc_i.dir['launch']

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

.. code:: python

    nwchem_i.set_resource(res_calc)

.. code:: python

    pprint(nwchem_i.properties['scratch'])

.. code:: python

    nwchem_i.make_dir()

.. code:: python

    print nwchem_i.dir['launch']

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

.. code:: python

    os.chdir(nwchem_i.dir['launch'])

.. code:: python

    nwchem_i.push()

.. code:: python

    nwchem_i.run()

Okay we have a couple calculations now, so let’s check their status

.. code:: python

    Methane_example.check()

.. code:: python

    nwchem_i.analysis()

.. code:: python

    print nwchem_i.properties['alpha_energies']

.. code:: python

    print nwchem_i.properties['N_alpha_occ']

.. code:: python

    Methane_example.store()

.. code:: python

    Methane_example.pull()

.. code:: python

    os.chdir(nwchem_i.dir['home'])
    Methane_example.export_json()

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

.. code:: python

    nwchem_opt.set_resource(res_calc)

.. code:: python

    pprint(nwchem_opt.properties['scratch'])

.. code:: python

    nwchem_opt.make_dir()

.. code:: python

    print nwchem_opt.dir['launch']

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

.. code:: python

    Methane_example.check()

.. code:: python

    Methane_example.store()

.. code:: python

    Methane_example.pull()

.. code:: python

    nwchem_opt.analysis()

.. code:: python

    print nwchem_opt.strucC.write_xyz_str()

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

.. code:: python

    Methane_example.add_calc(nwchem_esp)

.. code:: python

    nwchem_esp.strucC = ME_OPT

.. code:: python

    print nwchem_esp.tag

.. code:: python

    nwchem_esp.set_resource(res_calc)

.. code:: python

    pprint(nwchem_esp.properties['scratch'])

.. code:: python

    nwchem_esp.make_dir()

.. code:: python

    print nwchem_esp.dir['launch']

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
    file_name = "nwchem.nw"
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
    nwchem_esp.properties['task'] = 'esp'
    nwchem_esp.properties['coord'] = nwchem_esp.strucC.write_coord()

.. code:: python

    pprint(nwchem_esp.properties)

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

.. code:: python

    os.chdir(nwchem_esp.dir['launch'])

.. code:: python

    nwchem_esp.push()

.. code:: python

    nwchem_esp.run()

.. code:: python

    print nwchem_esp.tag,nwchem_esp.files['output']

.. code:: python

    nwchem_esp.check()

.. code:: python

    print nwchem_esp.meta

.. code:: python

    Methane_example.check()

.. code:: python

    Methane_example.store()

.. code:: python

    Methane_example.pull()

.. code:: python

    nwchem_opt.analysis()

Now we have an optimized molecular geometry and ESP charges
