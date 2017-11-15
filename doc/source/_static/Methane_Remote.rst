.. _Methane_Remote:
  
Methane_Remote
===============
 

.. code:: ipython2

    import os 
    from pprint import pprint

.. code:: ipython2

    import numpy as np
    import decimal

.. code:: ipython2

    from pathlib2 import Path
    import os

.. code:: ipython2

    import streamm

In this getting started example we will calculate the electronic
properties of methane with NWChem and the structural properties with
LAMMPS

Now let’s create project and resource to keep track of our work

.. code:: ipython2

    Methane_example = streamm.Project('Methane_example')

And a resource object to keep track of where our files are

.. code:: ipython2

    need_files = ['local_res.json','remote_res.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run resource_example.ipynb")
            os.system("jupyter nbconvert --to python  resource_example.ipynb")
            os.system("python resource_example.py")

Load resources from resources example

.. code:: ipython2

    res_local = streamm.Resource('local')

The calc resource can be changed to local or remote host resource

.. code:: ipython2

    res_calc = streamm.Resource('res_calc')

.. code:: ipython2

    res_local.import_json()
    res_calc.import_json()

.. code:: ipython2

    pprint(res_calc.properties)

.. code:: ipython2

    Methane_example.set_resource(res_calc)

.. code:: ipython2

    res_calc.make_dir()

Create .xyz file using a molecular viewer, such as Avogadro
(https://avogadro.cc/) or explicitly as in the structure.ipynb example.

.. code:: ipython2

    ME = streamm.Buildingblock('methane')

.. code:: ipython2

    ME.read_xyz()

.. code:: ipython2

    print(ME.write_xyz_str())

Looks good let’s set up some calculations

.. code:: ipython2

    calc_i = streamm.Gaussian('g_methane_HF')

.. code:: ipython2

    Methane_example.add_calc(calc_i) # Add it to the project 

.. code:: ipython2

    calc_i.strucC = ME               # set the strucC to the structure container 

.. code:: ipython2

    print calc_i.tag

Let’s use the remote resource to run this calculation

.. code:: ipython2

    calc_i.set_resource(res_calc)

-  home : directory is the root directory for the calculation/project
-  template : directory of template input and run files to be modified
   to run the calculation
-  materials : directory to store structure files (.xyz)
-  launch : directory to temporarily store files before they are copied
   to the remote resource
-  scratch : directory to run the calculation
-  storage : directory to store completed calculation data

.. code:: ipython2

    pprint(calc_i.dir)

.. code:: ipython2

    calc_i.make_dir()

.. code:: ipython2

    print calc_i.dir['launch']

.. code:: ipython2

    os.chdir(calc_i.dir['launch'])

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "gaussian_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    file_type = 'templates'
    file_key = 'com'
    file_name = "gaussian.com"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    calc_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    calc_i.load_str('templates','com')        
    calc_i.load_str('templates','run')

.. code:: ipython2

    calc_i.properties['commands'] = 'HF/3-21G SP'
    calc_i.properties['charge'] = 0
    calc_i.properties['spin_mult'] = 1
    calc_i.properties['coord'] = calc_i.strucC.write_coord()

.. code:: ipython2

    calc_i.replacewrite_prop('com','input','com','%s.com'%(calc_i.tag))

.. code:: ipython2

    calc_i.properties['input_com'] = calc_i.files['input']['com']
    calc_i.replacewrite_prop('run','scripts','run','%s.pbs'%(calc_i.tag))

.. code:: ipython2

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(calc_i.tag)
    calc_i.add_file(file_type,file_key,file_name)

.. code:: ipython2

    file_type = 'output'
    file_key = 'fchk'
    file_name = "%s.fchk"%(calc_i.tag)
    calc_i.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: ipython2

    os.chdir(calc_i.dir['home'])
    Methane_example.export_json()

.. code:: ipython2

    os.chdir(calc_i.dir['launch'])

.. code:: ipython2

    calc_i.push()

.. code:: ipython2

    calc_i.run()

Cool. While that is in the queue let’s setup some more jobs

Let’s also run a NWChem calculation

.. code:: ipython2

    nwchem_i = streamm.NWChem('nw_methane_HF')

.. code:: ipython2

    Methane_example.add_calc(nwchem_i)

.. code:: ipython2

    nwchem_i.strucC = ME

.. code:: ipython2

    print nwchem_i.tag

.. code:: ipython2

    nwchem_i.set_resource(res_calc)

.. code:: ipython2

    pprint(nwchem_i.properties['scratch'])

.. code:: ipython2

    nwchem_i.make_dir()

.. code:: ipython2

    print nwchem_i.dir['launch']

.. code:: ipython2

    os.chdir(nwchem_i.dir['launch'])

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    nwchem_i.load_str('templates','nw')        
    nwchem_i.load_str('templates','run')

.. code:: ipython2

    nwchem_i.properties['basis'] = '6-31g'
    nwchem_i.properties['method'] = 'UHF'
    nwchem_i.properties['charge'] = 0
    nwchem_i.properties['spin_mult'] = 1
    nwchem_i.properties['task'] = 'SCF '
    nwchem_i.properties['coord'] = nwchem_i.strucC.write_coord()

.. code:: ipython2

    pprint(nwchem_i.properties)

.. code:: ipython2

    nwchem_i.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_i.tag))

.. code:: ipython2

    nwchem_i.properties['input_nw'] = nwchem_i.files['input']['nw']
    nwchem_i.replacewrite_prop('run','scripts','run','%s.pbs'%(nwchem_i.tag))

.. code:: ipython2

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_i.tag)
    nwchem_i.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: ipython2

    os.chdir(nwchem_i.dir['home'])
    Methane_example.export_json()

.. code:: ipython2

    os.chdir(nwchem_i.dir['launch'])

.. code:: ipython2

    nwchem_i.push()

.. code:: ipython2

    nwchem_i.run()

Okay we have a couple calculations now, so let’s check their status

.. code:: ipython2

    Methane_example.check()

.. code:: ipython2

    nwchem_i.analysis()

.. code:: ipython2

    print nwchem_i.properties['alpha_energies']

.. code:: ipython2

    print nwchem_i.properties['N_alpha_occ']

.. code:: ipython2

    Methane_example.store()

.. code:: ipython2

    Methane_example.pull()

.. code:: ipython2

    os.chdir(nwchem_i.dir['home'])
    Methane_example.export_json()

Neat!

Now let’s optimize the structure and calculate the ESP charges

.. code:: ipython2

    nwchem_opt = streamm.NWChem('nw_methane_OPT')

.. code:: ipython2

    Methane_example.add_calc(nwchem_opt)

.. code:: ipython2

    nwchem_opt.strucC = ME

.. code:: ipython2

    print nwchem_opt.tag

.. code:: ipython2

    nwchem_opt.set_resource(res_calc)

.. code:: ipython2

    pprint(nwchem_opt.properties['scratch'])

.. code:: ipython2

    nwchem_opt.make_dir()

.. code:: ipython2

    print nwchem_opt.dir['launch']

.. code:: ipython2

    os.chdir(nwchem_opt.dir['launch'])

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_opt.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_opt.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    nwchem_opt.load_str('templates','nw')        
    nwchem_opt.load_str('templates','run')

.. code:: ipython2

    nwchem_opt.properties['basis'] = '6-31g'
    nwchem_opt.properties['method'] = 'UHF'
    nwchem_opt.properties['charge'] = 0
    nwchem_opt.properties['spin_mult'] = 1
    nwchem_opt.properties['task'] = 'SCF optimize'
    nwchem_opt.properties['coord'] = nwchem_opt.strucC.write_coord()

.. code:: ipython2

    pprint(nwchem_opt.properties)

.. code:: ipython2

    nwchem_opt.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_opt.tag))

.. code:: ipython2

    nwchem_opt.properties['input_nw'] = nwchem_opt.files['input']['nw']
    nwchem_opt.replacewrite_prop('run','scripts','run','%s.pbs'%(nwchem_opt.tag))

.. code:: ipython2

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_opt.tag)
    nwchem_opt.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: ipython2

    os.chdir(nwchem_opt.dir['home'])
    Methane_example.export_json()

.. code:: ipython2

    os.chdir(nwchem_opt.dir['launch'])

.. code:: ipython2

    nwchem_opt.push()

.. code:: ipython2

    nwchem_opt.run()

.. code:: ipython2

    nwchem_opt.check()

.. code:: ipython2

    print nwchem_opt.meta

.. code:: ipython2

    Methane_example.check()

.. code:: ipython2

    Methane_example.store()

.. code:: ipython2

    Methane_example.pull()

.. code:: ipython2

    nwchem_opt.analysis()

.. code:: ipython2

    print nwchem_opt.strucC.write_xyz_str()

.. code:: ipython2

    os.chdir(nwchem_opt.dir['materials'])

.. code:: ipython2

    nwchem_opt.strucC.tag = '{}_{}'.format(nwchem_opt.strucC.tag,nwchem_opt.tag)

.. code:: ipython2

    nwchem_opt.strucC.write_xyz()

.. code:: ipython2

    nwchem_esp = streamm.NWChem('nw_methane_ESP')

.. code:: ipython2

    ME_OPT = streamm.Buildingblock('methane_nw_methane_OPT')

.. code:: ipython2

    ME_OPT.read_xyz()

.. code:: ipython2

    print(ME.write_xyz_str())

.. code:: ipython2

    Methane_example.add_calc(nwchem_esp)

.. code:: ipython2

    nwchem_esp.strucC = ME_OPT

.. code:: ipython2

    print nwchem_esp.tag

.. code:: ipython2

    nwchem_esp.set_resource(res_calc)

.. code:: ipython2

    pprint(nwchem_esp.properties['scratch'])

.. code:: ipython2

    nwchem_esp.make_dir()

.. code:: ipython2

    print nwchem_esp.dir['launch']

.. code:: ipython2

    os.chdir(nwchem_esp.dir['launch'])

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem_remote.pbs"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_esp.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'launch'
    nwchem_esp.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    nwchem_esp.load_str('templates','nw')        
    nwchem_esp.load_str('templates','run')

.. code:: ipython2

    nwchem_esp.properties['basis'] = '6-31g'
    nwchem_esp.properties['method'] = 'UHF'
    nwchem_esp.properties['charge'] = 0
    nwchem_esp.properties['spin_mult'] = 1
    nwchem_esp.properties['task'] = 'esp'
    nwchem_esp.properties['coord'] = nwchem_esp.strucC.write_coord()

.. code:: ipython2

    pprint(nwchem_esp.properties)

.. code:: ipython2

    nwchem_esp.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_esp.tag))

.. code:: ipython2

    nwchem_esp.properties['input_nw'] = nwchem_esp.files['input']['nw']
    nwchem_esp.replacewrite_prop('run','scripts','run','%s.pbs'%(nwchem_esp.tag))

.. code:: ipython2

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_esp.tag)
    nwchem_esp.add_file(file_type,file_key,file_name)

Save details in .json files

.. code:: ipython2

    os.chdir(nwchem_esp.dir['home'])
    Methane_example.export_json()

.. code:: ipython2

    os.chdir(nwchem_esp.dir['launch'])

.. code:: ipython2

    nwchem_esp.push()

.. code:: ipython2

    nwchem_esp.run()

.. code:: ipython2

    print nwchem_esp.tag,nwchem_esp.files['output']

.. code:: ipython2

    nwchem_esp.check()

.. code:: ipython2

    print nwchem_esp.meta

.. code:: ipython2

    Methane_example.check()

.. code:: ipython2

    Methane_example.store()

.. code:: ipython2

    Methane_example.pull()

.. code:: ipython2

    nwchem_opt.analysis()

Now we have an optimized molecular geometry and ESP charges
