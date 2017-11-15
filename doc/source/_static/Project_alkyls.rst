.. _Project_alkyls:
  
Project_alkyls
===============
 

In this example we will create alkyl chains of various lengths, run
quantum chemical analysis on each and then replicate them into a
simulation cell for an MD simulation

.. code:: ipython2

    import os 
    from pprint import pprint

Check that output from other examples has been generated

.. code:: ipython2

    from pathlib2 import Path

.. code:: ipython2

    need_files = ['methane.xyz']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run structures_example.ipynb")
            os.system("jupyter nbconvert --to python  structures_example.ipynb")
            os.system("python structures_example.py")


.. parsed-literal::

    Need to run structures.ipynb


.. code:: ipython2

    need_files = ['ethane.xyz']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run buildingblocks_example.ipynb")
            os.system("jupyter nbconvert --to python  buildingblocks_example.ipynb")
            os.system("python buildingblocks_example.py")


.. parsed-literal::

    Need to run buildingblocks.ipynb


.. code:: ipython2

    need_files = ['oplsaa_param.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run forcefields_example.ipynb")
            os.system("jupyter nbconvert --to python  forcefields_example.ipynb")
            os.system("python forcefields_example.py")

.. code:: ipython2

    import streamm

Now let’s create project and resource to keep track of our work

.. code:: ipython2

    alkyl_example = streamm.Project('alkyl_example')
    res_local = streamm.Resource('local')

Update relative location of templates directory

.. code:: ipython2

    res_local.dir['templates'] =  os.path.join(res_local.dir['home'],'..','templates','')

Make sure this is the location of the templates directory that comes
with the streamm git repository https://github.com/NREL/streamm-tools

.. code:: ipython2

    print res_local.dir['templates']


.. parsed-literal::

    /Users/tkemper/Development/streamm-tools/examples/../templates/


Create the local directories that will store our files

.. code:: ipython2

    res_local.make_dir()

Tell the project about our directories

.. code:: ipython2

    alkyl_example.set_resource(res_local)

Read in the methane.xyz file created in the structure_example.ipynb
example

.. code:: ipython2

    methane = streamm.Buildingblock('methane')

.. code:: ipython2

    methane.read_xyz()

Create the neighbor list

.. code:: ipython2

    methane.bonded_nblist = methane.guess_nblist(0,radii_buffer=1.25)

and the bonded interactions

.. code:: ipython2

    methane.bonded_bonds()
    methane.bonded_angles()
    methane.bonded_dih()

.. code:: ipython2

    print methane.n_particles


.. parsed-literal::

    5


.. code:: ipython2

    print methane.print_properties()


.. parsed-literal::

     n_particles:5 
     n_bonds:4
     n_angles:6
     n_dihedrals:0
     n_impropers:0


Set the paramkeys so we can identify force field parameters later on

.. code:: ipython2

    for pkey,p in methane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC'

.. code:: ipython2

    for pk,p in methane.particles.iteritems():
        p.residue = 1
        p.resname = 'METH'

Set some rsites to be able to join molecules together

.. code:: ipython2

    methane.particles[1].rsite = 'RH'
    methane.particles[2].rsite = 'RH'

.. code:: ipython2

    methane.find_rsites()

.. code:: ipython2

    print methane.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom H (H) index:2 n_bonds:1] 
    


Read in ethane.xyz from the buildinblock_example.ipynb example

.. code:: ipython2

    ethane = streamm.Buildingblock('ethane')

.. code:: ipython2

    ethane.read_xyz()

Guess bonded neighbor list based on ``bonded_radii``

.. code:: ipython2

    ethane.bonded_nblist = ethane.guess_nblist(0,radii_buffer=1.25)

.. code:: ipython2

    ethane.bonded_bonds()
    ethane.bonded_angles()
    ethane.bonded_dih()

.. code:: ipython2

    print ethane.print_properties()


.. parsed-literal::

     n_particles:8 
     n_bonds:7
     n_angles:12
     n_dihedrals:9
     n_impropers:0


Set the ``paramkey``\ ’s as described in the force field example

.. code:: ipython2

    for pkey,p in ethane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC'

Set the ``resname`` of each particle to ``ETH``

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


.. parsed-literal::

    rsite:RH[ paticle:atom H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom H (H) index:5 n_bonds:1] 
    


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

Print new ``rsite``\ ’s

.. code:: ipython2

    print octane.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom H (H) index:23 n_bonds:1] 
    


Find the 4th carbon to attach an ethane

.. code:: ipython2

    print octane.particles[14].symbol


.. parsed-literal::

    H


.. code:: ipython2

    octane.particles[14].rsite = 'R2'

.. code:: ipython2

    octane.find_rsites()

Attach the ethane to the fourth carbon to make 4-ethyloctane

.. code:: ipython2

    ethyl_octane = attach(octane,ethane,'R2',0,'RH',0)

.. code:: ipython2

    ethyl_octane.tag = '4-ethyloctane'

.. code:: ipython2

    ethyl_octane.write_xyz()

Read in pickled oplsaa parameters from forcefield example

.. code:: ipython2

    oplsaa = streamm.forcefields.parameters.Parameters('oplsaa')

.. code:: ipython2

    oplsaa.import_json()

.. code:: ipython2

    print oplsaa


.. parsed-literal::

    
        Parameters 
          LJ parameters 2 
          Bond parameters 2 
          Angle parameters 2 
          Dihedral parameters 1 
          Improper Dihedral parameters 0 
    


.. code:: ipython2

    for pk,ptypes in oplsaa.particletypes.iteritems():
        print ptypes.fftype1


.. parsed-literal::

    CT
    HC


Create NWChem Calculation object

.. code:: ipython2

    nwchem_i = streamm.NWChem('nw_ethane_HF')

Add calculation to project

.. code:: ipython2

    alkyl_example.add_calc(nwchem_i)

Set the structure of the calculation to ethane

.. code:: ipython2

    nwchem_i.strucC = ethane

Set the resource to be local

.. code:: ipython2

    nwchem_i.set_resource(res_local)

Make the local directories

.. code:: ipython2

    nwchem_i.make_dir()

Change to the ``scratch`` directory

.. code:: ipython2

    os.chdir(nwchem_i.dir['scratch'])

Copy the template files to the scratch direcotry

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    nwchem_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    nwchem_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Read in the template files and add them to the ``str`` dictionary

.. code:: ipython2

    nwchem_i.load_str('templates','nw')        
    nwchem_i.load_str('templates','run')

Set the properties dictionary to desired calculation details

.. code:: ipython2

    nwchem_i.properties['basis'] = '6-31g'
    nwchem_i.properties['method'] = 'UHF'
    nwchem_i.properties['charge'] = 0
    nwchem_i.properties['spin_mult'] = 1
    nwchem_i.properties['task'] = 'SCF '
    nwchem_i.properties['coord'] = nwchem_i.strucC.write_coord()

.. code:: ipython2

    pprint(nwchem_i.properties)


.. parsed-literal::

    {u'allocation': u'',
     u'basis': '6-31g',
     u'charge': 0,
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     'coord': u'     C      64.69950113      18.57911440      15.46230149 \n     H      65.01275226      19.04910573      14.40353718 \n     H      64.61181440      17.38666991      15.36038405 \n     H      63.64100621      19.03038260      15.80280919 \n     C      65.62960687      18.88193760      16.37816251 \n     H      65.31635574      18.41194627      17.43692682 \n     H      65.71729360      20.07438209      16.48007995 \n     H      66.68810179      18.43066940      16.03765481 \n',
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
     'scratch': u'/Users/tkemper/Development/streamm-tools/examples/scratch/nw_ethane_HF/',
     u'spin_mult': 1,
     u'task': 'SCF ',
     'uncompress': 'tar -xzf ',
     u'walltime': 24}


Replace the keys in the template strings and write the input files

.. code:: ipython2

    nwchem_i.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_i.tag))

Add the input file to the properties to be written into the run file

.. code:: ipython2

    nwchem_i.properties['input_nw'] = nwchem_i.files['input']['nw']
    nwchem_i.replacewrite_prop('run','scripts','run','%s.sh'%(nwchem_i.tag))

Add the log file to the files dictionary

.. code:: ipython2

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_i.tag)
    nwchem_i.add_file(file_type,file_key,file_name)

Change back to the root directory and write a json file

.. code:: ipython2

    os.chdir(nwchem_i.dir['home'])
    alkyl_example.export_json()




.. parsed-literal::

    {u'calculations': {'gaus_ethane_HF': u'gaussian', 'nw_ethane_HF': u'nwchem'},
     u'meta': {'date': '2017-11-14T18:07:49.218138',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['local']}



Change back to scratch

.. code:: ipython2

    os.chdir(nwchem_i.dir['scratch'])

Run the bash script for the calculation or submit the job to the cluster

.. code:: ipython2

    nwchem_i.run()

Check the status of all the calculations in the project

.. code:: ipython2

    alkyl_example.check()


.. parsed-literal::

    Calculation nw_ethane_HF has status running
    Calculation gaus_ethane_HF has status running


Run the analysis

.. code:: ipython2

    nwchem_i.analysis()


.. parsed-literal::

    File nw_ethane_HF.log not found 


Tar and zip the results and copy them to a storage location

.. code:: ipython2

    nwchem_i.store()

Save json in home directory

.. code:: ipython2

    os.chdir(nwchem_i.dir['home'])
    alkyl_example.export_json()




.. parsed-literal::

    {u'calculations': {'gaus_ethane_HF': u'gaussian', 'nw_ethane_HF': u'nwchem'},
     u'meta': {'date': '2017-11-14T18:07:49.218138',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['local']}



Create a Gaussian Calculation object

.. code:: ipython2

    gaussian_i = streamm.Gaussian('gaus_ethane_HF')

Add the calculation to the project

.. code:: ipython2

    alkyl_example.add_calc(gaussian_i)

Set the structure of the calculation to ethane

.. code:: ipython2

    gaussian_i.strucC = ethane

Set the resource to be local

.. code:: ipython2

    gaussian_i.set_resource(res_local)

Make the local directories

.. code:: ipython2

    gaussian_i.make_dir()

Copy the template files to the scratch direcotry

.. code:: ipython2

    os.chdir(gaussian_i.dir['scratch'])

Copy the template files to the scratch direcotry

.. code:: ipython2

    file_type = 'templates'
    file_key = 'run'
    file_name = "gaussian.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    gaussian_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    file_type = 'templates'
    file_key = 'com'
    file_name = "gaussian.com"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    gaussian_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

Read in the template files and add them to the ``str`` dictionary

.. code:: ipython2

    gaussian_i.load_str('templates','com')        
    gaussian_i.load_str('templates','run')

Set the properties dictionary to desired calculation details

.. code:: ipython2

    gaussian_i.properties['commands'] = 'HF/3-21G SP'
    gaussian_i.properties['method'] = 'UHF'
    gaussian_i.properties['charge'] = 0
    gaussian_i.properties['spin_mult'] = 1
    gaussian_i.properties['coord'] = gaussian_i.strucC.write_coord()

.. code:: ipython2

    pprint(gaussian_i.properties)


.. parsed-literal::

    {u'allocation': u'',
     'charge': 0,
     'commands': 'HF/3-21G SP',
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     'coord': u'     C      64.69950113      18.57911440      15.46230149 \n     H      65.01275226      19.04910573      14.40353718 \n     H      64.61181440      17.38666991      15.36038405 \n     H      63.64100621      19.03038260      15.80280919 \n     C      65.62960687      18.88193760      16.37816251 \n     H      65.31635574      18.41194627      17.43692682 \n     H      65.71729360      20.07438209      16.48007995 \n     H      66.68810179      18.43066940      16.03765481 \n',
     u'exe_command': u'./',
     u'feature': u'24core',
     u'finish_str': u'Normal termination of Gaussian',
     'method': 'UHF',
     u'nodes': 1,
     u'nproc': 1,
     u'pmem': 1500,
     u'ppn': 1,
     u'queue': u'batch',
     'scratch': u'/Users/tkemper/Development/streamm-tools/examples/scratch/gaus_ethane_HF/',
     'spin_mult': 1,
     'uncompress': 'tar -xzf ',
     u'walltime': 24}


Replace the keys in the template strings and write the input files

.. code:: ipython2

    gaussian_i.replacewrite_prop('com','input','com','%s.com'%(gaussian_i.tag))

Add the input file to the properties to be written into the run file

.. code:: ipython2

    gaussian_i.properties['input_com'] = gaussian_i.files['input']['com']
    gaussian_i.replacewrite_prop('run','scripts','run','%s.sh'%(gaussian_i.tag))

Add the log file to the files dictionary

.. code:: ipython2

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(gaussian_i.tag)
    gaussian_i.add_file(file_type,file_key,file_name)

Change back to the root directory and write a json file

.. code:: ipython2

    os.chdir(gaussian_i.dir['home'])
    alkyl_example.export_json()




.. parsed-literal::

    {u'calculations': {'gaus_ethane_HF': u'gaussian', 'nw_ethane_HF': u'nwchem'},
     u'meta': {'date': '2017-11-14T18:07:49.218138',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['local']}



Change back to scratch

.. code:: ipython2

    os.chdir(gaussian_i.dir['scratch'])

Run the bash script for the calculation or submit the job to the cluster

.. code:: ipython2

    gaussian_i.run()

Check the status of all the calculations in the project

.. code:: ipython2

    alkyl_example.check()


.. parsed-literal::

    Calculation nw_ethane_HF has status running
    Calculation gaus_ethane_HF has status running


Run the analysis

.. code:: ipython2

    os.chdir(alkyl_example.dir['home'])
    alkyl_example.export_json()




.. parsed-literal::

    {u'calculations': {'gaus_ethane_HF': u'gaussian', 'nw_ethane_HF': u'nwchem'},
     u'meta': {'date': '2017-11-14T18:07:49.218138',
      'software': u'streamm_proj',
      'status': 'written'},
     u'resources': ['local']}



Create a LAMMPS Calculation object

.. code:: ipython2

    lmp_alkyl = streamm.LAMMPS('lmp_alkyl')

Turn periodic boundries on in all three directions

.. code:: ipython2

    lmp_alkyl.strucC.lat.pbcs = [True,True,True]

Run the ``add_struc()`` function to create 10 randomly placed
4-ethyloctane molecules

.. code:: ipython2

    seed = 92734
    lmp_alkyl.strucC = streamm.add_struc(lmp_alkyl.strucC,ethyl_octane,10,seed)


.. parsed-literal::

    No overlap found adding structure 0
    No overlap found adding structure 1
    No overlap found adding structure 2
    No overlap found adding structure 3
    No overlap found adding structure 4
    No overlap found adding structure 5
    No overlap found adding structure 6
    No overlap found adding structure 7
    No overlap found adding structure 8
    No overlap found adding structure 9


The ``add_struc()`` function randomly places each molecule in a space
defined by the lattice of the lmp_alkyl.strucC, then randomly rotates
it.

Then the function checks to make sure it does not overlap any other
particles that are already in the lmp_alkyl.strucC.

If an overlap is found a new position and rotation is chosen until the
max placements are exceeded, then the entire system is cleared, and the
placement starts again. If the maximum restarts are exceeded, then the
size of the lattice is increased, until all the molecules have been
added.

Check the lattice see if it expanded

.. code:: ipython2

    print lmp_alkyl.strucC.lat


.. parsed-literal::

    100.000000 0.000000 0.000000
    0.000000 100.000000 0.000000
    0.000000 0.000000 100.000000


Find the maximum molecule index

.. code:: ipython2

    print lmp_alkyl.strucC.n_molecules()


.. parsed-literal::

    9


.. code:: ipython2

    print ethyl_octane.tag


.. parsed-literal::

    4-ethyloctane


Update the structure tag

.. code:: ipython2

    lmp_alkyl.strucC.tag = ethyl_octane.tag + '_x10'

Write the structure to an xyz file

.. code:: ipython2

    lmp_alkyl.strucC.write_xyz()

Add 10 ethane to the structure container

.. code:: ipython2

    seed = 283674
    lmp_alkyl.strucC = streamm.add_struc(lmp_alkyl.strucC,ethane,10,seed)


.. parsed-literal::

    No overlap found adding structure 0
    No overlap found adding structure 1
    No overlap found adding structure 2
    No overlap found adding structure 3
    No overlap found adding structure 4
    No overlap found adding structure 5
    No overlap found adding structure 6
    No overlap found adding structure 7
    No overlap found adding structure 8
    No overlap found adding structure 9


.. code:: ipython2

    print lmp_alkyl.strucC.n_molecules()


.. parsed-literal::

    19


Update tag

.. code:: ipython2

    lmp_alkyl.strucC.tag += '_ethane_x10'

Add 50 methane to structure container using the ``add_struc_grid()``
which places solvent on grid

.. code:: ipython2

    lmp_alkyl.strucC = streamm.add_struc_grid(lmp_alkyl.strucC,methane,50)

Check to see if the lattice was expanded

.. code:: ipython2

    print lmp_alkyl.strucC.lat


.. parsed-literal::

    100.000000 0.000000 0.000000
    0.000000 100.000000 0.000000
    0.000000 0.000000 100.000000


Update tag

.. code:: ipython2

    lmp_alkyl.strucC.tag += '_methane_x50'

.. code:: ipython2

    lmp_alkyl.strucC.write_xyz()

Print all the particles in the structure container

.. code:: ipython2

    for pk,p in lmp_alkyl.strucC.particles.iteritems():
        print p,p.paramkey,p.mol,p.residue,p.resname


.. parsed-literal::

    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom H (H) HC 0 2 OCT
    atom C (C) CT 0 1 ETH
    atom H (H) HC 0 1 ETH
    atom H (H) HC 0 1 ETH
    atom C (C) CT 0 1 ETH
    atom H (H) HC 0 1 ETH
    atom H (H) HC 0 1 ETH
    atom H (H) HC 0 1 ETH
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom H (H) HC 1 2 OCT
    atom C (C) CT 1 1 ETH
    atom H (H) HC 1 1 ETH
    atom H (H) HC 1 1 ETH
    atom C (C) CT 1 1 ETH
    atom H (H) HC 1 1 ETH
    atom H (H) HC 1 1 ETH
    atom H (H) HC 1 1 ETH
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom H (H) HC 2 2 OCT
    atom C (C) CT 2 1 ETH
    atom H (H) HC 2 1 ETH
    atom H (H) HC 2 1 ETH
    atom C (C) CT 2 1 ETH
    atom H (H) HC 2 1 ETH
    atom H (H) HC 2 1 ETH
    atom H (H) HC 2 1 ETH
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom H (H) HC 3 2 OCT
    atom C (C) CT 3 1 ETH
    atom H (H) HC 3 1 ETH
    atom H (H) HC 3 1 ETH
    atom C (C) CT 3 1 ETH
    atom H (H) HC 3 1 ETH
    atom H (H) HC 3 1 ETH
    atom H (H) HC 3 1 ETH
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom H (H) HC 4 2 OCT
    atom C (C) CT 4 1 ETH
    atom H (H) HC 4 1 ETH
    atom H (H) HC 4 1 ETH
    atom C (C) CT 4 1 ETH
    atom H (H) HC 4 1 ETH
    atom H (H) HC 4 1 ETH
    atom H (H) HC 4 1 ETH
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom H (H) HC 5 2 OCT
    atom C (C) CT 5 1 ETH
    atom H (H) HC 5 1 ETH
    atom H (H) HC 5 1 ETH
    atom C (C) CT 5 1 ETH
    atom H (H) HC 5 1 ETH
    atom H (H) HC 5 1 ETH
    atom H (H) HC 5 1 ETH
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom H (H) HC 6 2 OCT
    atom C (C) CT 6 1 ETH
    atom H (H) HC 6 1 ETH
    atom H (H) HC 6 1 ETH
    atom C (C) CT 6 1 ETH
    atom H (H) HC 6 1 ETH
    atom H (H) HC 6 1 ETH
    atom H (H) HC 6 1 ETH
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom H (H) HC 7 2 OCT
    atom C (C) CT 7 1 ETH
    atom H (H) HC 7 1 ETH
    atom H (H) HC 7 1 ETH
    atom C (C) CT 7 1 ETH
    atom H (H) HC 7 1 ETH
    atom H (H) HC 7 1 ETH
    atom H (H) HC 7 1 ETH
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom H (H) HC 8 2 OCT
    atom C (C) CT 8 1 ETH
    atom H (H) HC 8 1 ETH
    atom H (H) HC 8 1 ETH
    atom C (C) CT 8 1 ETH
    atom H (H) HC 8 1 ETH
    atom H (H) HC 8 1 ETH
    atom H (H) HC 8 1 ETH
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom H (H) HC 9 2 OCT
    atom C (C) CT 9 1 ETH
    atom H (H) HC 9 1 ETH
    atom H (H) HC 9 1 ETH
    atom C (C) CT 9 1 ETH
    atom H (H) HC 9 1 ETH
    atom H (H) HC 9 1 ETH
    atom H (H) HC 9 1 ETH
    atom C (C) CT 10 1 ETH
    atom H (H) HC 10 1 ETH
    atom H (H) HC 10 1 ETH
    atom H (H) HC 10 1 ETH
    atom C (C) CT 10 1 ETH
    atom H (H) HC 10 1 ETH
    atom H (H) HC 10 1 ETH
    atom H (H) HC 10 1 ETH
    atom C (C) CT 11 1 ETH
    atom H (H) HC 11 1 ETH
    atom H (H) HC 11 1 ETH
    atom H (H) HC 11 1 ETH
    atom C (C) CT 11 1 ETH
    atom H (H) HC 11 1 ETH
    atom H (H) HC 11 1 ETH
    atom H (H) HC 11 1 ETH
    atom C (C) CT 12 1 ETH
    atom H (H) HC 12 1 ETH
    atom H (H) HC 12 1 ETH
    atom H (H) HC 12 1 ETH
    atom C (C) CT 12 1 ETH
    atom H (H) HC 12 1 ETH
    atom H (H) HC 12 1 ETH
    atom H (H) HC 12 1 ETH
    atom C (C) CT 13 1 ETH
    atom H (H) HC 13 1 ETH
    atom H (H) HC 13 1 ETH
    atom H (H) HC 13 1 ETH
    atom C (C) CT 13 1 ETH
    atom H (H) HC 13 1 ETH
    atom H (H) HC 13 1 ETH
    atom H (H) HC 13 1 ETH
    atom C (C) CT 14 1 ETH
    atom H (H) HC 14 1 ETH
    atom H (H) HC 14 1 ETH
    atom H (H) HC 14 1 ETH
    atom C (C) CT 14 1 ETH
    atom H (H) HC 14 1 ETH
    atom H (H) HC 14 1 ETH
    atom H (H) HC 14 1 ETH
    atom C (C) CT 15 1 ETH
    atom H (H) HC 15 1 ETH
    atom H (H) HC 15 1 ETH
    atom H (H) HC 15 1 ETH
    atom C (C) CT 15 1 ETH
    atom H (H) HC 15 1 ETH
    atom H (H) HC 15 1 ETH
    atom H (H) HC 15 1 ETH
    atom C (C) CT 16 1 ETH
    atom H (H) HC 16 1 ETH
    atom H (H) HC 16 1 ETH
    atom H (H) HC 16 1 ETH
    atom C (C) CT 16 1 ETH
    atom H (H) HC 16 1 ETH
    atom H (H) HC 16 1 ETH
    atom H (H) HC 16 1 ETH
    atom C (C) CT 17 1 ETH
    atom H (H) HC 17 1 ETH
    atom H (H) HC 17 1 ETH
    atom H (H) HC 17 1 ETH
    atom C (C) CT 17 1 ETH
    atom H (H) HC 17 1 ETH
    atom H (H) HC 17 1 ETH
    atom H (H) HC 17 1 ETH
    atom C (C) CT 18 1 ETH
    atom H (H) HC 18 1 ETH
    atom H (H) HC 18 1 ETH
    atom H (H) HC 18 1 ETH
    atom C (C) CT 18 1 ETH
    atom H (H) HC 18 1 ETH
    atom H (H) HC 18 1 ETH
    atom H (H) HC 18 1 ETH
    atom C (C) CT 19 1 ETH
    atom H (H) HC 19 1 ETH
    atom H (H) HC 19 1 ETH
    atom H (H) HC 19 1 ETH
    atom C (C) CT 19 1 ETH
    atom H (H) HC 19 1 ETH
    atom H (H) HC 19 1 ETH
    atom H (H) HC 19 1 ETH
    atom C (C) CT 20 1 METH
    atom H (H) HC 20 1 METH
    atom H (H) HC 20 1 METH
    atom H (H) HC 20 1 METH
    atom H (H) HC 20 1 METH
    atom C (C) CT 21 1 METH
    atom H (H) HC 21 1 METH
    atom H (H) HC 21 1 METH
    atom H (H) HC 21 1 METH
    atom H (H) HC 21 1 METH
    atom C (C) CT 22 1 METH
    atom H (H) HC 22 1 METH
    atom H (H) HC 22 1 METH
    atom H (H) HC 22 1 METH
    atom H (H) HC 22 1 METH
    atom C (C) CT 23 1 METH
    atom H (H) HC 23 1 METH
    atom H (H) HC 23 1 METH
    atom H (H) HC 23 1 METH
    atom H (H) HC 23 1 METH
    atom C (C) CT 24 1 METH
    atom H (H) HC 24 1 METH
    atom H (H) HC 24 1 METH
    atom H (H) HC 24 1 METH
    atom H (H) HC 24 1 METH
    atom C (C) CT 25 1 METH
    atom H (H) HC 25 1 METH
    atom H (H) HC 25 1 METH
    atom H (H) HC 25 1 METH
    atom H (H) HC 25 1 METH
    atom C (C) CT 26 1 METH
    atom H (H) HC 26 1 METH
    atom H (H) HC 26 1 METH
    atom H (H) HC 26 1 METH
    atom H (H) HC 26 1 METH
    atom C (C) CT 27 1 METH
    atom H (H) HC 27 1 METH
    atom H (H) HC 27 1 METH
    atom H (H) HC 27 1 METH
    atom H (H) HC 27 1 METH
    atom C (C) CT 28 1 METH
    atom H (H) HC 28 1 METH
    atom H (H) HC 28 1 METH
    atom H (H) HC 28 1 METH
    atom H (H) HC 28 1 METH
    atom C (C) CT 29 1 METH
    atom H (H) HC 29 1 METH
    atom H (H) HC 29 1 METH
    atom H (H) HC 29 1 METH
    atom H (H) HC 29 1 METH
    atom C (C) CT 30 1 METH
    atom H (H) HC 30 1 METH
    atom H (H) HC 30 1 METH
    atom H (H) HC 30 1 METH
    atom H (H) HC 30 1 METH
    atom C (C) CT 31 1 METH
    atom H (H) HC 31 1 METH
    atom H (H) HC 31 1 METH
    atom H (H) HC 31 1 METH
    atom H (H) HC 31 1 METH
    atom C (C) CT 32 1 METH
    atom H (H) HC 32 1 METH
    atom H (H) HC 32 1 METH
    atom H (H) HC 32 1 METH
    atom H (H) HC 32 1 METH
    atom C (C) CT 33 1 METH
    atom H (H) HC 33 1 METH
    atom H (H) HC 33 1 METH
    atom H (H) HC 33 1 METH
    atom H (H) HC 33 1 METH
    atom C (C) CT 34 1 METH
    atom H (H) HC 34 1 METH
    atom H (H) HC 34 1 METH
    atom H (H) HC 34 1 METH
    atom H (H) HC 34 1 METH
    atom C (C) CT 35 1 METH
    atom H (H) HC 35 1 METH
    atom H (H) HC 35 1 METH
    atom H (H) HC 35 1 METH
    atom H (H) HC 35 1 METH
    atom C (C) CT 36 1 METH
    atom H (H) HC 36 1 METH
    atom H (H) HC 36 1 METH
    atom H (H) HC 36 1 METH
    atom H (H) HC 36 1 METH
    atom C (C) CT 37 1 METH
    atom H (H) HC 37 1 METH
    atom H (H) HC 37 1 METH
    atom H (H) HC 37 1 METH
    atom H (H) HC 37 1 METH
    atom C (C) CT 38 1 METH
    atom H (H) HC 38 1 METH
    atom H (H) HC 38 1 METH
    atom H (H) HC 38 1 METH
    atom H (H) HC 38 1 METH
    atom C (C) CT 39 1 METH
    atom H (H) HC 39 1 METH
    atom H (H) HC 39 1 METH
    atom H (H) HC 39 1 METH
    atom H (H) HC 39 1 METH
    atom C (C) CT 40 1 METH
    atom H (H) HC 40 1 METH
    atom H (H) HC 40 1 METH
    atom H (H) HC 40 1 METH
    atom H (H) HC 40 1 METH
    atom C (C) CT 41 1 METH
    atom H (H) HC 41 1 METH
    atom H (H) HC 41 1 METH
    atom H (H) HC 41 1 METH
    atom H (H) HC 41 1 METH
    atom C (C) CT 42 1 METH
    atom H (H) HC 42 1 METH
    atom H (H) HC 42 1 METH
    atom H (H) HC 42 1 METH
    atom H (H) HC 42 1 METH
    atom C (C) CT 43 1 METH
    atom H (H) HC 43 1 METH
    atom H (H) HC 43 1 METH
    atom H (H) HC 43 1 METH
    atom H (H) HC 43 1 METH
    atom C (C) CT 44 1 METH
    atom H (H) HC 44 1 METH
    atom H (H) HC 44 1 METH
    atom H (H) HC 44 1 METH
    atom H (H) HC 44 1 METH
    atom C (C) CT 45 1 METH
    atom H (H) HC 45 1 METH
    atom H (H) HC 45 1 METH
    atom H (H) HC 45 1 METH
    atom H (H) HC 45 1 METH
    atom C (C) CT 46 1 METH
    atom H (H) HC 46 1 METH
    atom H (H) HC 46 1 METH
    atom H (H) HC 46 1 METH
    atom H (H) HC 46 1 METH
    atom C (C) CT 47 1 METH
    atom H (H) HC 47 1 METH
    atom H (H) HC 47 1 METH
    atom H (H) HC 47 1 METH
    atom H (H) HC 47 1 METH
    atom C (C) CT 48 1 METH
    atom H (H) HC 48 1 METH
    atom H (H) HC 48 1 METH
    atom H (H) HC 48 1 METH
    atom H (H) HC 48 1 METH
    atom C (C) CT 49 1 METH
    atom H (H) HC 49 1 METH
    atom H (H) HC 49 1 METH
    atom H (H) HC 49 1 METH
    atom H (H) HC 49 1 METH
    atom C (C) CT 50 1 METH
    atom H (H) HC 50 1 METH
    atom H (H) HC 50 1 METH
    atom H (H) HC 50 1 METH
    atom H (H) HC 50 1 METH
    atom C (C) CT 51 1 METH
    atom H (H) HC 51 1 METH
    atom H (H) HC 51 1 METH
    atom H (H) HC 51 1 METH
    atom H (H) HC 51 1 METH
    atom C (C) CT 52 1 METH
    atom H (H) HC 52 1 METH
    atom H (H) HC 52 1 METH
    atom H (H) HC 52 1 METH
    atom H (H) HC 52 1 METH
    atom C (C) CT 53 1 METH
    atom H (H) HC 53 1 METH
    atom H (H) HC 53 1 METH
    atom H (H) HC 53 1 METH
    atom H (H) HC 53 1 METH
    atom C (C) CT 54 1 METH
    atom H (H) HC 54 1 METH
    atom H (H) HC 54 1 METH
    atom H (H) HC 54 1 METH
    atom H (H) HC 54 1 METH
    atom C (C) CT 55 1 METH
    atom H (H) HC 55 1 METH
    atom H (H) HC 55 1 METH
    atom H (H) HC 55 1 METH
    atom H (H) HC 55 1 METH
    atom C (C) CT 56 1 METH
    atom H (H) HC 56 1 METH
    atom H (H) HC 56 1 METH
    atom H (H) HC 56 1 METH
    atom H (H) HC 56 1 METH
    atom C (C) CT 57 1 METH
    atom H (H) HC 57 1 METH
    atom H (H) HC 57 1 METH
    atom H (H) HC 57 1 METH
    atom H (H) HC 57 1 METH
    atom C (C) CT 58 1 METH
    atom H (H) HC 58 1 METH
    atom H (H) HC 58 1 METH
    atom H (H) HC 58 1 METH
    atom H (H) HC 58 1 METH
    atom C (C) CT 59 1 METH
    atom H (H) HC 59 1 METH
    atom H (H) HC 59 1 METH
    atom H (H) HC 59 1 METH
    atom H (H) HC 59 1 METH
    atom C (C) CT 60 1 METH
    atom H (H) HC 60 1 METH
    atom H (H) HC 60 1 METH
    atom H (H) HC 60 1 METH
    atom H (H) HC 60 1 METH
    atom C (C) CT 61 1 METH
    atom H (H) HC 61 1 METH
    atom H (H) HC 61 1 METH
    atom H (H) HC 61 1 METH
    atom H (H) HC 61 1 METH
    atom C (C) CT 62 1 METH
    atom H (H) HC 62 1 METH
    atom H (H) HC 62 1 METH
    atom H (H) HC 62 1 METH
    atom H (H) HC 62 1 METH
    atom C (C) CT 63 1 METH
    atom H (H) HC 63 1 METH
    atom H (H) HC 63 1 METH
    atom H (H) HC 63 1 METH
    atom H (H) HC 63 1 METH
    atom C (C) CT 64 1 METH
    atom H (H) HC 64 1 METH
    atom H (H) HC 64 1 METH
    atom H (H) HC 64 1 METH
    atom H (H) HC 64 1 METH
    atom C (C) CT 65 1 METH
    atom H (H) HC 65 1 METH
    atom H (H) HC 65 1 METH
    atom H (H) HC 65 1 METH
    atom H (H) HC 65 1 METH
    atom C (C) CT 66 1 METH
    atom H (H) HC 66 1 METH
    atom H (H) HC 66 1 METH
    atom H (H) HC 66 1 METH
    atom H (H) HC 66 1 METH
    atom C (C) CT 67 1 METH
    atom H (H) HC 67 1 METH
    atom H (H) HC 67 1 METH
    atom H (H) HC 67 1 METH
    atom H (H) HC 67 1 METH
    atom C (C) CT 68 1 METH
    atom H (H) HC 68 1 METH
    atom H (H) HC 68 1 METH
    atom H (H) HC 68 1 METH
    atom H (H) HC 68 1 METH
    atom C (C) CT 69 1 METH
    atom H (H) HC 69 1 METH
    atom H (H) HC 69 1 METH
    atom H (H) HC 69 1 METH
    atom H (H) HC 69 1 METH


Set ff parameters for all the bonds, bond angles and dihedrals in the
structure container

.. code:: ipython2

    lmp_alkyl.paramC = oplsaa

.. code:: ipython2

    lmp_alkyl.set_ffparam()

Add template files to calculations

.. code:: ipython2

    file_type = 'templates'
    file_key = 'in'
    file_name = "lammps_spneut.in"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    lmp_alkyl.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: ipython2

    pprint("Calculation:{} has status:{}".format(lmp_alkyl.tag,lmp_alkyl.meta['status']))


.. parsed-literal::

    'Calculation:lmp_alkyl has status:written'


Calculate the center mass of structure

.. code:: ipython2

    lmp_alkyl.strucC.calc_center_mass()

Create groups out of the molecules

.. code:: ipython2

    groupset_i = streamm.Groups('mol',lmp_alkyl.strucC)
    groupset_i.group_prop('mol','group_mol')

Caculate the ceneter of mass, radius and asphericity of each group

.. code:: ipython2

    groupset_i.calc_cent_mass()
    groupset_i.calc_radius_asphericity()
    groupset_i.calc_dl()

Write the center of mass of each group to an .xyz file for visulization

.. code:: ipython2

    groupset_i.write_cm_xyz()

.. code:: ipython2

    import numpy as np

.. code:: ipython2

    print np.mean(groupset_i.radius),groupset_i.strucC.unit_conf['length']


.. parsed-literal::

    1.79932546227 ang


.. code:: ipython2

    print groupset_i.strucC.lat.pbcs


.. parsed-literal::

    [True, True, True]


Create a neighbor list of groups

.. code:: ipython2

    groupset_i.group_nblist.radii_nblist(groupset_i.strucC.lat,groupset_i.cent_mass,groupset_i.radius,radii_buffer=5.25)

Apply periodic boundries to all the groups, so the molecules are not
split across pbc’s

.. code:: ipython2

    groupset_i.group_pbcs()

Loop over each group, shift the group to the center of the simulation
cell and write an .xyz file that includes the neighbors of the group.

.. code:: ipython2

    for gk_i,g_i in groupset_i.groups.iteritems():
        if( len(g_i.pkeys) == 32 ):
            print g_i.tag,groupset_i.group_nblist.calc_nnab(gk_i),g_i.mol 
            print g_i.cent_mass
            list_i = []
            for g_j in groupset_i.group_nblist.getnbs(gk_i):
                list_i += groupset_i.groups[g_j].pkeys
            groupset_i.strucC.shift_pos(-1.0*g_i.cent_mass)  # Place center of mass at origin
            groupset_i.strucC.write_xyz_list(list_i,xyz_file='{}_blob.xyz'.format(g_i.tag))
            groupset_i.strucC.shift_pos(g_i.cent_mass)  # Return center of mass 
            


.. parsed-literal::

    group_mol_0 16 0
    [ 76.653846  10.30594    0.623454]
    group_mol_1 13 1
    [ 73.198925  53.63841   74.637648]
    group_mol_2 18 2
    [  6.854253  86.197885  25.639487]
    group_mol_3 14 3
    [  4.771635  68.191421  89.863634]
    group_mol_4 16 4
    [ 25.9855    36.503382  59.348536]
    group_mol_5 14 5
    [ 98.939799  17.188399  88.616208]
    group_mol_6 12 6
    [ 34.513904   5.223367  49.861893]
    group_mol_7 9 7
    [ 66.215706  71.928804  48.970751]
    group_mol_8 15 8
    [ 98.116022  81.682274   1.600691]
    group_mol_9 13 9
    [ 36.312087  57.242611  93.212946]


Fancy aye!
