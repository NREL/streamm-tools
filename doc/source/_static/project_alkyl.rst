.. _Project_alkyl:

Project_alkyl
=============

.. code:: python

    %load_ext autoreload
    %autoreload 2

In this example we will create alkyl chains of various lengths, run
quantum chemical analysis on each and then replicate them into a
simulation cell for an MD simulation

.. code:: python

    from __future__ import division, unicode_literals

.. code:: python

    import os 
    from pprint import pprint

.. code:: python

    import streamm

Now let's create project and resource to keep track of our work

.. code:: python

    alkyl_example = streamm.Project('alkyl_example')
    res_local = streamm.Resource('local')

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

    res_local.make_dir()

Tell the project about our directories

.. code:: python

    alkyl_example.set_resource(res_local)

Read in the methane.xyz file created in the structure.ipynb example

.. code:: python

    methane = streamm.Buildingblock('methane')

.. code:: python

    methane.read_xyz()

Create the neighbor list and the bonded interactions

.. code:: python

    methane.bonded_nblist = methane.guess_nblist(0,radii_buffer=1.25)

.. code:: python

    methane.bonded_bonds()
    methane.bonded_angles()
    methane.bonded_dih()

.. code:: python

    print methane.n_particles


.. parsed-literal::

    5


.. code:: python

    print methane.print_properties()


.. parsed-literal::

     n_particles:5 
     n_bonds:4
     n_angles:6
     n_dihedrals:0
     n_impropers:0


Set the paramkeys so we can identify force field paramters later on

.. code:: python

    for pkey,p in methane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC'

.. code:: python

    for pk,p in methane.particles.iteritems():
        p.residue = 1
        p.resname = 'METH'

Set some rsites to be able to join molecules together

.. code:: python

    methane.particles[1].rsite = 'RH'
    methane.particles[2].rsite = 'RH'

.. code:: python

    methane.find_rsites()

.. code:: python

    print methane.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom[1] H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom[2] H (H) index:2 n_bonds:1] 
    


Read in ethane.xyz from the buildinblock.ipynb example

.. code:: python

    ethane = streamm.Buildingblock('ethane')

.. code:: python

    ethane.read_xyz()

.. code:: python

    ethane.bonded_nblist = ethane.guess_nblist(0,radii_buffer=1.25)

.. code:: python

    ethane.bonded_bonds()
    ethane.bonded_angles()
    ethane.bonded_dih()

.. code:: python

    print ethane.n_particles


.. parsed-literal::

    8


.. code:: python

    print ethane.print_properties()


.. parsed-literal::

     n_particles:8 
     n_bonds:7
     n_angles:12
     n_dihedrals:9
     n_impropers:0


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

.. code:: python

    ethane.particles[1].rsite = 'RH'
    ethane.particles[5].rsite = 'RH'

.. code:: python

    ethane.find_rsites()

.. code:: python

    print ethane.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom[1] H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom[5] H (H) index:5 n_bonds:1] 
    


.. code:: python

    import copy

.. code:: python

    octane = copy.deepcopy(ethane)

.. code:: python

    from streamm.structures.buildingblock import attach

.. code:: python

    for i in range(3):
        octane = attach(octane,ethane,'RH',1,'RH',0)

.. code:: python

    octane.tag = 'octane'

.. code:: python

    for pk,p in octane.particles.iteritems():
        p.residue = 2
        p.resname = "OCT"
     

.. code:: python

    octane.write_xyz()

.. code:: python

    print octane.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom[1] H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom[23] H (H) index:23 n_bonds:1] 
    


.. code:: python

    print octane.particles[14].symbol


.. parsed-literal::

    H


.. code:: python

    octane.particles[14].rsite = 'R2'

.. code:: python

    octane.find_rsites()

.. code:: python

    ethyl_octane = attach(octane,ethane,'R2',0,'RH',0)

.. code:: python

    ethyl_octane.tag = 'ethyl_octane'

.. code:: python

    ethyl_octane.write_xyz()

.. code:: python

    for pk,p in ethyl_octane.particles.iteritems():
        print p.residue,p.resname
        #p.residue = 1
        #p.resname = 'ETHHEX'


.. parsed-literal::

    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    2 OCT
    1 ETH
    1 ETH
    1 ETH
    1 ETH
    1 ETH
    1 ETH
    1 ETH


.. code:: python

    oplsaa = streamm.forcefields.parameters.read_pickle('oplsaa')

.. code:: python

    print oplsaa


.. parsed-literal::

    
        Parameters 
          LJ parameters 2 
          Bond parameters 3 
          Angle parameters 2 
          Dihedral parameters 1 
          Imporper Dihedral parameters 0 
    


.. code:: python

    nwchem_i = streamm.NWChem('nw_ethane_HF')

.. code:: python

    alkyl_example.add_calc(nwchem_i)

.. code:: python

    nwchem_i.strucC = ethane

.. code:: python

    nwchem_i.set_resource(res_local)

.. code:: python

    nwchem_i.make_dir()

.. code:: python

    os.chdir(nwchem_i.dir['scratch'])

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "nwchem.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    nwchem_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'nw'
    file_name = "nwchem.nw"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
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

    {u'allocation': u'',
     u'basis': u'6-31g',
     u'charge': 0,
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     u'coord': u'     C       1.34000000      -0.00000000       0.00000000 \n     H       1.74000000      -0.00000000      -1.13137084 \n     H       1.74000000       0.97979589       0.56568542 \n     H       1.74000000      -0.97979589       0.56568542 \n     C       0.00000000       0.00000000       0.00000000 \n     H      -0.40000000       0.00000000       1.13137084 \n     H      -0.40000000      -0.97979589      -0.56568542 \n     H      -0.40000000       0.97979589      -0.56568542 \n',
     u'exe_command': u'./',
     u'feature': u'24core',
     u'finish_str': u'Total times  cpu:',
     u'method': u'UHF',
     u'nodes': 1,
     u'nproc': 1,
     u'pmem': 1500,
     u'ppn': 1,
     u'queue': u'batch',
     u'scratch': u'/Users/tkemper/Development/streamm-tools/examples/scratch/nw_ethane_HF/',
     u'spin_mult': 1,
     u'task': u'SCF ',
     'uncompress': 'tar -xzf ',
     u'walltime': 24}


.. code:: python

    nwchem_i.replacewrite_prop('nw','input','nw','%s.nw'%(nwchem_i.tag))

.. code:: python

    nwchem_i.properties['input_nw'] = nwchem_i.files['input']['nw']
    nwchem_i.replacewrite_prop('run','scripts','run','%s.sh'%(nwchem_i.tag))

.. code:: python

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(nwchem_i.tag)
    nwchem_i.add_file(file_type,file_key,file_name)

.. code:: python

    os.chdir(nwchem_i.dir['home'])
    alkyl_example.dump_json()

.. code:: python

    os.chdir(nwchem_i.dir['scratch'])

.. code:: python

    nwchem_i.run()

.. code:: python

    alkyl_example.check()


.. parsed-literal::

    Calculation nw_ethane_HF has status running


.. code:: python

    nwchem_i.analysis()

.. code:: python

    nwchem_i.store()

.. code:: python

    os.chdir(nwchem_i.dir['home'])
    alkyl_example.dump_json()

.. code:: python

    gaussian_i = streamm.Gaussian('gaus_ethane_HF')

.. code:: python

    alkyl_example.add_calc(gaussian_i)

.. code:: python

    gaussian_i.strucC = ethane

.. code:: python

    gaussian_i.set_resource(res_local)

.. code:: python

    gaussian_i.make_dir()

.. code:: python

    os.chdir(gaussian_i.dir['scratch'])

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "gaussian.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    gaussian_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'com'
    file_name = "gaussian.com"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    gaussian_i.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    gaussian_i.load_str('templates','com')        
    gaussian_i.load_str('templates','run')

.. code:: python

    gaussian_i.properties['commands'] = 'HF/3-21G SP'
    gaussian_i.properties['method'] = 'UHF'
    gaussian_i.properties['charge'] = 0
    gaussian_i.properties['spin_mult'] = 1
    gaussian_i.properties['coord'] = gaussian_i.strucC.write_coord()

.. code:: python

    pprint(gaussian_i.properties)


.. parsed-literal::

    {u'allocation': u'',
     u'charge': 0,
     u'commands': u'HF/3-21G SP',
     'comp_key': 'compressed',
     'compress': 'tar -czf ',
     'compress_sufix': 'tgz',
     u'coord': u'     C       1.34000000      -0.00000000       0.00000000 \n     H       1.74000000      -0.00000000      -1.13137084 \n     H       1.74000000       0.97979589       0.56568542 \n     H       1.74000000      -0.97979589       0.56568542 \n     C       0.00000000       0.00000000       0.00000000 \n     H      -0.40000000       0.00000000       1.13137084 \n     H      -0.40000000      -0.97979589      -0.56568542 \n     H      -0.40000000       0.97979589      -0.56568542 \n',
     u'exe_command': u'./',
     u'feature': u'24core',
     u'finish_str': u'Normal termination of Gaussian',
     u'method': u'UHF',
     u'nodes': 1,
     u'nproc': 1,
     u'pmem': 1500,
     u'ppn': 1,
     u'queue': u'batch',
     u'scratch': u'/Users/tkemper/Development/streamm-tools/examples/scratch/gaus_ethane_HF/',
     u'spin_mult': 1,
     'uncompress': 'tar -xzf ',
     u'walltime': 24}


.. code:: python

    gaussian_i.replacewrite_prop('com','input','com','%s.com'%(gaussian_i.tag))

.. code:: python

    gaussian_i.properties['input_com'] = gaussian_i.files['input']['com']
    gaussian_i.replacewrite_prop('run','scripts','run','%s.sh'%(gaussian_i.tag))

.. code:: python

    file_type = 'output'
    file_key = 'log'
    file_name = "%s.log"%(gaussian_i.tag)
    gaussian_i.add_file(file_type,file_key,file_name)

.. code:: python

    os.chdir(gaussian_i.dir['home'])
    alkyl_example.dump_json()

.. code:: python

    os.chdir(gaussian_i.dir['scratch'])

.. code:: python

    gaussian_i.run()

.. code:: python

    alkyl_example.check()


.. parsed-literal::

    Calculation nw_ethane_HF has status running
    Calculation gaus_ethane_HF has status running


.. code:: python

    os.chdir(alkyl_example.dir['home'])
    alkyl_example.dump_json()

.. code:: python

    lmp_alkyl = streamm.LAMMPS('lmp_alkyl')

.. code:: python

    alkyl_example.add_calc(lmp_alkyl)

.. code:: python

    lmp_alkyl.set_resource(res_local)

.. code:: python

    lmp_alkyl.make_dir()

.. code:: python

    print lmp_alkyl.dir['scratch']


.. parsed-literal::

    /Users/tkemper/Development/streamm-tools/examples/scratch/lmp_alkyl/


.. code:: python

    os.chdir(lmp_alkyl.dir['scratch'])

.. code:: python

    lmp_alkyl.paramC = oplsaa

.. code:: python

    import streamm

.. code:: python

    lmp_alkyl.strucC =  streamm.Buildingblock(matrix=[50.0,0.0,0.0,0.0,50.0,0.0,0.0,0.0,50.0])

.. code:: python

    lmp_alkyl.strucC.lat.pbcs = [True,True,True]

.. code:: python

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


.. code:: python

    print lmp_alkyl.strucC.lat


.. parsed-literal::

    50.000000 0.000000 0.000000
    0.000000 50.000000 0.000000
    0.000000 0.000000 50.000000


.. code:: python

    print lmp_alkyl.strucC.n_molecules()


.. parsed-literal::

    9


.. code:: python

    print ethyl_octane.tag


.. parsed-literal::

    ethyl_octane


.. code:: python

    lmp_alkyl.strucC.tag = ethyl_octane.tag + '_x10'

.. code:: python

    lmp_alkyl.strucC.write_xyz()

.. code:: python

    seed = 283674
    lmp_alkyl.strucC = streamm.add_struc(lmp_alkyl.strucC,ethane,10,seed)


.. parsed-literal::

    No overlap found adding structure 0
    No overlap found adding structure 1
    Max placments 10 exceeded resetting to original system 
    No overlap found adding structure 0
    No overlap found adding structure 1
    No overlap found adding structure 2
    No overlap found adding structure 3
    No overlap found adding structure 4
    Max placments 10 exceeded resetting to original system 
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


.. code:: python

    print lmp_alkyl.strucC.n_molecules()


.. parsed-literal::

    19


.. code:: python

    lmp_alkyl.strucC.tag += '_ethane_x10'

.. code:: python

    lmp_alkyl.strucC = streamm.add_struc_grid(lmp_alkyl.strucC,methane,50)

.. code:: python

    print lmp_alkyl.strucC.lat


.. parsed-literal::

    66.550000 0.000000 0.000000
    0.000000 66.550000 0.000000
    0.000000 0.000000 66.550000


.. code:: python

    lmp_alkyl.strucC.tag += '_methane_x50'

.. code:: python

    lmp_alkyl.strucC.write_xyz()

.. code:: python

    for pk,p in lmp_alkyl.strucC.particles.iteritems():
        print p,p.mol,p.residue,p.resname


.. parsed-literal::

    atom[0] C (C) 0 2 OCT
    atom[1] H (H) 0 2 OCT
    atom[2] H (H) 0 2 OCT
    atom[3] H (H) 0 2 OCT
    atom[4] C (C) 0 2 OCT
    atom[5] H (H) 0 2 OCT
    atom[6] H (H) 0 2 OCT
    atom[7] C (C) 0 2 OCT
    atom[8] H (H) 0 2 OCT
    atom[9] H (H) 0 2 OCT
    atom[10] C (C) 0 2 OCT
    atom[11] H (H) 0 2 OCT
    atom[12] H (H) 0 2 OCT
    atom[13] C (C) 0 2 OCT
    atom[14] H (H) 0 2 OCT
    atom[15] C (C) 0 2 OCT
    atom[16] H (H) 0 2 OCT
    atom[17] H (H) 0 2 OCT
    atom[18] C (C) 0 2 OCT
    atom[19] H (H) 0 2 OCT
    atom[20] H (H) 0 2 OCT
    atom[21] C (C) 0 2 OCT
    atom[22] H (H) 0 2 OCT
    atom[23] H (H) 0 2 OCT
    atom[24] H (H) 0 2 OCT
    atom[25] C (C) 0 1 ETH
    atom[26] H (H) 0 1 ETH
    atom[27] H (H) 0 1 ETH
    atom[28] C (C) 0 1 ETH
    atom[29] H (H) 0 1 ETH
    atom[30] H (H) 0 1 ETH
    atom[31] H (H) 0 1 ETH
    atom[32] C (C) 1 2 OCT
    atom[33] H (H) 1 2 OCT
    atom[34] H (H) 1 2 OCT
    atom[35] H (H) 1 2 OCT
    atom[36] C (C) 1 2 OCT
    atom[37] H (H) 1 2 OCT
    atom[38] H (H) 1 2 OCT
    atom[39] C (C) 1 2 OCT
    atom[40] H (H) 1 2 OCT
    atom[41] H (H) 1 2 OCT
    atom[42] C (C) 1 2 OCT
    atom[43] H (H) 1 2 OCT
    atom[44] H (H) 1 2 OCT
    atom[45] C (C) 1 2 OCT
    atom[46] H (H) 1 2 OCT
    atom[47] C (C) 1 2 OCT
    atom[48] H (H) 1 2 OCT
    atom[49] H (H) 1 2 OCT
    atom[50] C (C) 1 2 OCT
    atom[51] H (H) 1 2 OCT
    atom[52] H (H) 1 2 OCT
    atom[53] C (C) 1 2 OCT
    atom[54] H (H) 1 2 OCT
    atom[55] H (H) 1 2 OCT
    atom[56] H (H) 1 2 OCT
    atom[57] C (C) 1 1 ETH
    atom[58] H (H) 1 1 ETH
    atom[59] H (H) 1 1 ETH
    atom[60] C (C) 1 1 ETH
    atom[61] H (H) 1 1 ETH
    atom[62] H (H) 1 1 ETH
    atom[63] H (H) 1 1 ETH
    atom[64] C (C) 2 2 OCT
    atom[65] H (H) 2 2 OCT
    atom[66] H (H) 2 2 OCT
    atom[67] H (H) 2 2 OCT
    atom[68] C (C) 2 2 OCT
    atom[69] H (H) 2 2 OCT
    atom[70] H (H) 2 2 OCT
    atom[71] C (C) 2 2 OCT
    atom[72] H (H) 2 2 OCT
    atom[73] H (H) 2 2 OCT
    atom[74] C (C) 2 2 OCT
    atom[75] H (H) 2 2 OCT
    atom[76] H (H) 2 2 OCT
    atom[77] C (C) 2 2 OCT
    atom[78] H (H) 2 2 OCT
    atom[79] C (C) 2 2 OCT
    atom[80] H (H) 2 2 OCT
    atom[81] H (H) 2 2 OCT
    atom[82] C (C) 2 2 OCT
    atom[83] H (H) 2 2 OCT
    atom[84] H (H) 2 2 OCT
    atom[85] C (C) 2 2 OCT
    atom[86] H (H) 2 2 OCT
    atom[87] H (H) 2 2 OCT
    atom[88] H (H) 2 2 OCT
    atom[89] C (C) 2 1 ETH
    atom[90] H (H) 2 1 ETH
    atom[91] H (H) 2 1 ETH
    atom[92] C (C) 2 1 ETH
    atom[93] H (H) 2 1 ETH
    atom[94] H (H) 2 1 ETH
    atom[95] H (H) 2 1 ETH
    atom[96] C (C) 3 2 OCT
    atom[97] H (H) 3 2 OCT
    atom[98] H (H) 3 2 OCT
    atom[99] H (H) 3 2 OCT
    atom[100] C (C) 3 2 OCT
    atom[101] H (H) 3 2 OCT
    atom[102] H (H) 3 2 OCT
    atom[103] C (C) 3 2 OCT
    atom[104] H (H) 3 2 OCT
    atom[105] H (H) 3 2 OCT
    atom[106] C (C) 3 2 OCT
    atom[107] H (H) 3 2 OCT
    atom[108] H (H) 3 2 OCT
    atom[109] C (C) 3 2 OCT
    atom[110] H (H) 3 2 OCT
    atom[111] C (C) 3 2 OCT
    atom[112] H (H) 3 2 OCT
    atom[113] H (H) 3 2 OCT
    atom[114] C (C) 3 2 OCT
    atom[115] H (H) 3 2 OCT
    atom[116] H (H) 3 2 OCT
    atom[117] C (C) 3 2 OCT
    atom[118] H (H) 3 2 OCT
    atom[119] H (H) 3 2 OCT
    atom[120] H (H) 3 2 OCT
    atom[121] C (C) 3 1 ETH
    atom[122] H (H) 3 1 ETH
    atom[123] H (H) 3 1 ETH
    atom[124] C (C) 3 1 ETH
    atom[125] H (H) 3 1 ETH
    atom[126] H (H) 3 1 ETH
    atom[127] H (H) 3 1 ETH
    atom[128] C (C) 4 2 OCT
    atom[129] H (H) 4 2 OCT
    atom[130] H (H) 4 2 OCT
    atom[131] H (H) 4 2 OCT
    atom[132] C (C) 4 2 OCT
    atom[133] H (H) 4 2 OCT
    atom[134] H (H) 4 2 OCT
    atom[135] C (C) 4 2 OCT
    atom[136] H (H) 4 2 OCT
    atom[137] H (H) 4 2 OCT
    atom[138] C (C) 4 2 OCT
    atom[139] H (H) 4 2 OCT
    atom[140] H (H) 4 2 OCT
    atom[141] C (C) 4 2 OCT
    atom[142] H (H) 4 2 OCT
    atom[143] C (C) 4 2 OCT
    atom[144] H (H) 4 2 OCT
    atom[145] H (H) 4 2 OCT
    atom[146] C (C) 4 2 OCT
    atom[147] H (H) 4 2 OCT
    atom[148] H (H) 4 2 OCT
    atom[149] C (C) 4 2 OCT
    atom[150] H (H) 4 2 OCT
    atom[151] H (H) 4 2 OCT
    atom[152] H (H) 4 2 OCT
    atom[153] C (C) 4 1 ETH
    atom[154] H (H) 4 1 ETH
    atom[155] H (H) 4 1 ETH
    atom[156] C (C) 4 1 ETH
    atom[157] H (H) 4 1 ETH
    atom[158] H (H) 4 1 ETH
    atom[159] H (H) 4 1 ETH
    atom[160] C (C) 5 2 OCT
    atom[161] H (H) 5 2 OCT
    atom[162] H (H) 5 2 OCT
    atom[163] H (H) 5 2 OCT
    atom[164] C (C) 5 2 OCT
    atom[165] H (H) 5 2 OCT
    atom[166] H (H) 5 2 OCT
    atom[167] C (C) 5 2 OCT
    atom[168] H (H) 5 2 OCT
    atom[169] H (H) 5 2 OCT
    atom[170] C (C) 5 2 OCT
    atom[171] H (H) 5 2 OCT
    atom[172] H (H) 5 2 OCT
    atom[173] C (C) 5 2 OCT
    atom[174] H (H) 5 2 OCT
    atom[175] C (C) 5 2 OCT
    atom[176] H (H) 5 2 OCT
    atom[177] H (H) 5 2 OCT
    atom[178] C (C) 5 2 OCT
    atom[179] H (H) 5 2 OCT
    atom[180] H (H) 5 2 OCT
    atom[181] C (C) 5 2 OCT
    atom[182] H (H) 5 2 OCT
    atom[183] H (H) 5 2 OCT
    atom[184] H (H) 5 2 OCT
    atom[185] C (C) 5 1 ETH
    atom[186] H (H) 5 1 ETH
    atom[187] H (H) 5 1 ETH
    atom[188] C (C) 5 1 ETH
    atom[189] H (H) 5 1 ETH
    atom[190] H (H) 5 1 ETH
    atom[191] H (H) 5 1 ETH
    atom[192] C (C) 6 2 OCT
    atom[193] H (H) 6 2 OCT
    atom[194] H (H) 6 2 OCT
    atom[195] H (H) 6 2 OCT
    atom[196] C (C) 6 2 OCT
    atom[197] H (H) 6 2 OCT
    atom[198] H (H) 6 2 OCT
    atom[199] C (C) 6 2 OCT
    atom[200] H (H) 6 2 OCT
    atom[201] H (H) 6 2 OCT
    atom[202] C (C) 6 2 OCT
    atom[203] H (H) 6 2 OCT
    atom[204] H (H) 6 2 OCT
    atom[205] C (C) 6 2 OCT
    atom[206] H (H) 6 2 OCT
    atom[207] C (C) 6 2 OCT
    atom[208] H (H) 6 2 OCT
    atom[209] H (H) 6 2 OCT
    atom[210] C (C) 6 2 OCT
    atom[211] H (H) 6 2 OCT
    atom[212] H (H) 6 2 OCT
    atom[213] C (C) 6 2 OCT
    atom[214] H (H) 6 2 OCT
    atom[215] H (H) 6 2 OCT
    atom[216] H (H) 6 2 OCT
    atom[217] C (C) 6 1 ETH
    atom[218] H (H) 6 1 ETH
    atom[219] H (H) 6 1 ETH
    atom[220] C (C) 6 1 ETH
    atom[221] H (H) 6 1 ETH
    atom[222] H (H) 6 1 ETH
    atom[223] H (H) 6 1 ETH
    atom[224] C (C) 7 2 OCT
    atom[225] H (H) 7 2 OCT
    atom[226] H (H) 7 2 OCT
    atom[227] H (H) 7 2 OCT
    atom[228] C (C) 7 2 OCT
    atom[229] H (H) 7 2 OCT
    atom[230] H (H) 7 2 OCT
    atom[231] C (C) 7 2 OCT
    atom[232] H (H) 7 2 OCT
    atom[233] H (H) 7 2 OCT
    atom[234] C (C) 7 2 OCT
    atom[235] H (H) 7 2 OCT
    atom[236] H (H) 7 2 OCT
    atom[237] C (C) 7 2 OCT
    atom[238] H (H) 7 2 OCT
    atom[239] C (C) 7 2 OCT
    atom[240] H (H) 7 2 OCT
    atom[241] H (H) 7 2 OCT
    atom[242] C (C) 7 2 OCT
    atom[243] H (H) 7 2 OCT
    atom[244] H (H) 7 2 OCT
    atom[245] C (C) 7 2 OCT
    atom[246] H (H) 7 2 OCT
    atom[247] H (H) 7 2 OCT
    atom[248] H (H) 7 2 OCT
    atom[249] C (C) 7 1 ETH
    atom[250] H (H) 7 1 ETH
    atom[251] H (H) 7 1 ETH
    atom[252] C (C) 7 1 ETH
    atom[253] H (H) 7 1 ETH
    atom[254] H (H) 7 1 ETH
    atom[255] H (H) 7 1 ETH
    atom[256] C (C) 8 2 OCT
    atom[257] H (H) 8 2 OCT
    atom[258] H (H) 8 2 OCT
    atom[259] H (H) 8 2 OCT
    atom[260] C (C) 8 2 OCT
    atom[261] H (H) 8 2 OCT
    atom[262] H (H) 8 2 OCT
    atom[263] C (C) 8 2 OCT
    atom[264] H (H) 8 2 OCT
    atom[265] H (H) 8 2 OCT
    atom[266] C (C) 8 2 OCT
    atom[267] H (H) 8 2 OCT
    atom[268] H (H) 8 2 OCT
    atom[269] C (C) 8 2 OCT
    atom[270] H (H) 8 2 OCT
    atom[271] C (C) 8 2 OCT
    atom[272] H (H) 8 2 OCT
    atom[273] H (H) 8 2 OCT
    atom[274] C (C) 8 2 OCT
    atom[275] H (H) 8 2 OCT
    atom[276] H (H) 8 2 OCT
    atom[277] C (C) 8 2 OCT
    atom[278] H (H) 8 2 OCT
    atom[279] H (H) 8 2 OCT
    atom[280] H (H) 8 2 OCT
    atom[281] C (C) 8 1 ETH
    atom[282] H (H) 8 1 ETH
    atom[283] H (H) 8 1 ETH
    atom[284] C (C) 8 1 ETH
    atom[285] H (H) 8 1 ETH
    atom[286] H (H) 8 1 ETH
    atom[287] H (H) 8 1 ETH
    atom[288] C (C) 9 2 OCT
    atom[289] H (H) 9 2 OCT
    atom[290] H (H) 9 2 OCT
    atom[291] H (H) 9 2 OCT
    atom[292] C (C) 9 2 OCT
    atom[293] H (H) 9 2 OCT
    atom[294] H (H) 9 2 OCT
    atom[295] C (C) 9 2 OCT
    atom[296] H (H) 9 2 OCT
    atom[297] H (H) 9 2 OCT
    atom[298] C (C) 9 2 OCT
    atom[299] H (H) 9 2 OCT
    atom[300] H (H) 9 2 OCT
    atom[301] C (C) 9 2 OCT
    atom[302] H (H) 9 2 OCT
    atom[303] C (C) 9 2 OCT
    atom[304] H (H) 9 2 OCT
    atom[305] H (H) 9 2 OCT
    atom[306] C (C) 9 2 OCT
    atom[307] H (H) 9 2 OCT
    atom[308] H (H) 9 2 OCT
    atom[309] C (C) 9 2 OCT
    atom[310] H (H) 9 2 OCT
    atom[311] H (H) 9 2 OCT
    atom[312] H (H) 9 2 OCT
    atom[313] C (C) 9 1 ETH
    atom[314] H (H) 9 1 ETH
    atom[315] H (H) 9 1 ETH
    atom[316] C (C) 9 1 ETH
    atom[317] H (H) 9 1 ETH
    atom[318] H (H) 9 1 ETH
    atom[319] H (H) 9 1 ETH
    atom[320] C (C) 10 1 ETH
    atom[321] H (H) 10 1 ETH
    atom[322] H (H) 10 1 ETH
    atom[323] H (H) 10 1 ETH
    atom[324] C (C) 10 1 ETH
    atom[325] H (H) 10 1 ETH
    atom[326] H (H) 10 1 ETH
    atom[327] H (H) 10 1 ETH
    atom[328] C (C) 11 1 ETH
    atom[329] H (H) 11 1 ETH
    atom[330] H (H) 11 1 ETH
    atom[331] H (H) 11 1 ETH
    atom[332] C (C) 11 1 ETH
    atom[333] H (H) 11 1 ETH
    atom[334] H (H) 11 1 ETH
    atom[335] H (H) 11 1 ETH
    atom[336] C (C) 12 1 ETH
    atom[337] H (H) 12 1 ETH
    atom[338] H (H) 12 1 ETH
    atom[339] H (H) 12 1 ETH
    atom[340] C (C) 12 1 ETH
    atom[341] H (H) 12 1 ETH
    atom[342] H (H) 12 1 ETH
    atom[343] H (H) 12 1 ETH
    atom[344] C (C) 13 1 ETH
    atom[345] H (H) 13 1 ETH
    atom[346] H (H) 13 1 ETH
    atom[347] H (H) 13 1 ETH
    atom[348] C (C) 13 1 ETH
    atom[349] H (H) 13 1 ETH
    atom[350] H (H) 13 1 ETH
    atom[351] H (H) 13 1 ETH
    atom[352] C (C) 14 1 ETH
    atom[353] H (H) 14 1 ETH
    atom[354] H (H) 14 1 ETH
    atom[355] H (H) 14 1 ETH
    atom[356] C (C) 14 1 ETH
    atom[357] H (H) 14 1 ETH
    atom[358] H (H) 14 1 ETH
    atom[359] H (H) 14 1 ETH
    atom[360] C (C) 15 1 ETH
    atom[361] H (H) 15 1 ETH
    atom[362] H (H) 15 1 ETH
    atom[363] H (H) 15 1 ETH
    atom[364] C (C) 15 1 ETH
    atom[365] H (H) 15 1 ETH
    atom[366] H (H) 15 1 ETH
    atom[367] H (H) 15 1 ETH
    atom[368] C (C) 16 1 ETH
    atom[369] H (H) 16 1 ETH
    atom[370] H (H) 16 1 ETH
    atom[371] H (H) 16 1 ETH
    atom[372] C (C) 16 1 ETH
    atom[373] H (H) 16 1 ETH
    atom[374] H (H) 16 1 ETH
    atom[375] H (H) 16 1 ETH
    atom[376] C (C) 17 1 ETH
    atom[377] H (H) 17 1 ETH
    atom[378] H (H) 17 1 ETH
    atom[379] H (H) 17 1 ETH
    atom[380] C (C) 17 1 ETH
    atom[381] H (H) 17 1 ETH
    atom[382] H (H) 17 1 ETH
    atom[383] H (H) 17 1 ETH
    atom[384] C (C) 18 1 ETH
    atom[385] H (H) 18 1 ETH
    atom[386] H (H) 18 1 ETH
    atom[387] H (H) 18 1 ETH
    atom[388] C (C) 18 1 ETH
    atom[389] H (H) 18 1 ETH
    atom[390] H (H) 18 1 ETH
    atom[391] H (H) 18 1 ETH
    atom[392] C (C) 19 1 ETH
    atom[393] H (H) 19 1 ETH
    atom[394] H (H) 19 1 ETH
    atom[395] H (H) 19 1 ETH
    atom[396] C (C) 19 1 ETH
    atom[397] H (H) 19 1 ETH
    atom[398] H (H) 19 1 ETH
    atom[399] H (H) 19 1 ETH
    atom[400] C (C) 20 1 METH
    atom[401] H (H) 20 1 METH
    atom[402] H (H) 20 1 METH
    atom[403] H (H) 20 1 METH
    atom[404] H (H) 20 1 METH
    atom[405] C (C) 21 1 METH
    atom[406] H (H) 21 1 METH
    atom[407] H (H) 21 1 METH
    atom[408] H (H) 21 1 METH
    atom[409] H (H) 21 1 METH
    atom[410] C (C) 22 1 METH
    atom[411] H (H) 22 1 METH
    atom[412] H (H) 22 1 METH
    atom[413] H (H) 22 1 METH
    atom[414] H (H) 22 1 METH
    atom[415] C (C) 23 1 METH
    atom[416] H (H) 23 1 METH
    atom[417] H (H) 23 1 METH
    atom[418] H (H) 23 1 METH
    atom[419] H (H) 23 1 METH
    atom[420] C (C) 24 1 METH
    atom[421] H (H) 24 1 METH
    atom[422] H (H) 24 1 METH
    atom[423] H (H) 24 1 METH
    atom[424] H (H) 24 1 METH
    atom[425] C (C) 25 1 METH
    atom[426] H (H) 25 1 METH
    atom[427] H (H) 25 1 METH
    atom[428] H (H) 25 1 METH
    atom[429] H (H) 25 1 METH
    atom[430] C (C) 26 1 METH
    atom[431] H (H) 26 1 METH
    atom[432] H (H) 26 1 METH
    atom[433] H (H) 26 1 METH
    atom[434] H (H) 26 1 METH
    atom[435] C (C) 27 1 METH
    atom[436] H (H) 27 1 METH
    atom[437] H (H) 27 1 METH
    atom[438] H (H) 27 1 METH
    atom[439] H (H) 27 1 METH
    atom[440] C (C) 28 1 METH
    atom[441] H (H) 28 1 METH
    atom[442] H (H) 28 1 METH
    atom[443] H (H) 28 1 METH
    atom[444] H (H) 28 1 METH
    atom[445] C (C) 29 1 METH
    atom[446] H (H) 29 1 METH
    atom[447] H (H) 29 1 METH
    atom[448] H (H) 29 1 METH
    atom[449] H (H) 29 1 METH
    atom[450] C (C) 30 1 METH
    atom[451] H (H) 30 1 METH
    atom[452] H (H) 30 1 METH
    atom[453] H (H) 30 1 METH
    atom[454] H (H) 30 1 METH
    atom[455] C (C) 31 1 METH
    atom[456] H (H) 31 1 METH
    atom[457] H (H) 31 1 METH
    atom[458] H (H) 31 1 METH
    atom[459] H (H) 31 1 METH
    atom[460] C (C) 32 1 METH
    atom[461] H (H) 32 1 METH
    atom[462] H (H) 32 1 METH
    atom[463] H (H) 32 1 METH
    atom[464] H (H) 32 1 METH
    atom[465] C (C) 33 1 METH
    atom[466] H (H) 33 1 METH
    atom[467] H (H) 33 1 METH
    atom[468] H (H) 33 1 METH
    atom[469] H (H) 33 1 METH
    atom[470] C (C) 34 1 METH
    atom[471] H (H) 34 1 METH
    atom[472] H (H) 34 1 METH
    atom[473] H (H) 34 1 METH
    atom[474] H (H) 34 1 METH
    atom[475] C (C) 35 1 METH
    atom[476] H (H) 35 1 METH
    atom[477] H (H) 35 1 METH
    atom[478] H (H) 35 1 METH
    atom[479] H (H) 35 1 METH
    atom[480] C (C) 36 1 METH
    atom[481] H (H) 36 1 METH
    atom[482] H (H) 36 1 METH
    atom[483] H (H) 36 1 METH
    atom[484] H (H) 36 1 METH
    atom[485] C (C) 37 1 METH
    atom[486] H (H) 37 1 METH
    atom[487] H (H) 37 1 METH
    atom[488] H (H) 37 1 METH
    atom[489] H (H) 37 1 METH
    atom[490] C (C) 38 1 METH
    atom[491] H (H) 38 1 METH
    atom[492] H (H) 38 1 METH
    atom[493] H (H) 38 1 METH
    atom[494] H (H) 38 1 METH
    atom[495] C (C) 39 1 METH
    atom[496] H (H) 39 1 METH
    atom[497] H (H) 39 1 METH
    atom[498] H (H) 39 1 METH
    atom[499] H (H) 39 1 METH
    atom[500] C (C) 40 1 METH
    atom[501] H (H) 40 1 METH
    atom[502] H (H) 40 1 METH
    atom[503] H (H) 40 1 METH
    atom[504] H (H) 40 1 METH
    atom[505] C (C) 41 1 METH
    atom[506] H (H) 41 1 METH
    atom[507] H (H) 41 1 METH
    atom[508] H (H) 41 1 METH
    atom[509] H (H) 41 1 METH
    atom[510] C (C) 42 1 METH
    atom[511] H (H) 42 1 METH
    atom[512] H (H) 42 1 METH
    atom[513] H (H) 42 1 METH
    atom[514] H (H) 42 1 METH
    atom[515] C (C) 43 1 METH
    atom[516] H (H) 43 1 METH
    atom[517] H (H) 43 1 METH
    atom[518] H (H) 43 1 METH
    atom[519] H (H) 43 1 METH
    atom[520] C (C) 44 1 METH
    atom[521] H (H) 44 1 METH
    atom[522] H (H) 44 1 METH
    atom[523] H (H) 44 1 METH
    atom[524] H (H) 44 1 METH
    atom[525] C (C) 45 1 METH
    atom[526] H (H) 45 1 METH
    atom[527] H (H) 45 1 METH
    atom[528] H (H) 45 1 METH
    atom[529] H (H) 45 1 METH
    atom[530] C (C) 46 1 METH
    atom[531] H (H) 46 1 METH
    atom[532] H (H) 46 1 METH
    atom[533] H (H) 46 1 METH
    atom[534] H (H) 46 1 METH
    atom[535] C (C) 47 1 METH
    atom[536] H (H) 47 1 METH
    atom[537] H (H) 47 1 METH
    atom[538] H (H) 47 1 METH
    atom[539] H (H) 47 1 METH
    atom[540] C (C) 48 1 METH
    atom[541] H (H) 48 1 METH
    atom[542] H (H) 48 1 METH
    atom[543] H (H) 48 1 METH
    atom[544] H (H) 48 1 METH
    atom[545] C (C) 49 1 METH
    atom[546] H (H) 49 1 METH
    atom[547] H (H) 49 1 METH
    atom[548] H (H) 49 1 METH
    atom[549] H (H) 49 1 METH
    atom[550] C (C) 50 1 METH
    atom[551] H (H) 50 1 METH
    atom[552] H (H) 50 1 METH
    atom[553] H (H) 50 1 METH
    atom[554] H (H) 50 1 METH
    atom[555] C (C) 51 1 METH
    atom[556] H (H) 51 1 METH
    atom[557] H (H) 51 1 METH
    atom[558] H (H) 51 1 METH
    atom[559] H (H) 51 1 METH
    atom[560] C (C) 52 1 METH
    atom[561] H (H) 52 1 METH
    atom[562] H (H) 52 1 METH
    atom[563] H (H) 52 1 METH
    atom[564] H (H) 52 1 METH
    atom[565] C (C) 53 1 METH
    atom[566] H (H) 53 1 METH
    atom[567] H (H) 53 1 METH
    atom[568] H (H) 53 1 METH
    atom[569] H (H) 53 1 METH
    atom[570] C (C) 54 1 METH
    atom[571] H (H) 54 1 METH
    atom[572] H (H) 54 1 METH
    atom[573] H (H) 54 1 METH
    atom[574] H (H) 54 1 METH
    atom[575] C (C) 55 1 METH
    atom[576] H (H) 55 1 METH
    atom[577] H (H) 55 1 METH
    atom[578] H (H) 55 1 METH
    atom[579] H (H) 55 1 METH
    atom[580] C (C) 56 1 METH
    atom[581] H (H) 56 1 METH
    atom[582] H (H) 56 1 METH
    atom[583] H (H) 56 1 METH
    atom[584] H (H) 56 1 METH
    atom[585] C (C) 57 1 METH
    atom[586] H (H) 57 1 METH
    atom[587] H (H) 57 1 METH
    atom[588] H (H) 57 1 METH
    atom[589] H (H) 57 1 METH
    atom[590] C (C) 58 1 METH
    atom[591] H (H) 58 1 METH
    atom[592] H (H) 58 1 METH
    atom[593] H (H) 58 1 METH
    atom[594] H (H) 58 1 METH
    atom[595] C (C) 59 1 METH
    atom[596] H (H) 59 1 METH
    atom[597] H (H) 59 1 METH
    atom[598] H (H) 59 1 METH
    atom[599] H (H) 59 1 METH
    atom[600] C (C) 60 1 METH
    atom[601] H (H) 60 1 METH
    atom[602] H (H) 60 1 METH
    atom[603] H (H) 60 1 METH
    atom[604] H (H) 60 1 METH
    atom[605] C (C) 61 1 METH
    atom[606] H (H) 61 1 METH
    atom[607] H (H) 61 1 METH
    atom[608] H (H) 61 1 METH
    atom[609] H (H) 61 1 METH
    atom[610] C (C) 62 1 METH
    atom[611] H (H) 62 1 METH
    atom[612] H (H) 62 1 METH
    atom[613] H (H) 62 1 METH
    atom[614] H (H) 62 1 METH
    atom[615] C (C) 63 1 METH
    atom[616] H (H) 63 1 METH
    atom[617] H (H) 63 1 METH
    atom[618] H (H) 63 1 METH
    atom[619] H (H) 63 1 METH
    atom[620] C (C) 64 1 METH
    atom[621] H (H) 64 1 METH
    atom[622] H (H) 64 1 METH
    atom[623] H (H) 64 1 METH
    atom[624] H (H) 64 1 METH
    atom[625] C (C) 65 1 METH
    atom[626] H (H) 65 1 METH
    atom[627] H (H) 65 1 METH
    atom[628] H (H) 65 1 METH
    atom[629] H (H) 65 1 METH
    atom[630] C (C) 66 1 METH
    atom[631] H (H) 66 1 METH
    atom[632] H (H) 66 1 METH
    atom[633] H (H) 66 1 METH
    atom[634] H (H) 66 1 METH
    atom[635] C (C) 67 1 METH
    atom[636] H (H) 67 1 METH
    atom[637] H (H) 67 1 METH
    atom[638] H (H) 67 1 METH
    atom[639] H (H) 67 1 METH
    atom[640] C (C) 68 1 METH
    atom[641] H (H) 68 1 METH
    atom[642] H (H) 68 1 METH
    atom[643] H (H) 68 1 METH
    atom[644] H (H) 68 1 METH
    atom[645] C (C) 69 1 METH
    atom[646] H (H) 69 1 METH
    atom[647] H (H) 69 1 METH
    atom[648] H (H) 69 1 METH
    atom[649] H (H) 69 1 METH


.. code:: python

    lmp_alkyl.set_ffparam()

.. code:: python

    file_type = 'templates'
    file_key = 'in'
    file_name = "lammps_sp.in"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    lmp_alkyl.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    file_type = 'templates'
    file_key = 'run'
    file_name = "lammps.sh"
    from_dirkey = 'templates'
    to_dirkey = 'scratch'
    lmp_alkyl.cp_file(file_type,file_key,file_name,from_dirkey,to_dirkey)

.. code:: python

    os.chdir(lmp_alkyl.dir['scratch'])

.. code:: python

    lmp_alkyl.load_str('templates','in')
    lmp_alkyl.load_str('templates','run')

.. code:: python

    lmp_alkyl.write_data()

.. code:: python

    lmp_alkyl.replacewrite_prop('in','input','in','%s.in'%(lmp_alkyl.tag))

.. code:: python

    lmp_alkyl.properties['input_in'] = lmp_alkyl.files['input']['in']
    lmp_alkyl.replacewrite_prop('run','scripts','run','%s.sh'%(lmp_alkyl.tag))

.. code:: python

    os.chdir(lmp_alkyl.dir['home'])
    lmp_alkyl.dump_json()

.. code:: python

    lmp_alkyl.run()

.. code:: python

    os.chdir(lmp_alkyl.dir['scratch'])
    lmp_alkyl.check()

.. code:: python

    pprint("Calculation:{} has status:{}".format(lmp_alkyl.tag,lmp_alkyl.meta['status']))


.. parsed-literal::

    u'Calculation:lmp_alkyl has status:written'


.. code:: python

    lmp_alkyl.strucC.calc_center_mass()

.. code:: python

    groupset_i = streamm.Groups('mol',lmp_alkyl.strucC)
    groupset_i.group_prop('mol','group_mol')

.. code:: python

    groupset_i.calc_cent_mass()
    groupset_i.calc_radius_asphericity()
    groupset_i.calc_dl()

.. code:: python

    groupset_i.write_cm_xyz()

.. code:: python

    import numpy as np

.. code:: python

    print np.mean(groupset_i.radius),groupset_i.strucC.unit_conf['length']


.. parsed-literal::

    1.79932546227 ang


.. code:: python

    print groupset_i.strucC.lat.pbcs


.. parsed-literal::

    [True, True, True]


.. code:: python

    groupset_i.group_nblist.radii_nblist(groupset_i.strucC.lat,groupset_i.cent_mass,groupset_i.radius,radii_buffer=5.25)

.. code:: python

    groupset_i.group_pbcs()

.. code:: python

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

    group_mol_0 47 0
    [ 38.326923   5.15297    0.311727]
    group_mol_1 35 1
    [ 36.599462  26.819205  37.318824]
    group_mol_2 35 2
    [  3.427126  43.098942  12.819743]
    group_mol_3 34 3
    [  2.385817  34.09571   44.931817]
    group_mol_4 36 4
    [ 12.99275   18.251691  29.674268]
    group_mol_5 45 5
    [ 17.256952   2.611683  24.930946]
    group_mol_6 34 6
    [ 33.107853  35.964402  24.485375]
    group_mol_7 42 7
    [ 18.156043  28.621305  46.606473]
    group_mol_8 39 8
    [ 29.279355  48.702058  37.262101]
    group_mol_9 39 9
    [ 46.718605  28.275142  16.173422]


Fancy aye!
