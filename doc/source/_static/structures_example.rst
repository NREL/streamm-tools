.. _structures_example:
  
structures_example
========================
 

This notebook imports the fundamental objects of the streamm.structure
module and goes through the functionality of each

.. code:: python

    from pprint import pprint 
    import copy

Set up a log file for this example so we can read what exactly streamm
is doing, if we feel like it.

.. code:: python

    import logging
    logging.basicConfig(filename='structures_example.log',level=logging.INFO)

Let’s start with the Particle object

.. code:: python

    from streamm.structures.particle import Particle

Create a particle object with label ‘C1’

.. code:: python

    p_i = Particle(label='C1')

.. code:: python

    print(p_i)


.. parsed-literal::

    atom C1 (C1)


Assign the carbon element to the particle

.. code:: python

    p_i.set_element('C')

Let’s oxidize the carbon just to make the charge non-zero

.. code:: python

    p_i.charge = -1.0

Check that the element properties were set to the particle

.. code:: python

    print p_i.show_attributes()


.. parsed-literal::

     type:atom 
     label:C1
     symbol:C
     mass:12.0107 (amu)
     charge:-1.0 (e)
     bonded_radius:0.67 (ang)
     nonbonded_radius:1.7 (ang)


Say we want to change the units to SI

Let’s look at the current units of the particle instance

.. code:: python

    default_unit_conf = copy.deepcopy(p_i.unit_conf)
    pprint(default_unit_conf)


.. parsed-literal::

    {u'amount': u'atom',
     u'angle': u'degree',
     u'capacitance': u'F',
     u'charge': u'e',
     u'conductance': u'S',
     u'current': u'A',
     u'density': u'amu_nm^3',
     u'electric_dipole_moment': u'D',
     u'emf': u'V',
     u'energy': u'Ha',
     u'force': u'GN',
     u'frequency': u'Hz',
     u'harm_bond_coeff': u'kCalmolsqang',
     u'intensity': u'cd',
     u'length': u'ang',
     u'magnetic_flux': u'Wb',
     u'mass': u'amu',
     u'memory': u'Kb',
     u'power': u'GW',
     u'pressure': u'KPa',
     u'resistance': u'ohm',
     u'temperature': u'K',
     u'time': u'ns',
     u'volume': u'nm^3'}


Create a dictionary with new units

.. code:: python

    new_unit_conf = {'length':'m','mass':'kg','charge':'C'}

.. code:: python

    p_i.update_units(new_unit_conf)

.. code:: python

    print p_i.show_attributes()


.. parsed-literal::

     type:atom 
     label:C1
     symbol:C
     mass:1.99442362477e-26 (kg)
     charge:-1.6021766208e-19 (C)
     bonded_radius:6.7e-11 (m)
     nonbonded_radius:1.7e-10 (m)


That’s cool, but we should stick with the default units, so let’s change
them back

.. code:: python

    p_i.update_units(default_unit_conf)

.. code:: python

    print p_i.show_attributes()


.. parsed-literal::

     type:atom 
     label:C1
     symbol:C
     mass:12.0107 (amu)
     charge:-1.0 (e)
     bonded_radius:0.67 (ang)
     nonbonded_radius:1.7 (ang)


Let’s create another particle and set the element to hydrogen

.. code:: python

    p_j = Particle(symbol='H')

.. code:: python

    print p_j.show_attributes()


.. parsed-literal::

     type:atom 
     label:H
     symbol:H
     mass:1.00794 (amu)
     charge:0.0 (e)
     bonded_radius:0.53 (ang)
     nonbonded_radius:1.2 (ang)


Let’s make an empty structure container

.. code:: python

    from streamm.structures.structure import Structure

.. code:: python

    mol_i = Structure('methane')

Now let’s construct a molecule

We can add the carbon at the origin using the ``add_partpos()``
function.

.. code:: python

    pos_i = [0.0,0.0,0.0]
    mol_i.add_partpos(p_i,pos_i)

.. code:: python

    for p_index,particle_i in mol_i.particles.iteritems():
        if( particle_i.symbol == 'H' ):
            particle_i.residue = 1
    
            h_cnt += 1
            

.. code:: python

    for p_index,particle_i in mol_i.particles.iteritems():
        print p_index,particle_i


.. parsed-literal::

    0 atom C1 (C)


.. code:: python

    print("Now the structure container has {} particle ".format(mol_i.n_particles))


.. parsed-literal::

    Now the structure container has 1 particle 


Find the positions of the hydrogens to give a tetrahedral molecular
geometry

.. code:: python

    import numpy as np
    import decimal

.. code:: python

    bond_length = p_i.bonded_radius + p_j.bonded_radius 

.. code:: python

    print bond_length,mol_i.unit_conf['length']


.. parsed-literal::

    1.2 ang


.. code:: python

    tet_a = bond_length/np.sqrt(3)

.. code:: python

    print tet_a


.. parsed-literal::

    0.692820323028


Add hydrogens

.. code:: python

    pos_j = [tet_a,tet_a,tet_a]
    mol_i.add_partpos(p_j,pos_j)

.. code:: python

    for p_index,particle_i in mol_i.particles.iteritems():
        print p_index,particle_i


.. parsed-literal::

    0 atom C1 (C)
    1 atom H (H)


We can add the subsequent hydrogens using the same particle object since
add_partpos makes a deepcopy of the object when adding to the structure
container

.. code:: python

    pos_j = [-tet_a,-tet_a,tet_a]
    mol_i.add_partpos(p_j,pos_j)

.. code:: python

    pos_j = [-tet_a,tet_a,-tet_a]
    mol_i.add_partpos(p_j,pos_j)

.. code:: python

    pos_j = [tet_a,-tet_a,-tet_a]
    mol_i.add_partpos(p_j,pos_j)

Check the position array

.. code:: python

    print mol_i.positions


.. parsed-literal::

    [[ 0.          0.          0.        ]
     [ 0.69282032  0.69282032  0.69282032]
     [-0.69282032 -0.69282032  0.69282032]
     [-0.69282032  0.69282032 -0.69282032]
     [ 0.69282032 -0.69282032 -0.69282032]]


The particles instance variable of the structure container is a
dictionary, so we can just loop over that using the iteritems()
function.

.. code:: python

    for p_index,particle_i in mol_i.particles.iteritems():
        print p_index,particle_i


.. parsed-literal::

    0 atom C1 (C)
    1 atom H (H)
    2 atom H (H)
    3 atom H (H)
    4 atom H (H)


Hum, let’s fix the labels of the hydrogens…

.. code:: python

    h_cnt = 1
    for p_index,particle_i in mol_i.particles.iteritems():
        if( particle_i.symbol == 'H' ):
            particle_i.label = 'H{}'.format(h_cnt)
    
            h_cnt += 1
            

.. code:: python

    for p_index,particle_i in mol_i.particles.iteritems():
        print p_index,particle_i 


.. parsed-literal::

    0 atom C1 (C)
    1 atom H1 (H)
    2 atom H2 (H)
    3 atom H3 (H)
    4 atom H4 (H)


Okay, that looks better

Print .xyz file and check geometry with a molecular viewer such as
Avogadro (https://avogadro.cc/)

.. code:: python

    mol_i.write_xyz()

Looks good, you should have the geometry of a methane molecule with a
C-H bond length of 1.2 Angstroms

However, we have not told streamm about the bonds. There are a few ways
to do this, let’s do it explicitly with the Bond object fist.

.. code:: python

    from streamm.structures.bond import Bond

based on the particle index values

.. code:: python

    b_ij = Bond(0,1)

Now add the bond to the bonds dictionary in the structure container

.. code:: python

    mol_i.add_bond(b_ij)

.. code:: python

    print("Now the structure container has {} particle/s and {} bond/s".format(mol_i.n_particles,mol_i.n_bonds))


.. parsed-literal::

    Now the structure container has 5 particle/s and 1 bond/s


Neat, but adding all the bonds, bond angles and dihedrals explicitly
would be pretty tedious, so let’s use some functions to do that.

First, let’s guess the ``bonded_nblist`` of the molecule based on the
``bonded_radius`` of each particle (atom)

.. code:: python

    mol_i.bonded_nblist = mol_i.guess_nblist(0,radii_buffer=1.35)

.. code:: python

    print mol_i.bonded_nblist


.. parsed-literal::

     NBlist of 5 particles with 8 connections


Let’s take a look at the neighbor lists ``list`` and ``index`` instance
variables

.. code:: python

    print mol_i.bonded_nblist.list 
    print mol_i.bonded_nblist.index 


.. parsed-literal::

    [1, 2, 3, 4, 0, 0, 0, 0]
    [0, 4, 5, 6, 7, 8]


Looking at the ``index`` for particle 0, we get that it has neighbors in
the ``list`` from 0:3 (index[0]:index[0+1]-1). Therefore we know
particle 0 has [1, 2, 3, 4] for neighbors.

.. code:: python

    print mol_i.bonded_nblist.calc_nnab(0)


.. parsed-literal::

    4


Now we can use the bonded neighbor list to construct the bonds, bond
angles and dihedrals

.. code:: python

    mol_i.bonded_bonds()
    mol_i.bonded_angles()
    mol_i.bonded_dih()

.. code:: python

    print(mol_i.print_properties())


.. parsed-literal::

     n_particles:5 
     n_bonds:4
     n_angles:6
     n_dihedrals:0
     n_impropers:0


A little easier than adding everything by hand

We can dump a json file ``methane_struc.json`` for checkpointing the
entire structure

.. code:: python

    mol_json = mol_i.export_json()

.. code:: python

    print mol_i.tag


.. parsed-literal::

    methane


We can reload the molecule as another instance of the Structure object

.. code:: python

    mol_test = Structure('methane')

.. code:: python

    mol_test.import_json()

Then we can use the structure’s **eq** and **ne** functions to make sure
they are the same

.. code:: python

    if( mol_test != mol_i ):
        print "Molecules not equal"
    elif( mol_test == mol_i ):
        print "Molecules equal"
    else:
        print "Uh Oh"


.. parsed-literal::

    Molecules equal


Now let’s set some groups. This is a little unnecessary for methane, but
it will come in super handy if you have a large simulation of thousands
of molecules.

To do this we will set the residue variable for each particle.

.. code:: python

    mol_i.particles[0].residue = 0
    for p_index,particle_i in mol_i.particles.iteritems():
        if( particle_i.symbol == 'H' ):
            particle_i.residue = 1
        print particle_i, particle_i.residue


.. parsed-literal::

    atom C1 (C) 0
    atom H1 (H) 1
    atom H2 (H) 1
    atom H3 (H) 1
    atom H4 (H) 1


.. code:: python

    import streamm.structures.group as group

.. code:: python

    groups_i = group.Groups('methane_residues',mol_i)

Find groups based on residue variable

.. code:: python

    groups_i.group_prop('residue',groups_i.tag)

.. code:: python

    for g_index,group_i in groups_i.groups.iteritems():
        print group_i.pkeys
        print group_i.write_coord()


.. parsed-literal::

    [0]
         C       0.00000000       0.00000000       0.00000000 
    
    [1, 2, 3, 4]
         H       0.69282032       0.69282032       0.69282032 
         H      -0.69282032      -0.69282032       0.69282032 
         H      -0.69282032       0.69282032      -0.69282032 
         H       0.69282032      -0.69282032      -0.69282032 
    


Looks good. We have two groups in the group container, the first with
the carbon particle index 0 and the rest are the hyrdogens.

.. code:: python

    res_json = groups_i.export_json()

Now let’s change the units

.. code:: python

    mol_i.update_units({'length':'pm'})

Check the positions

.. code:: python

    print mol_i.positions


.. parsed-literal::

    [[  0.          0.          0.       ]
     [ 69.2820323  69.2820323  69.2820323]
     [-69.2820323 -69.2820323  69.2820323]
     [-69.2820323  69.2820323 -69.2820323]
     [ 69.2820323 -69.2820323 -69.2820323]]


Check the particle bond radii

.. code:: python

    for p_index,particle_i in mol_i.particles.iteritems():
        print particle_i,particle_i.bonded_radius


.. parsed-literal::

    atom C1 (C) 67.0
    atom H1 (H) 53.0
    atom H2 (H) 53.0
    atom H3 (H) 53.0
    atom H4 (H) 53.0


Cool beans bro!
