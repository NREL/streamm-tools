.. _buildingblocks:

buildingblocks
==============


This notebook imports the fundamental objects of the
streamm.buildingblocks module and goes through the functionality of each

Let's start with a methane molecule object

.. code:: python

    import streamm.structures.buildingblock as bb

.. code:: python

    import math

Create a particle object with tag methane

.. code:: python

    mol_i = bb.Buildingblock('methane')

.. code:: python

    print(mol_i)


.. parsed-literal::

     methane


You can read in the .xyz file from the structures example or create a
methane geometry using a molecular viewer such as Avogadro
(https://avogadro.cc/)

.. code:: python

    mol_i.read_xyz()

.. code:: python

    mol_i.bonded_nblist = mol_i.guess_nblist(0,radii_buffer=1.25)

Check that all the particles have been read in

.. code:: python

    print mol_i.n_particles


.. parsed-literal::

    5


Check that the neighbor list was set correctly

.. code:: python

    print mol_i.bonded_nblist


.. parsed-literal::

     NBlist of 5 particle with 8 connections


Looks good, you should have the geometry of a methane molecule with a
C-H bond length of 1.2 Angstroms

We want to use the functionality of the buildingblock object to join two
methane's together to create alkyl chains of any length

So let's set two of the hydrogens to be reactive sites (rsites).

You can view the numerical order of the atoms in Avogadro by setting the
label to "atom number," however, avogadro labels for 1-N, while stream
uses 0 - N-1

We will choose the first two hydrogens and set their rsite variable to
'RH'. It does not matter what this identifier is, as long as the same
identifier is passed to the attach() function later. Also, if the
identifiers are not unique, the order in which it appears in the
particles list will also be used.

.. code:: python

    mol_i.particles[1].rsite = 'RH'

.. code:: python

    mol_i.particles[2].rsite = 'RH'

Now use the find\_rsites() function to create the dictionary of lists to
be used by the attach() function

.. code:: python

    mol_i.find_rsites()

.. code:: python

    print mol_i.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom[1] H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom[2] H (H) index:2 n_bonds:1] 
    


Pass the molecule to the attach function and set the rsite id's and the
list positions of the rsites

.. code:: python

    mol_j = bb.attach(mol_i,mol_i,'RH',0,'RH',1,tag='ethane')

Write the .xyz to file to be viewed with a molecular viewer.

.. code:: python

    mol_j.write_xyz()

While the ethane molecule was generated, the hydrogens are eclipsed
rather than staggered.

We can avoid this by using the prepattach() function to orient the
molecule and remove the reactive site

.. code:: python

    mol_k = mol_i.prepattach('RH',0,dir=-1,yangle=90.0)

Then apply a shift to set the bond length

.. code:: python

    CC_bl = mol_i.particles[0].bonded_radius*2.0
    mol_k.shift_pos([CC_bl,0.0,0.0])

Then apply a rotation to set the conformation to staggered. Use a 180.0
degree rotation to place the reactive site in the correct orientation
for subsequent attachments.

.. code:: python

    angle_rad = 180.0*math.pi/180.0 
    mol_k.rotate_yz(angle_rad)

.. code:: python

    mol_l = mol_i.prepattach('RH',1,dir=1)

.. code:: python

    mol_m = bb.attachprep(mol_k,mol_l)

.. code:: python

    mol_m.tag = 'ethane'

.. code:: python

    mol_m.write_xyz()

.. code:: python

    print mol_m.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom[1] H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom[5] H (H) index:5 n_bonds:1] 
    


Attachments can also be done in a loop

.. code:: python

    alkly_n = (12-1)/2 # Number of ethanes to add to get a dodecyl 

.. code:: python

    print alkly_n


.. parsed-literal::

    5


.. code:: python

    mol_n = mol_m 

.. code:: python

    mol_n.find_rsites()

.. code:: python

    print mol_n.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom[1] H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom[5] H (H) index:5 n_bonds:1] 
    


.. code:: python

    for i in range(alkly_n):
        mol_n = bb.attach(mol_n,mol_m,'RH',1,'RH',0)

.. code:: python

    mol_n.tag = 'dodecyl'

.. code:: python

    mol_n.write_xyz()

Oh, so alkyl!
