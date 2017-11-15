.. _buildingblocks_example:
  
buildingblocks_example
===============
 

This notebook imports the fundamental objects of the
streamm.buildingblocks module and goes through the functionality of each

.. code:: ipython2

    import streamm.structures.buildingblock as bb

.. code:: ipython2

    import math

.. code:: ipython2

    from pathlib2 import Path
    import os

Create a Buildingblock object with tag methane

.. code:: ipython2

    mol_i = bb.Buildingblock('methane')

.. code:: ipython2

    print(mol_i.print_properties())

If methane.xyz is not around run the structures example

.. code:: ipython2

    need_files = ['methane.xyz']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run structures_example.ipynb")
            os.system("jupyter nbconvert --to python  structures_example.ipynb")
            os.system("python structures_example.py")

.. code:: ipython2

    mol_i.read_xyz()

.. code:: ipython2

    mol_i.bonded_nblist = mol_i.guess_nblist(0,radii_buffer=1.25)

Check that all the particles have been read in

.. code:: ipython2

    print mol_i.n_particles

Check that the neighbor list was set correctly

.. code:: ipython2

    print mol_i.bonded_nblist

Looks good, you should have the geometry of a methane molecule with a
C-H bond length of 1.2 Angstroms

We want to use the functionality of the buildingblock object to join two
methanes together to create alkyl chains of any length

So let’s set two of the hydrogens to be reactive sites (rsites).

You can view the numerical order of the atoms in Avogadro by setting the
label to “atom number,” however, Avogadro labels atoms from 1 to N,
while streamm uses 0 to N-1

We will choose the first two hydrogens and set their rsite variable to
‘RH’. It does not matter what this identifier is, as long as the same
identifier is passed to the attach() function later. Also, if the
identifiers are not unique, the order in which it appears in the
particles list will also be used.

.. code:: ipython2

    mol_i.particles[1].rsite = 'RH'

.. code:: ipython2

    mol_i.particles[2].rsite = 'RH'

Now use the find_rsites() function to create the dictionary of lists to
be used by the attach() function

.. code:: ipython2

    mol_i.find_rsites()

.. code:: ipython2

    print mol_i.show_rsites()

Pass the molecule to the attach function and set the rsite id’s and the
list positions of the rsites

.. code:: ipython2

    mol_j = bb.attach(mol_i,mol_i,'RH',0,'RH',1,tag='ethane')

Write the .xyz to file to be viewed with a molecular viewer.

.. code:: ipython2

    mol_j.write_xyz()

While the ethane molecule was generated, the hydrogens are eclipsed
rather than staggered.

We can avoid this by using the prepattach() function to orient the
molecule and remove the reactive site

.. code:: ipython2

    mol_k = mol_i.prepattach('RH',0,dir=-1,yangle=90.0)

Then apply a shift to set the bond length

.. code:: ipython2

    CC_bl = mol_i.particles[0].bonded_radius*2.0
    mol_k.shift_pos([CC_bl,0.0,0.0])

Then apply a rotation to set the conformation to staggered. Use a 180.0
degree rotation to place the reactive site in the correct orientation
for subsequent attachments.

.. code:: ipython2

    angle_rad = 180.0*math.pi/180.0 
    mol_k.rotate_yz(angle_rad)

.. code:: ipython2

    mol_l = mol_i.prepattach('RH',1,dir=1)

.. code:: ipython2

    mol_m = bb.attachprep(mol_k,mol_l)

.. code:: ipython2

    mol_m.tag = 'ethane'

.. code:: ipython2

    for pk,p in mol_m.particles.iteritems():
        print pk,p

.. code:: ipython2

    print mol_m.bonded_nblist.list 
    print mol_m.bonded_nblist.index 

.. code:: ipython2

    mol_m.write_xyz()

.. code:: ipython2

    print mol_m.show_rsites()

.. code:: ipython2

    mol_m.bonded_bonds()
    mol_m.bonded_angles()
    mol_m.bonded_dih()

.. code:: ipython2

    mol_json = mol_m.export_json()

Attachments can also be done in a loop

.. code:: ipython2

    alkly_n = (12-1)/2 # Number of ethanes to add to get a dodecyl 

.. code:: ipython2

    print alkly_n

.. code:: ipython2

    mol_n = mol_m 

.. code:: ipython2

    mol_n.find_rsites()

.. code:: ipython2

    print mol_n.show_rsites()

.. code:: ipython2

    for i in range(alkly_n):
        mol_n = bb.attach(mol_n,mol_m,'RH',1,'RH',0)

.. code:: ipython2

    mol_n.tag = 'dodecyl'

.. code:: ipython2

    mol_n.write_xyz()

Oh, so alkyl!
