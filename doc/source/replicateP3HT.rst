.. _replicateP3HT:


This entire example can be executed by running::

   replicateP3HT.sh

in tools/examples


Replicate a P3HT for an MD run
-------------------------------------------------------

In this example a we use the ``replicate_data.py`` script to read in
the structure from a `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_ file and
replicate the structure in space in two ways. 

The first method is the default and chooses a random point within the
box dimensions specified in the `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_  file and places a copied
structure at that point after randomly rotating the structure along
the 2 axes.  Then distance between that added particles of the
copied structure and all the particles in the system is calculated
to be sure no 2 structures overlap each other. 

This example requires output from previous examples

* :ref:`P3HT <P3HT>` 

Run::

   replicate_data.py --mol_n 10 
   --in_data mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data 
   --out_data  mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.data 
   --out_xyz mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.xyz 


this creates the files::

    mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.data
    mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.xyz

The new `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_  file with 10
molecules randomly placed and randomly oriented. 

The second method is to place the molecular structures on a grid. 

Run::

    replicate_data.py --mol_n 10 --grid
    --in_data mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data 
    --out_data mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10_g.data 
    --out_xyz mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10_g.xyz 


this creates the files::

    mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10_g.data
    mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10_g.xyz

The new `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_  file with 10
molecules placed on grid in the box defined by the original `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_ file. The grid
method has the advantage on not need to calculate the inter-particle
distances, so it is much faster.  
