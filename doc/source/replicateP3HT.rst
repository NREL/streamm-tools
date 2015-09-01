.. _replicateP3HT:


Replicate a P3HT for an MD run
-------------------------------------------------------

In this example a we use the ``replicate_data.py`` script to read in
the structure from a `LAMMPS <http://lammps.sandia.gov/>`_ input file and
replicate the structure in space in two ways. 

The first method is the default and chooses a random point within the
box dimensions specified in the `LAMMPS <http://lammps.sandia.gov/>`_ input file and places a copied
structure at that point after randomly rotating the structure along
2 axes.  Then the distance between the added particles of the
copied structure and all the particles in the system is calculated
to be sure structures are not overlapping each other. 


Run::

   replicate_data.py --mol_n 10 
   --in_data mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data 
   --out_data  mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.data 
   --out_xyz mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.xyz 


this creates the files::

    mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.data
    mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5x10.xyz

The new `LAMMPS <http://lammps.sandia.gov/>`_ input file with 10
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

The new `LAMMPS <http://lammps.sandia.gov/>`_ input  file with 10
molecules placed on grid in the box defined by the original `LAMMPS <http://lammps.sandia.gov/>`_ input  file. The grid method has the
advantage on not need to calculate the inter-particle distances, so it is much faster.  

.. note::


   This entire example can be executed by running::

      replicateP3HT.sh
   
   in `tools/examples`. Needed files and scripts are 

   *  replicate_data.py

   This example requires output from previous examples

   * :ref:`P3HT <P3HT>` 
