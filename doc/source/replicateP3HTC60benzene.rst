.. _replicateP3HTC60benzene:


Replicate a P3HT and C60 molecule with benzene as a solvent for an MD run
--------------------------------------------------------------------------------------------------------------

In this example a we use the replicate_multimol.py script to read in
the structure from a `LAMMPS <http://lammps.sandia.gov/>`_ input files of a  :ref:`P3HT <P3HT>`,  :ref:`C60 <C60>` and :ref:`benzene  <benzene>`.

Run::

       ./replicate_multimol.py   -v 
       --mol1_data  mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data  
       --mol1_n 10 
       --mol2_data mols/C60/acc1_C60_n1.data 
       --mol2_n 10 
       --mol3_data mols/benzene/acc1_benzene_n1.data    
       --mol3_n 50
       --out_data  mols/P3HTC60benzene_mix1.data    
       --out_xyz mols/P3HTC60benzene_mix1.xyz 

This reads in the structures of P3HT (thiophene_R_hexane) with five
repeat units and replicates it randomly in space 10 times. Similarly
for C60 is replicated randomly in space 10 times to make a 1:1 mixed
phase. Benzene is treated as a solvent and is replicated on a grid for
faster replicatation. For each benzene added the distance between the
benzene and the already added P3HT and C60 molecules is calculated,
and if benzene is within the cutoff it is not added to the grid
point. The new mixed phase `LAMMPS <http://lammps.sandia.gov/>`_ input file::

    mols/P3HTC60benzene_mix1.data   
    mols/P3HTC60benzene_mix1.xyz 

.. note::


   This entire example can be executed by running::

      replicateP3HTC60benzene.sh

   in `tools/examples`. Needed files and scripts are 

   *  replicate_multimol.py

   This example requires output from previous examples

   * :ref:`P3HT <P3HT>` 
   * :ref:`C60 <C60>` 
   * :ref:`benzene  <benzene>` 


