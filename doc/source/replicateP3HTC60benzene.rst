.. _replicateP3HTC60benzene:

This entire example can be executed by running::

   replicateP3HTC60benzene.sh

in `tools/examples`. Needed files and scripts are 

*  replicate_multimol.py

This example requires output from previous examples

* :ref:`P3HT <P3HT>` 
* :ref:`C60 <C60>` 
* :ref:`benzene  <benzene>` 


Replicate a P3HT and C60 molecule with benzene as a solvent for an MD run
--------------------------------------------------------------------------------------------------------------

TRAVIS: fix all nearby text
In this example a we use the replicate_multimol.py script to read in
the structure from a `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_ files of a P3HT,  C60 and benzene.

Run::

       ./replicate_multimol.py   -v 
       --mol1_data  mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data  
       --mol1_n 10 
       --mol2_data mols/C60/acc1_C60_n1.data --mol2_n 10 
       --mol3_data mols/benzene/acc1_benzene_n1.data    --mol3_n 50
       --out_data  mols/P3HTC60benzene_mix1.data    
       --out_xyz mols/P3HTC60benzene_mix1.xyz 

this creates the files::

    mols/P3HTC60benzene_mix1.data   
    mols/P3HTC60benzene_mix1.xyz 

The new `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_  file has 10
P3HT molecules and 10 C60's placed randomly in a simulation cell and
50 benzene molecules placed on a grid. 
