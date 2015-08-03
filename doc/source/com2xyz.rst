.. _com2xyz:


Convert `Gaussian <http://www.gaussian.com/>`_  input generated :ref:`Example 3 <molgenex3>` into xyz file 
================================================================================

So this first example is a bit redundant, in that there is already an
.xyz file generated using the :ref:`donoracceptorsystems.py<donoracceptorsystems>`
scripts, in :ref:`Example 3 <molgenex3>`. However, it is often useful to generate .xyz files to view
structures. 

Run::

    python com2xyz.py 
    --in_com mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.com 
    --out_xyz mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.xyz


this creates the file::

   mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.xyz

which is just a file containing the cartesian coordinates of the P3HT molecule.  You can view with your favorite viewer. 

* :ref:`xyz2gromacs`
* :ref:`xyz2data`
* :ref:`replicate_data`

