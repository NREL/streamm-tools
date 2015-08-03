.. _xyz2gromacs:


Convert a structure generated :ref:`Example 3 <molgenex3>` into input files for `GROMACS <http://www.gromacs.org/>`_
===============================================================================================================================================================================================


In this example a we use the ``xyz2gromacs.py`` example script to generate input files for `GROMACS <http://www.gromacs.org/>`_. We use the .xyz file generated in :ref:`Example 3 <molgenex3>` as an input structure. The .xyz file is read in to get the atomic positions and atom types, and a `.itp <http://www.gromacs.org/Documentation/File_Formats/.itp_File>`_ file "conj.itp"  is read in to get a set of reference Force Field parameters. The conj.itp file contains parameters from the OPLSaa force-field that is included in the `GROMACS release <http://www.gromacs.org/Downloads>`_.  And input files for a `GROMACS <http://www.gromacs.org/>`_ run are generated. 

Run::

   python xyz2gromacs.py --in_itp conj.itp 
   --in_com mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.com 
   --out_gro mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.gro 
   --out_itp acc1_thiophene_R_hexane__n5.itp 
   --out_top mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.top


this creates the files::

   mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.gro
   mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.itp 
   mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.top

The `.gro <http://manual.gromacs.org/current/online/gro.html>`_ file includes the structural information, the `.top <http://manual.gromacs.org/current/online/top.html>`_ file is the conectivity file and the new `.itp <http://www.gromacs.org/Documentation/File_Formats/.itp_File>`_ file contains the Force-field parameters for the molecule. 


* :ref:`xyz2data`
* :ref:`replicate_data`
