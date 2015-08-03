.. _xyz2data:


Convert a structure generated :ref:`Example 3 <molgenex3>` into input files for `LAMMPS <http://lammps.sandia.gov/>`_
================================================================================================================================================

In this example a we use the ``xyz2data.py`` example script to
generate input files for `LAMMPS <http://lammps.sandia.gov/>`_. We use
the .xyz file generated in :ref:`Example 3 <molgenex3>` as an input
structure. The .xyz file is read in to get the atomic positions and
atom types, and a `.itp
<http://www.gromacs.org/Documentation/File_Formats/.itp_File>`_ file
"conj.itp"  is read in to get a set of reference Force Field
parameters. The conj.itp file contains parameters from the OPLSaa
force-field that is included in the `GROMACS release
<http://www.gromacs.org/Downloads>`_.  And an input `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_  for a
`LAMMPS <http://lammps.sandia.gov/>`_ run is generated. 


Run::

  python xyz2data.py --in_itp conj.itp 
  --in_com mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.com 
  --out_data mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data
 

this creates the files::

    mols/thiophene_R_hexane_/acc1_thiophene_R_hexane__n5.data

* :ref:`replicate_data`
