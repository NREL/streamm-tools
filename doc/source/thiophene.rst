.. _thiophene:


This entire example can be exicuted by running::

   thiophene.sh

in tools/examples


Thiophene
-------------------------------------------------------

*Generate structure and quantum input*

Generate a .xyz file, a `Gaussian <http://www.gaussian.com/>`_  .com input file and a submission
script ".pbs"  for thiophene by running ::

   donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release  -r 1 

this creates the file::

   mols/thiophene/acc1_thiophene_n1.xyz

which is just a file containing the cartesian coordinates of a
thiophene molecule.  You can view with your favorite viewer. The -r option is set to 1 to generate a single molecule rather than an oligomer::

   mols/thiophene/acc1_thiophene_n1.com

which is the `Gaussian <http://www.gaussian.com/>`_  input file based on "donoracceptor.com.template"::

   mols/thiophene/acc1_thiophene_n1.pbs

which is the pbs script file based on `donoracceptor.pbs.template`


*Generate topology  files*

The .xyz file is read in to get the atomic positions and
atom types, and a `.itp
<http://www.gromacs.org/Documentation/File_Formats/.itp_File>`_ file
"conj.itp"  is read in to get a set of reference Force Field
parameters. The conj.itp file contains parameters from the OPLSaa
force-field that is included in the `GROMACS release
<http://www.gromacs.org/Downloads>`_.  

For `LAMMPS <http://lammps.sandia.gov/>`_ input generation run `xyz2data.py`::

  xyz2data.py  --in_itp conj.itp 
    --in_xyz  mols/thiophene/acc1_thiophene_n1.xyz 
    --out_data  mols/thiophene/acc1_thiophene_n1.data

this creates the files::

    mols/thiophene/acc1_thiophene_n1.data

And an input `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_  for a
`LAMMPS <http://lammps.sandia.gov/>`_ run is generated. 


And input files for a `GROMACS <http://www.gromacs.org/>`_ run are
generated `xyz2gromacs.py`::

   xyz2gromacs.py --in_itp conj.itp 
      --in_xyz  mols/thiophene/acc1_thiophene_n1.xyz 
      --out_gro mols/thiophene/acc1_thiophene_n1.gro 
      --out_top mols/thiophene/acc1_thiophene_n1.top
      --out_itp  acc1_thiophene_n1.itp 

this creates the files::

      mols/thiophene/acc1_thiophene_n1.gro 
      mols/thiophene/acc1_thiophene_n1.top
      acc1_thiophene_n1.itp 

The `.gro <http://manual.gromacs.org/current/online/gro.html>`_ file includes the structural information, the `.top <http://manual.gromacs.org/current/online/top.html>`_ file is the conectivity file and the new `.itp <http://www.gromacs.org/Documentation/File_Formats/.itp_File>`_ file contains the Force-field parameters for the molecule. 

