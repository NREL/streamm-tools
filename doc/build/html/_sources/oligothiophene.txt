.. _oligothiophene:


This entire example can be executed by running::

   oligothiophene.sh

in `tools/examples`. Needed files and scripts are 

*  donoracceptor.com.template
*  donoracceptor.pbs.template 
*  conj.itp      
*  xyz2data.py
*  xyz2gromacs.py

Oligo-thiophene
-------------------------------------------------------

*Generate structure and quantum input*

Generate a .xyz file, a `Gaussian <http://www.gaussian.com/>`_ .com input file and a submission
script ".pbs"  for thiophene by running ::

   donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release  -r 5

will generate oligomer of thiophene with n=1-5 repeat units::

   mols/thiophene/acc1_thiophene_n1.xyz
   mols/thiophene/acc1_thiophene_n2.xyz
   mols/thiophene/acc1_thiophene_n3.xyz
   mols/thiophene/acc1_thiophene_n4.xyz
   mols/thiophene/acc1_thiophene_n5.xyz

you will notice that the angles between repeat units is rather
random. In order to control the angle between units use the -p
function::

   donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release  
      -r 5 -p "0"

sets the inter-ring dihedral angle to zero making the sulfurs of
thiophene in the cis configuration ::

   donoracceptorsystems.py  "thiophene" -b BuildingBlocks-release  
      -r 5  -p "180 0 "

sets the inter-ring dihedral angle to alternate between zero and 180
making the sulfurs of thiophene in the trans configuration. The
`Gaussian <http://www.gaussian.com/>`_  input files are also created for each oligomer::

   mols/thiophene/acc1_thiophene_n1.com
   mols/thiophene/acc1_thiophene_n1.pbs
   ...
   mols/thiophene/acc1_thiophene_n5.com
   mols/thiophene/acc1_thiophene_n5.pbs
 


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
    --in_xyz  mols/thiophene/acc1_thiophene_n5.xyz 
    --out_data  mols/thiophene/acc1_thiophene_n5.data

this creates the files::

    mols/thiophene/acc1_thiophene_n5.data

And an input `.data
<http://lammps.sandia.gov/doc/2001/data_format.html>`_  for a
`LAMMPS <http://lammps.sandia.gov/>`_ run is generated. 


And input files for a `GROMACS <http://www.gromacs.org/>`_ run are
generated `xyz2gromacs.py`::

   xyz2gromacs.py --in_itp conj.itp 
      --in_xyz  mols/thiophene/acc1_thiophene_n5.xyz 
      --out_gro mols/thiophene/acc1_thiophene_n5.gro 
      --out_top mols/thiophene/acc1_thiophene_n5.top
      --out_itp  acc1_thiophene_n5.itp 

this creates the files::

      mols/thiophene/acc1_thiophene_n5.gro 
      mols/thiophene/acc1_thiophene_n5.top
      acc1_thiophene_n5.itp 

The `.gro <http://manual.gromacs.org/current/online/gro.html>`_ file includes the structural information, the `.top <http://manual.gromacs.org/current/online/top.html>`_ file is the connectivity file and the new `.itp <http://www.gromacs.org/Documentation/File_Formats/.itp_File>`_ file contains the Force-field parameters for the molecule. 

