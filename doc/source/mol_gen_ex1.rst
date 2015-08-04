.. _molgenex1:


Example 1: Gaussian input file for thiophene
================================================================================

Generate a gaussian input file for thiophene by running ::

   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hydrogen  ) " -b ../../BuildingBlocks-release  -r 1 


this creates the file::

   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.xyz

which is just a file containing the cartesian coordinates of a
thiophene molecule.  You can view with your favorite viewer. The -r option is set to 1 to generate a single molecule rather than an oligomer 

To create a gaussian input file you can add the gaussian templates to the current directory::
   
   cp ../../donoracceptor.* ./
   
   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hydrogen  ) " -b ../../BuildingBlocks-release  -r 1 

creates the file::

	mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.com    

which is the gaussian input file based on "donoracceptor.com.template"::

	mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.pbs    

which is the pbs script file based on "donoracceptor.pbs.template"::

	mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.json   

which is the json file contain the molecular information ::

	mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.meta   
	mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.r1.com 
	mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.r2.com

These templates can be modified for your HPC system or what ever properties of the molecule you wish to calculate 
