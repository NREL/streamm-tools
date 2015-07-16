.. _examples:

*************************
Examples
*************************

Selected tools-tests
========================================================
These examples are closely related to regression tests in the
tools-tests repo. 

.. toctree::
   :maxdepth: 1

   example_subStructure.rst


IPython notebook downloads
========================================================
Each of the selections above is available as an IPython notebook, so you may explore these examples further.
After clicking on the download link, drag the desired file to a location on your computer.

:download:`example_subStructure.ipynb <example_subStructure.ipynb>`

The IPython web interface for a notebook can be started from the command line by e.g. ::

   ipython notebook example_subStructure.ipynb



.. Removing this code since OPV-database description is out for release
.. More complicated example with simulation objects
.. ========================================================
.. Functional/basis 'spamming' description (SCOTT)


Molecular generation
========================================================
 Examples of running the donoracceptorsystems.py script to build complex molecules 

Files 
donoracceptorsystems.py 
  Concatinates molecular cply files based on conectivety tags 
BuildingBlocks
  Must contain the directories 
     acceptors         donors            functional_groups spacers           terminals
  The backbone of the molecule can be composed of  acceptors, spacers or donors, and can be 
  decorated with molecules from the   functional_groups directory 
  and terminated with molecules from the terminals directory 
  Run "python ../da_builder/donoracceptorsystems.py -h" to see a list of the commands 


 Example 1
   Generate a gaussian input file for thiohpene by ::
   
   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hydrogen  ) " -b ../../BuildingBlocks_release  -r 1 


 creates the file:
   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.xyz
which is just a file containing the cartision cordinates of a thiophene molecule
which you can view with your favorate viewer 

 -r 1 is set to 1 to generate a single molecule rather than an oligomer 

To create a gaussian input file you can add the gaussian templates to the current directory

cp ../../donoracceptor.* ./
python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hydrogen  ) " -b ../../BuildingBlocks_release  -r 1 

 creates the file:
   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.com    
which is the gaussian input file based on "donoracceptor.com.template"
   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.pbs    
which is the pbs script file based on "donoracceptor.pbs.template"
   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.json   
which is the json file contain the molecular information 
   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.meta   
   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.r1.com 
   mols/thiophene_R_hydrogen_/acc1_thiophene_R_hydrogen__n1.r2.com

These templates can be modified for your HPC system or what ever properties of the molecule you wish to calculate 

Example 2
   Generate a gaussian input file for oligo methylthiol 

python ../da_builder/donoracceptorsystems.py  "thiophene ( R_methane  ) " -b ../../BuildingBlocks_release  -r 5

will generate a rigioregular oligomer of methylthiol n=1-5 

python ../da_builder/donoracceptorsystems.py  "thiophene ( R_methane  ) " -b ../../BuildingBlocks_release  -r 5  -p "0"

sets the inter-ring dihedral angle to zero making the sulfurs of thiophene in the cis configuration 

python ../da_builder/donoracceptorsystems.py  "thiophene ( R_methane ) " -b ../../BuildingBlocks_release  -r 5  -p "180 0 "

sets the inter-ring dihedral angle to alternate between sero and 180 making the sulfurs of thiophene in the trans configuration 


Example 3
   Generate a gaussian input file for oligo hexylthiol (P3HT)

python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks_release  -r 5

will generate a rigioregular oligomer of methylthiol n=1-5 

python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks_release  -r 5  -p "0"

sets the inter-ring dihedral angle to zero making the sulfurs of thiophene in the cis configuration 

python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks_release  -r 5  -p "180 0 "

sets the inter-ring dihedral angle to alternate between sero and 180 making the sulfurs of thiophene in the trans configuration 
however, the alkyl chains end up overlapping 

python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks_release  -r 5  -p "140 40 "

increases the inter-ring dihedral angle to remove the overlap between alkyl chains


Coming soon

- running MD sim (LAMMPS, Gromacs)
- calculation of RDFsd

- setup of system inputs
- running MD sim (LAMMPS, Gromacs)
- calculation of RDFs
