.. _molgenex2:


Example 2: Gaussian input file for oligo methylthiophene
============================================================================================

Generate a gaussian input file for oligo methylthiophene by running ::
   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_methane  ) " -b ../../BuildingBlocks-release  -r 5

will generate a rigioregular oligomer of methylthiol n=1-5 ::

   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_methane  ) " -b ../../BuildingBlocks-release  -r 5  -p "0"

sets the inter-ring dihedral angle to zero making the sulfurs of thiophene in the cis configuration ::

   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_methane ) " -b ../../BuildingBlocks-release  -r 5  -p "180 0 "

sets the inter-ring dihedral angle to alternate between zero and 180 making the sulfurs of thiophene in the trans configuration 
