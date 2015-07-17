.. _example_da_3:


Example 3: Gaussian input file for oligo hexylthiophene (P3HT)
========================================================================================

Generate a gaussian input file for oligo hexylthiophene (P3HT) ::

   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks-release  -r 5

will generate a rigioregular oligomer of methylthiol n=1-5 ::

   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks-release  -r 5  -p "0"

sets the inter-ring dihedral angle to zero making the sulfurs of thiophene in the cis configuration ::

   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks-release  -r 5  -p "180 0 "

sets the inter-ring dihedral angle to alternate between sero and 180 making the sulfurs of thiophene in the trans configuration 
however, the alkyl chains end up overlapping ::

   python ../da_builder/donoracceptorsystems.py  "thiophene ( R_hexane  ) " -b ../../BuildingBlocks-release  -r 5  -p "140 40 "

increases the inter-ring dihedral angle to remove the overlap between alkyl chains 
