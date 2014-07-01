mkdir opvtest
cd opvtest

PYTHONPATH=/Users/tkemper/Scripts/streamm/tools/src

python ../../../opv-project/src/donoracceptorsystems.py -f testmol -b ../../../BuildingBlocks-private/   -r 1  " D100 "
python ../../../opv-project/src/donoracceptorsystems.py -f testmol -b ../../../BuildingBlocks-private/   -r 1  "D1  ( R2 R2 0 0 )   "    --make_ff  
python ../../../opv-project/src/donoracceptorsystems.py -f testmol -b ../../../BuildingBlocks-private/   -r 1  "D1  ( R2 R2 0 0 )   A2 ( R3 )"    --make_ff  




# BDT TPD 
# D1_R2R200_A2_R3_n1_R41n1R41n1R40n1


python ../../../opv-project/src/donoracceptorsystems.py -f testmol -b ../../../BuildingBlocks-private/   -r 1  "D1  ( R2 R2 0 0 )   A2 ( R3 )"    --make_ff  
python ~/Scripts/streamm/tools/scripts/rdf.py -j testmol/D1_R2R200_A2_R3_/acc1_D1_R2R200_A2_R3__n1.json 
 

# Need atomicpy for now 
PYTHONPATH=/Users/tkemper/Scripts/streamm/tools/AtomicPy/src/:$PYTHONPATH
python ../../../opv-project/src/donoracceptorsystems.py  -b ../../../BuildingBlocks-private/   -r 1  " D100 "   --make_ff  
