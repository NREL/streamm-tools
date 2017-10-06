

for t in test_*py
do
    echo $t
    
    
    mv $t  temp.py
    sed 's/HOME_DIR = os.getcwd()/from streamm_testutil import */'  temp.py > $t
    mv $t  temp.py
    sed 's/.*DIR.*//'  temp.py > $t    

    if [ ! `grep "setUp_streamm"  $t` ];
    then
        mv $t  temp.py
sed 's/    def setUp/    @setUp_streamm \
    def setUp/g' temp.py > $t 
    fi

    if [ ! `grep "tearDown_streamm"  $t` ];
    then
        mv $t  temp.py
sed 's/    def tearDown/    @tearDown_streamm \
    def tearDown/g' temp.py > $t 
    fi
done

rm temp.py
for t in test_*py
do
    echo $t

    mv $t  temp.py
    sed 's/file_i//'  temp.py > $t    
done


for t in test_*py
do
    echo $t

    mv $t  temp.py
    sed 's/.*os.remove.*//'  temp.py > $t    
done




for d in "buildingblocks" "forcefields" "mpi" "calculations" "util";
do
echo $d;
cp structures/tests/streamm_testutil.py  ${d}/tests/
git add ${d}/tests/*py
    
done





rm temp.py
for t in test_*py
do
    echo $t
    mv $t  temp.py
    sed 's/._matrix//'  temp.py > $t    
    sed 's/set_matrix(matrix_i)/matrix = matrix_i/g'  temp.py > $t    
done

