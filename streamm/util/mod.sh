

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

for d in "buildingblocks" "forcefields" "mpi" "structures" "util";
do
echo $d;
cp calculations/tests/streamm_testutil.py  ${d}/tests/
git add ${d}/tests/*py
    
done

