#
# Help message
#
usage() {
    echo " "
    echo "Usage: check.sh [arg]"
    echo "    arg -- new:        Runs tests and copies output to results directory"
    echo "        -- all:        Runs all tests and compares to files in results directory"
    echo "        -- 'test-name' Run compare for 'test-name'"
    echo " "
    echo " "
    echo "To add a serial test:"
    echo "   1. create test with name test-*.py "
    echo "      Note: test-*.py is assumed to print results to screen \n"
    echo "   2. do 'check.sh run' \n"
    echo "   3. double check that 'check.sh compare' gives no errors \n"
    echo "   4. do 'check.sh new' to put new results in /tools/tests/results \n"
    echo "   5. make sure to git add/commit/push new test and new results \n"
    echo " "
    echo "To add a parallel (np = 2) test:"
    echo "   1. create test with name test_n2-*.py "
    echo "      Note: test-*.py is assumed to print results to screen \n"
    echo "   2. Repeat steps 2-5 for serial"
    echo " "
    echo "Note: if parallel np > 2 tests needed script will need to be edited"
    echo " "
}


#
# Runs test comparison and checks with results in ./results and assumes
# form of the file testName.txt
#
compareTest() {
    testName=$1
    runCmd=$2
    echo "------ Running test $testName --------"
    $runCmd $testName > tmp
    diff tmp results/$testName.txt
    rm -rf tmp
    echo "-------------------------------------------------"
    echo " "
}

#
# Checks new results into repo
#
newTest() {

    testName=$1
    runCmd=$2
    echo "------ Checking in results for $testName --------"
    $runCmd $testName > results/$testName.txt
    echo "-------------------------------------------------"
    echo " "
}


if [ $# == 0 ]; then
    usage

elif [ $1 == "new" ]; then

    echo " "
    echo "Checking in new results for python test files"
    echo "Looking for tests named 'test-*.py'  "
    echo " "

    # Serial test
    testNames=`ls -1 test-*.py`
    for testName in $testNames; do
	newTest $testName
    done

    # Parallel test (np = 2)
    testNames=`ls -1 test_n2-*.py`
    for testName in $testNames; do
	newTest $testName 'mpirun -n 2'
    done


elif [ $1 == "all" ]; then

    echo " "
    echo "Comparing results against python test files"
    echo "Looking for tests named 'test-*.py'  "
    echo " " 

    # Serial test
    testNames=`ls -1 test-*.py`
    for testName in $testNames; do
	compareTest $testName
    done

    # Parallel test (np = 2)
    testNames=`ls -1 test_n2-*.py`
    for testName in $testNames; do
	compareTest $testName 'mpirun -n 2'
    done

    echo "If no output (other than status messages)... tests passed"
    echo " "

else
    compareTest $1
fi
