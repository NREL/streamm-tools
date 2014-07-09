#
# Help message
#
usage() {
    echo " "
    echo "Usage: check.sh [arg]"
    echo "    arg -- new:     Runs tests and copies output to results directory"
    echo "        -- compare: Runs tests and compares to files in results directory"
    echo " "
    echo "To add a test:"
    echo "  1. create test with name test-*.py "
    echo "  2. do 'check.sh run'"
    echo "  3. double check that 'check.sh compare' gives no errors"
    echo " "
}


#
# Runs test comparison and checks with results in ./results and assumes
# form of the file testName.txt
#
compareTest() {

    testName=$1
    echo "------ Running test $testName --------"
    $testName > tmp
    diff tmp results/$testName.txt
    rm -rf tmp
    echo "-------------------------------------------------"
    echo " "
}

#
# Checks new results in
#
newTest() {

    testName=$1
    echo "------ Checking in results for $testName --------"
    $testName > results/$testName.txt
    echo "-------------------------------------------------"
    echo " "
}


if [ $# == 0 ]; then
    usage

elif [ $1 == "new" ]; then

    echo " "
    echo "Checking in new results for python test files"
    echo " " 

    testNames=`ls -1 test-*.py`
    for testName in $testNames; do
	newTest $testName
    done


elif [ $1 == "compare" ]; then

    echo " "
    echo "Comparing results against python test files"
    echo " " 

    testNames=`ls -1 test-*.py`
    for testName in $testNames; do
	compareTest $testName
    done

    echo "If no output (other than status messages)... tests passed"
    echo " "


else
    echo "Argument not recognized"
    usage
fi
