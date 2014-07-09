#
# Help message
#
usage() {
    echo " "
    echo "Usage: check.sh [arg]"
    echo "    arg -- new:     Runs tests and copies output to results directory"
    echo "        -- compare: Runs tests and compares to files in results directory"
    echo " "
}


#
# Runs test comparison and checks with results in ./results and assumes
# form of the file testName.txt
#
compareTest() {

    testName=$1
    echo "------ Running test $testName --------"
    $testName.py > tmp
    diff tmp results/$testName.txt
    rm -rf tmp
    echo "-------------------------------------------------"
    echo " "
}

if [ $# == 0 ]; then
    usage

elif [ $1 == "new" ]; then

    echo " "
    echo "Checking in new results for python test files"
    echo " " 

    test-bondContainer.py            > results/test-bondContainer.txt
    test-particleConstructors.py     > results/test-particleConstructors.txt
    test-particleSetInfo.py          > results/test-particleSetInfo.txt
    test-searchTags.py               > results/test-searchTags.txt
    test-checkTypes.py               > results/test-checkTypes.txt
    test-particleContainer.py        > results/test-particleContainer.txt
    test-ptclContainerConstructor.py > results/test-ptclContainerConstructor.txt
    test-subStructure.py             > results/test-subStructure.txt
    test-strucDumpSave.py            > results/test-strucDumpSave.txt
    test-strucDumpSave.py            > results/test-strucDumpSave.txt
    test-strucAdd.py                 > results/test-strucAdd.txt

elif [ $1 == "compare" ]; then

    echo " "
    echo "Comparing results against python test files"
    echo " " 

    compareTest test-bondContainer
    compareTest test-particleConstructors
    compareTest test-particleSetInfo
    compareTest test-searchTags
    compareTest test-checkTypes
    compareTest test-particleContainer
    compareTest test-ptclContainerConstructor
    compareTest test-subStructure
    compareTest test-strucDumpSave
    compareTest test-strucAdd

    echo "If no output (other than status messages)... tests passed"
    echo " "

else
    echo "Argument not recognized"
    usage
fi
