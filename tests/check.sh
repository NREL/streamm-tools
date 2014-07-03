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

elif [ $1 == "compare" ]; then

    echo " "
    echo "Comparing results against python test files"
    echo " " 

    test-bondContainer.py > tmp            ; diff tmp results/test-bondContainer.txt        ; echo " "
    test-particleConstructors.py > tmp     ; diff tmp results/test-particleConstructors.txt ; echo " "
    test-particleSetInfo.py > tmp          ; diff tmp results/test-particleSetInfo.txt      ; echo " "
    test-searchTags.py      > tmp          ; diff tmp results/test-searchTags.txt           ; echo " "
    test-checkTypes.py      > tmp          ; diff tmp results/test-checkTypes.txt           ; echo " "
    test-particleContainer.py > tmp        ; diff tmp results/test-particleContainer.txt    ; echo " "
    test-ptclContainerConstructor.py > tmp ; diff tmp results/test-ptclContainerConstructor.txt ; echo " "
    test-subStructure.py  > tmp            ; diff tmp results/test-subStructure.txt             ; echo " "
    test-strucDumpSave.py > tmp            ; diff tmp results/test-strucDumpSave.txt            ; echo " "

    rm -rf tmp
    echo "If no output... tests passed"
    echo " "

else
    echo "Argument not recognized"
    usage
fi
