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

    tst1.py > results/tst1.dat
    tst2.py > results/tst2.dat
    tst3.py > results/tst3.dat
    tst4.py > results/tst4.dat
    tst6.py > results/tst6.dat
    tst7.py > results/tst7.dat
    tst8.py > results/tst8.dat


elif [ $1 == "compare" ]; then

    echo " "
    echo "Comparing results against python test files"
    echo " " 

    tst1.py > tmp ; diff tmp results/tst1.dat
    tst2.py > tmp ; diff tmp results/tst2.dat
    tst3.py > tmp ; diff tmp results/tst3.dat
    tst4.py > tmp ; diff tmp results/tst4.dat
    tst6.py > tmp ; diff tmp results/tst6.dat
    tst7.py > tmp ; diff tmp results/tst7.dat
    tst8.py > tmp ; diff tmp results/tst8.dat

    echo "If no output... tests passed"
    echo " "

else
    echo "Argument not recognized"
    usage
fi
