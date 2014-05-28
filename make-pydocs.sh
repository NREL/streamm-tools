echo " "
echo "Making pydoc HTML files"

#
# Clean up
#
rm -rf *.html

echo " "
pydoc -w basics
sleep 1

echo " "
pydoc -w fragments
sleep 1

echo " "
pydoc -w donoracceptorsystems
sleep 1

echo "Moving basics.html to --> cs.hpc.nrel.gov:static/prospect [y] yes [n]no"
read -e ans
if [ $ans == "y" ]; then
    echo " "
    scp -r basics.html cs.hpc.nrel.gov:static/opv_scr_basics
else
    echo " ";
fi

echo "Moving fragments.html docs to --> cs.hpc.nrel.gov:static/prospect [y] yes [n]no"
read -e ans
if [ $ans == "y" ]; then
    echo " "
    scp -r fragments.html cs.hpc.nrel.gov:static/opv_scr_fragments
else
    echo " ";
fi


echo "Moving donoracceptorsystems.html docs to --> cs.hpc.nrel.gov:static/prospect [y] yes [n]no"
read -e ans
if [ $ans == "y" ]; then
    echo " "
    scp -r donoracceptorsystems.html cs.hpc.nrel.gov:static/opv_scr_donoracceptorsystems
else
    echo " ";
fi


