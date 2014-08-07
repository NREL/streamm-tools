echo " "
echo "Making pydoc HTML files"


scpHTMLtoStatic() {

    modName=$1 # Module name prefix
    dirName=$2 # Directory on cs to copy to

    mkdir $dirName
    pydoc -w ../src/$modName.py
    sleep 0.5

    mv $modName.html $dirName

    # echo "Moving developer docs ($modName) to -->
    # cs.hpc.nrel.gov:static/$dirName [y] yes [n]no"
    # read -e ans
    #if [ $ans == "y" ]; then
    #	echo " "
    #	scp -r $modName.html cs.hpc.nrel.gov:static/$modName
    #else
    #	echo " ";
    #fi
}

#
# Clean up
#
rm -rf *.html
echo " "

scpHTMLtoStatic "mpiNREL"            "mpiNREL"
scpHTMLtoStatic "runjobs"            "runjobs"

scpHTMLtoStatic "particles"                  "STREAMM"
scpHTMLtoStatic "bonds"                      "STREAMM"
scpHTMLtoStatic "angles"                     "STREAMM"
scpHTMLtoStatic "structureContainer"         "STREAMM"
scpHTMLtoStatic "lammps1_StructureContainer" "STREAMM"
scpHTMLtoStatic "simulation"                 "STREAMM"

echo " "
echo "Run: sftp cs.hpc.nrel.gov:static"
echo "      put -r * ."

#sftp cs.hpc.nrel.gov:
#sftp>
#   *  put  -r static
#   * quit

