echo " "
echo "Making pydoc HTML files"

scpHTMLtoStatic() {

    modName=$1 # Module name prefix
    dirName=$2 # Directory on cs to copy to

    pydoc -w $modName
    sleep 0.5

    echo "Moving developer docs ($modName) to --> cs.hpc.nrel.gov:static/prospect [y] yes [n]no"
    read -e ans
    if [ $ans == "y" ]; then
	echo " "
	scp -r $modName.html cs.hpc.nrel.gov:static/$modName
    else
	echo " ";
    fi
}


#
# Clean up
#
rm -rf *.html
echo " "

scpHTMLtoStatic "atom_types"  "STREAMM"
scpHTMLtoStatic "cluster"     "STREAMM"
scpHTMLtoStatic "file_io"  "STREAMM"
scpHTMLtoStatic "gromacs"  "STREAMM"
scpHTMLtoStatic "jsonapy"  "STREAMM"
scpHTMLtoStatic "nwchem"  "STREAMM"
scpHTMLtoStatic "prop"  "STREAMM"
scpHTMLtoStatic "xmol"  "STREAMM"
scpHTMLtoStatic "calc_qm"  "STREAMM"
scpHTMLtoStatic "elements"  "STREAMM"
scpHTMLtoStatic "gaussian"  "STREAMM"
scpHTMLtoStatic "groups"  "STREAMM"
scpHTMLtoStatic "lammps"  "STREAMM"
scpHTMLtoStatic "pdb"  "STREAMM"
scpHTMLtoStatic "top"  "STREAMM"
