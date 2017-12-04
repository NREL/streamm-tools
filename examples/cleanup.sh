
for f in *.ipynb 
do
    echo $f 
    jupyter nbconvert --to rst $f
done

for f in *.rst
do
    echo $f 
    mv $f temp.rst
    sed 's/ipython2/python/g'  temp.rst > $f
done

rm getting_started.rst 
rm temp.rst
rm resource_peregrine.rst

for f in *.rst 
do
    tag=${f%.rst}
    echo $f $tag
    echo ".. _$tag:" > header.txt
    echo "  "  >> header.txt

    echo "$tag"  >> header.txt
    echo "========================"  >> header.txt
    echo " "  >> header.txt
    cat header.txt $f  > ../doc/source/_static/$f
    
done

rm -rf *.data *.dat *.xyz  *.json *.log *.pkl *.rst *nw *.in  *csv *tgz materials/ scratch/   scripts/   storage/ templates/
