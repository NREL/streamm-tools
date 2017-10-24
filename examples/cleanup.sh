
for f in *.ipynb 
do
    echo $f 
    jupyter nbconvert --to rst $f
done



for f in *.rst 
do
    tag=${f%.rst}
    echo $f $tag
    echo ".. _$tag:" > header.txt
    echo "  "  >> header.txt

    echo "$tag"  >> header.txt
    echo "==============="  >> header.txt
    echo " "  >> header.txt
    cat header.txt $f  > ../doc/source/_static/$f
    
done

rm -rf *.dat *.xyz  *.json *.log *pkl *rst materials/ scratch/   scripts/   storage/
