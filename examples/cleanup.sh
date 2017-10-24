
for f in *.ipynb 
do
    echo $f 
    jupyter nbconvert --to rst $f
done


rm -rf *.dat *.xyz  *.json *.log *pkl materials/ scratch/   scripts/   storage/
