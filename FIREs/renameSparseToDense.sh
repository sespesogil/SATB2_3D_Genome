#!/bin/bash
# Change output from HiC-Pro sparseTodense into FIRE required input
# author=sergio.espeso-gil

dir=/analysis/HiC-analysis/FIRE/pooled_10kb_ICED
# store the files:
find $dir -name \*_dense.matrix | sort | uniq > $dir/filestorename.txt
list=$dir/filestorename.txt

# rename and gzip files:
while IFS= read -r file
do
string=$file
suffix="_dense.matrix"
foo=${string%"$suffix"}
echo "${foo}"
mv $file $foo 
gzip $foo
done < "$list"

rm $dir/*/*GL* $dir/*/*JH* $dir/*/*chrM*
