#!/bin/bash
# Script to convert Mustache loops into FANC format 

dir=/gpfs/data/fs71524/SergioEG/analysis/HiC-analysis/Mustache.c.loops/pooled.loops/5kb/flxnbqx
toname=$(echo $dir)
res=loops.5kb

# find Mustache chromatin loops and create a list:
name=$(basename $toname) 
find $dir -name \*pv001.tsv > $dir/mustache.files.txt 
list=$dir/mustache.files.txt

# clean the headers 

while IFS= read -r file
do
 sed 1d $file > $file.formatted.txt
 awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"".""\t"$8}' $file.formatted.txt > $file.toFANC.bedpe
done < "$list"


# concatenate the files
cat $dir/*toFANC.bedpe > $dir/$name.$res.toFANC.bedpe

# Second pooled samples batch:
dir=/gpfs/data/fs71524/SergioEG/analysis/HiC-analysis/Mustache.diff.test
toname=$(echo $dir)
res=loops.5kb

# find Mustache chromatin loops and create a list:
name=$(basename $toname) 
find $dir -name \output.082020201.5kb.* > $dir/mustache.files.txt 
list=$dir/mustache.files.txt

# clean the headers 

while IFS= read -r file
do
 sed 1d $file > $file.formatted.txt
 awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"".""\t"$8}' $file.formatted.txt > $file.toFANC.bedpe
done < "$list"


# concatenate the files
cat $dir/*toFANC.bedpe > $dir/$name.$res.toFANC.bedpe
