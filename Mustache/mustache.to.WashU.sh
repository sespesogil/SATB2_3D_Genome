#!/bin/bash

# Script to create WashU files from mustache ouput
# author=sergio.espeso-gil

dir=/analysis/hic_results/data/
toname=$(echo $dir)

# find Mustache chromatin loops and create a list:
name=$(basename $toname) 
find $dir -name \*pv001.tsv > $dir/mustache.files.txt 
list=$dir/mustache.files.txt

# clean the headers:
while IFS= read -r file
do
 sed 1d $file > $file.formatted.txt
 awk -F" " '{print "chr"$1"\t"$2"\t"$3"\t""chr"$4":"$5"-"$6","$8"\t""."}' $file.formatted.txt > $file.WashU.score.txt
 awk -F" " '{a = -log($7)/log(10); print "chr"$1"\t"$2"\t"$3"\t""chr"$4":"$5"-"$6","a"\t""."}' $file.formatted.txt > $file.WashU.FDR.txt
done < "$list"

# concatenate the files
cat $dir/*WashU.score.txt > $dir/$name.WashU.score.txt
cat $dir/*WashU.FDR.txt > $dir/$name.WashU.FDR.txt

# differential loops :
dir=/analysis/HiC-analysis/Mustache/diff.loops
toname=$(echo $dir)

# find Mustache chromatin loops and create a list:
name=$(basename $toname) 
find $dir -name \output.082020201.5kb.* > $dir/mustache.files.txt 
list=$dir/mustache.files.txt

# clean the headers:
while IFS= read -r file
do
 sed 1d $file > $file.formatted.txt
 awk -F" " '{print "chr"$1"\t"$2"\t"$3"\t""chr"$4":"$5"-"$6","$8"\t""."}' $file.formatted.txt > $file.WashU.score.txt
 awk -F" " '{a = -log($7)/log(10); print "chr"$1"\t"$2"\t"$3"\t""chr"$4":"$5"-"$6","a"\t""."}' $file.formatted.txt > $file.WashU.FDR.txt
done < "$list"


# concatenate the files
cat $dir/*WashU.score.txt > $dir/$name.WashU.score.txt
cat $dir/*WashU.FDR.txt > $dir/$name.WashU.FDR.txt
