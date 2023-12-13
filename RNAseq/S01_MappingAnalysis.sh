#!/bin/bash
##
#SBATCH -N 1
#SBATCH --partition=skylake_0384
#SBATCH --qos=skylake_0384
#SBATCH --ntasks-per-node=48
#SBATCH -J RNAseq

source activate /SergioEG/tools/DeepTools

STAR=/tools/STAR-2.7.6a/bin/Linux_x86_64/STAR
bamCoverage=/tools/deepTools/bin/bamCoverage
featureCOUNTS=/tools/subread-2.0.1-source/bin/featureCounts
samtools=/tools/samtools-1.11/samtools
trimmonatic=/tools/trimmomatic/classes/trimmomatic.jar

in=/gpfs/data/fs71524/SergioEG/fastq_repository/RNAseq.rescue
in2=/gpfs/data/fs71524/SergioEG/fastq_repository
out=/gpfs/data/fs71524/SergioEG/analysis/RNA-seq.Rescue/mouse.rescue

genomeDir=/resources/genomes/mouse/GRCm38.p6.vM25/STAR

REF=/resources/genomes/mouse/GRCm38.p6.vM25/GRCm38.p6.genome.fa
GTF=resources/genomes/mouse/GRCm38.p6.vM25/gencode.vM25.annotation.gtf


find $in -name \*.fq.gz -printf '%f\n' |  awk -F '_E' '{print $1}' | sort | uniq > $in2/filenames.txt


cd $out/flx_1

$STAR --runThreadN 48 --genomeDir $genomeDir --sjdbGTFfile $GTF --sjdbOverhang 100 --readFilesIn /fastq_repository/flx_1/flx_1.R1.fastq /fastq_repository/RNAseq.rescue/flx_1.R2.fastq --twopassMode Basic --outSAMtype BAM SortedByCoordinate Unsorted --quantMode TranscriptomeSAM GeneCounts
$samtools index $out/flx_1/Aligned.sortedByCoord.out.bam
$samtools flagstat $out/flx_1/Aligned.sortedByCoord.out.bam > $out/flx_1/flx_1.report.txt
$samtools index $out/flx_1/Aligned.sortedByCoord.bam
bamCoverage -b $out/flx_1/Aligned.sortedByCoord.bam -o $out/flx_1/flx_1.bw -p max
mv $out/flx_1/Aligned.sortedByCoord.out.bam $out/flx_1/Aligned.sortedByCoord.bam


$featureCOUNTS -T 40 --ignoreDup -a $GTF -g gene_name -O -o $out/flx_1/flx_1.exon.gene_name.txt $out/flx_1/Aligned.sortedByCoord.bam
awk -F"\t" '{print $7}' $out/flx_1/flx_1.exon.gene_name.txt > $out/flx_1/flx_1.exon.gene_name.counts.txt
awk -F"\t" '{print $1}' $out/flx_1/flx_1.exon.gene_name.txt > $out/features.txt
