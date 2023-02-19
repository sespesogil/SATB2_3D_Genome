# Analysis overview

Pair-end data were mapped by using BWA-MEM (v0.7.17) 103 and the resulting SAM files were
converted into sorted and indexed BAM files (SAMtools v1.1) 104 . Reproducibility of replicate
ATAC-seq libraries was assessed by calculating correlation between BAM files by using
deepTools2 105 . Peaks were identified per experimental condition by employing Genrich tool
(v0.6.1). Additionally, a consensus peak file across all conditions was produced by using
bedtools “merge”. The number of read counts per region was calculated by using
featureCounts (subread v.2.0.1) 106 . Next, differential analysis was done by using edgeR 107 .
Principal component analysis showed a technical bias and samples were batch-normalized
using RUVr method of RUVseq 108 . Homer (v4.11) was used for motif analysis. Genomic
annotation and GO enrichment analysis was done by using ChiPseeker 29 and clusterProfiler 109 .
Footprinting analysis and occupancy prediction was performed by using Tobias (v0.12.10) 45
over 5 kb-targets of promoter-based loops. Proximal Flx and cKO dOCRs were annotated
using GREATs 59 basal plus extention mode (Proximal: 5.0 kb upstream, 1.0 kb downstream,
plus Distal: up to 2.5 kb).
