#!/bin/bash
##
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH -J mustache


/tools/Python-3.8.7/python.new.sql/bin/python3.8 /tools/mustache/mustache/diff_mustache.py -f1 /analysis/HiC-analysis/hic_results/data/cKO_BIC/cKO_BIC.allValidPairs.hic -f2 /analysis/HiC-analysis/hic_results/data/flx_BIC/flx_BIC.allValidPairs.hic -pt 0.05 -pt2 0.05 -o /analysis/HiC-analysis/diff.ckoBIC_vs_flxBIC/output.082020201.5kb -r 5000 -st 0.8
