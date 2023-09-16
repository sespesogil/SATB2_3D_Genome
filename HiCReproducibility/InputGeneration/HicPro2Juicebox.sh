#!/bin/bash
##
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH -J juicebox


### DO NOT FORGET CROSS-CHECK FRAG

HICPRO_PATH=/tools/HiC-Pro2.11.4.installed/HiC-Pro_2.11.4/bin/utils
hic_pro_path=/hic_results/data/cKO_BIC/cKO_BIC.allValidPairs
JUICEBOX=/tools/juicer/SLURM/scripts/Juicebox_1.11.08.jar
FRAG=/resources/genomes/mouse/GRCm38.p6.vM25/hicpro/mm10_arima.bed

cd /hic_results/data/cKO_BIC/

$HICPRO_PATH/hicpro2juicebox.sh -i $hic_pro_path -g mm10 -j $JUICEBOX -r $FRAG

