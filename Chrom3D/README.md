#  Chrom3D analysis

3D modelling was computed using: 

- 1. [https://github.com/sespesogil/automat_chrom3D](https://github.com/sespesogil/automat_chrom3D) : generation of Gtracks need to run [Chrom3D](https://github.com/Chrom3D/Chrom3D)
- 2. Chrom3D run: script shows example run for one of the samples. Briefly, Chrom3D was run with the nucleus parameter set, with a radius of 3 µm for 1 million iterations: -r 3.0 -n 1000000 -l 5000 --nucleus.
- 3. [https://github.com/sespesogil/automat_chrom3D_colors](https://github.com/sespesogil/automat_chrom3D_colors): used to map regions of interest.
 
# Installation 

```conda environment

conda create --prefix ./chrom3D python=2.7.17

conda activate /tools/chrom3D 
conda install -c bioconda pybedtools
conda install -c bioconda R
conda install -c bioconda statsmodels 
conda install -c bioconda bedtools 
conda install -c bioconda ggplot2 

module load boost/1.72.0-gcc-10.2.0-jlf5edv```
