# Differential compartment analysis using [dcHiC](https://github.com/ay-lab/dcHiC)

Analysis was performed using [dcHiC-v1](https://github.com/ay-lab/dcHiC/wiki/dcHiC-V1-Wiki-Documentation). 

Steps:

1. Preprocessing: HiC-Pro files were parsed to [Juicer]() hic file format due some issues while using dcHiC-v1 directly with HiC-Pro files. Then O/E was computed using 'preprocessing.py' script following dcHiC guidelines.  
2. dHiC analysis

Some extra script have been added used to associate differential compartments to mouse genomic annotation. 
