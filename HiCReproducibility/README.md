# HiC reproducibility 

In order to assess the reproducibility of our HiC samples, [HiCRep method](https://github.com/TaoYang-dev/hicrep) was wrapped into HiCReproducibility function. This wrapper and function modification was needed due some issues with the HiCrep package related to the straw package. 

This repo contains the main function, running script and SLURM launcher example. 

## Requirements: 

### Input: 

- Juicebox hic files: they can be produced using [hicpro2juicebox.sh](https://github.com/nservant/HiC-Pro/blob/master/bin/utils/hicpro2juicebox.sh) script. 
