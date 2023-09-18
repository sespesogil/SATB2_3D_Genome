#' TAD N analysis
#' A TAD+N analysis computes the interaction density within TADs and their 1,2,. . . ,N neigh- bours. 
#' This can be used to compare whether TADs in two samples interact differently with their neighbouring TADs.
#' Requires Processed interaction and insulation data at a given resolution. 

library(GENOVA)

# load data:
Dir<-"/hic_results/matrix/"
ProcessedData_10kb<-readRDS(paste0(Dir, "ProcessedData_10000.rds"))

# For insulation computed with GENOVA:
load(paste0(Dir, "/ProcessedData_10kb_insulation.rda"))
TADcalls <- call_TAD_insulation(ProcessedData_10kb_insulation)

TAD_N_flxBIC   <- intra_inter_TAD(list("flxBIC" = ProcessedData_10kb$flxBIC,
                                   'ckoBIC' = ProcessedData_10kb$ckoBIC),
                              tad_bed = TADcalls$flxBIC,
                              max_neighbour = 10)

saveRDS(TAD_N_flxBIC, paste0(Dir, "TAD_N_flxBIC_intra_interTAD.10kb.wTADscomputedGENOVA.rds"))

TAD_N_flxBICjitter<-visualise(TAD_N_flxBIC, geom = 'jitter')
TAD_N_flxBICviolin<-visualise(TAD_N_flxBIC, geom = 'violin')

save(TAD_N_flxBICjitter, TAD_N_flxBICviolin, file=paste0(Dir, "plots.TAD_N_flxBIC_intra_interTAD.10kb.wTADscomputedGENOVA.rda"))
