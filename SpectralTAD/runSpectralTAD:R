require(SpectralTAD)
library(dplyr)
library(data.table)

HiCproPath<-"/analysis/HiC-analysis/hic_results/matrix/"
files<-c("cKO_BIC","flx_BIC")
out.dir<-"/HiC-analysis/SpectralTAD/"
dir.create(out.dir)
res<-"25000"

lapply(seq_along(files), function(i){
mat<-read.table(paste0(HiCproPath, files[i], "/raw/", res, "/", files[i], "_", res, ".matrix"))
bed <- read.table(paste0(HiCproPath, files[i], "/raw/", res, "/", files[i], "_" , res, "_abs.bed" ))  
sparse_mats <- HiCcompare::hicpro2bedpe(mat,bed)
#remove chrM
sparse_mats$cis = sparse_mats$cis[-20]
TADs<-lapply(seq_along(sparse_mats$cis), function(j){
x <- sparse_mats$cis[[j]]
#Pull out chromosome
CHR <- x[, 1][1]
#Subset to make three column matrix
x <- x[, c(2, 5, 7)]
# Call TADs using SpectralTAD
bed_coords <- bind_rows(SpectralTAD(x, chr = CHR , levels = 1)) 
return(bed_coords)
 })
TADs.df<-rbindlist(TADs)
write.table(TADs.df, paste0(out.dir, files[i], "TADs.1_levels.txt"), quote=F, sep="\t", row.names=F)
})
