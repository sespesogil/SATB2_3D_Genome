# First steps / HiC matrix 

library("FIREcaller")
FIREpath<-"/analysis/HiC-analysis/FIRE/pooled_10kb_ICED"

fileList<-list.files(path=FIREpath, pattern=".gz", recursive=TRUE)
fileNames<-paste(unique(dirname(fileList)))

for (i in 1:length(fileNames)){
setwd(paste0(FIREpath, "/", fileNames[i]))
cat(paste0('Calculating FIRE calls for ', fileNames[i]))
sampleName<-paste0(fileNames[i], "_10000_chr") 
file.list <- paste0(sampleName, 1:19,'.gz') 
file.list<-append(file.list, paste0(sampleName, "X",'.gz'))
file.list<-append(file.list, paste0(sampleName, "Y",'.gz'))
map_file<-'/analysis/HiC-analysis/FIRE/resources/mm10_F_GC_M_arima_10Kb_el_auto.txt.txt'
FIREcaller(file.list, gb="mm10", map_file, binsize = 10000, upper_cis = 200000, normalized = FALSE, 
           rm_mhc = TRUE, rm_EBL = TRUE, rm_perc = 0.25, dist = 'poisson', alpha = 0.05, plots = FALSE, diff_fires = FALSE)
}
