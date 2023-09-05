#' Preparation of HiC matrix files and SCC score computing
#' Calculating SCC scores is computationally expensive.  Run the script by batches of chromosomes
#' @author sergio.espeso-gil

library(strawr)
library(hicrep)

source("/gpfs/data/fs71524/SergioEG/scripts/HiCPro/HiCRep/HiCreproducibility.R")

chromosomes<-c("1","2","3","4","5","6","7")   # select desidered chromosomes 

# function to create matrices from hic files
bed2mat <- function(bed, resol = "NONE"){
  # n the max number of bin
  if (resol == "NONE") {
    N = max(bed[,c(1:2)])
  } else {
    N = max(bed[,c(1:2)])/resol
    bed[,1] = bed[,1]/resol
    bed[,2] = bed[,2]/resol
  }

  mat = matrix(0, N, N)
  mat[bed[,c(1,2)]] = bed[,3]
  mat[bed[,c(2,1)]] = bed[,3]

  return(mat)
}

# indicate the resolution and path
dir="/analysis/HiC-analysis/HiCRep/hicfiles"
# list all hic files
hic_files<-list.files(dir, pattern=".hic", full.names=TRUE)

if(TRUE){
for (i in 1:length(hic_files)){
file<-hic_files[i]
filename1<-gsub(".allValidPairs.hic", "", basename(hic_files[i]))
chromosomes<-c("1","2","3","4","5","6","7")

  for (j in 1:length(chromosomes)){
    chrom<-chromosomes[j]
    filename2<-basename(chromosomes[j])
    tmp<-straw(matrix="observed", norm="NONE", file, chrom, chrom, unit="BP", 250000)
    bed<- as.matrix(tmp)
    mat <- bed2mat(bed, 25000)
    dir.create(paste(dir, "/chr", filename2, sep=""))
    out<-paste(dir, "/chr", filename2, sep="")
    write.table(mat, paste(out, "/", filename1, ".chr", filename2, ".matrix", sep=""), quote=FALSE)
    }

 }
}


# SCC computing: 
for (j in 1:length(chromosomes)){
        chrom<-chromosomes[j]
        filename2<-basename(chromosomes[j])
        out<-paste(dir, "/chr", filename2, sep="")
        scc<-HiCreproducibity(out, resol=250000, h=1 , lbr = 0, ubr = 5000000)
        write.table(scc, paste(out, "/", filename2, "reproducibility.txt", sep=""))
}
