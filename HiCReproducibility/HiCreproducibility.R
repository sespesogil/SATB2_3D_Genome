#' HiCreproducibility function
#' Computes SCC reproducibitlity scores
#' @param path directory were to find hic matrices 
#' @param resol desired resolution 
#' @param h smoothing parameter. 
#' @param lbr  lower bounds of the genomic distance between interaction loci.
#' @param ubr  upper bounds of the genomic distance between interaction loci.
#' @author sergio.espeso-gil
#' @export 

HiCreproducibity<-function(path, resol, h , lbr , ubr){
files<-list.files(path, pattern=".matrix")
for (i in 1:length(files)){
  file1<-files[i]
  filename1<-basename(gsub(".matrix", "", files[i]))
  if (i == 1){
    mat1<-read.table(paste(path, file1, sep="/"))
    for (j in 1:length(files)){
      if (j == 1) {
      file2<-files[j]
      filename2<-basename(gsub(".matrix", "", files[j]))
      mat2<-read.table(paste(path, file2, sep="/"))
      scc.out = get.scc(mat1, mat2, resol, h, lbr, ubr)
      SCC_rep<-as.data.frame(scc.out$scc)
      colnames(SCC_rep)<-c("reproducibility")
      SCC_rep$group<-paste(filename1, filename2, sep="-")
     } else {
      file2<-files[j]
      filename2<-basename(gsub(".matrix", "", files[j]))
      mat2<-read.table(paste(path, file2, sep="/"))
      tmp = get.scc(mat1, mat2, resol, h, lbr, ubr)
      tmp2<-as.data.frame(tmp$scc)
      colnames(tmp2)<-c("reproducibility")
      tmp2$group<-paste(filename1, filename2, sep="-")
      SCC_rep <- rbind(SCC_rep, tmp2)
    }
  }
  } else {
    mat1<-read.table(paste(path, file1, sep="/"))
    for (j in 1:length(files)){
      file2<-files[j]
      filename2<-basename(gsub(".matrix", "", files[j]))
      mat2<-read.table(paste(path, file2, sep="/"))
      tmp3 = get.scc(mat1, mat2, resol, h, lbr, ubr)
      tmp4<-as.data.frame(tmp3$scc)
      colnames(tmp4)<-c("reproducibility")
      tmp4$group<-paste(filename1, filename2, sep="-")
      SCC_rep <- rbind(SCC_rep, tmp4)
  }
}
}
return(SCC_rep)
}

