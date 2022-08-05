

DecayFreq<-function(ICED_normalized_path, resolution="100000"){

files<-dir(input, pattern=paste0("_",resolution,_iced.matrixCIS.txt"))

for (i in 1:length(files)){
# iterate through the different files:
  if (i == 1){
  todecay<-files[i]
  filename<-basename(todecay)
  file<-read.table(todecay, sep="\t", header=T)
  file$distance<-(file$T_end - file$A_end)
  decay<-aggregate(value ~ distance , file , mean)
  decay.transf<-data.frame(dis=log10(decay$distance), freq=log10(decay$value))
  decay.transf$group<-filename
  } else {
  todecay<-files[i]
  filename<-basename(todecay)
  file<-read.table(todecay, sep="\t", header=T)
  file$distance<-(file$T_end - file$A_end)
  decay<-aggregate(value ~ distance , file , mean)
  tmp.transf<-data.frame(dis=log10(decay$distance), freq=log10(decay$value))
  tmp.transf$group<-filename
  decay.transf<-rbind(decay.transf,  tmp.transf)
  }
}
write.table(decay.transf, paste0("decay.groups.", resolution, .txt") , sep="\t", col.names=T, quote=F, row.names=F)


}

