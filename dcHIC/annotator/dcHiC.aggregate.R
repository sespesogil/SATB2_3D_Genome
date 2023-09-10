options(echo=T)
args <- commandArgs(trailingOnly = TRUE)
annotation <- args[1]
name<-basename(gsub(".ann.txt", "", annotation))

header<-read.table("./header.txt", header=T)
ann<-read.table(annotation, header=F , stringsAsFactors=F)
colnames(ann)<-names(header)
ann.collapsed<-aggregate(genes ~., ann, toString)
write.table(ann.collapsed, paste(name,".annotation.txt", sep=""), col.names=TRUE, quote= FALSE, row.names=FALSE, sep="\t") 
