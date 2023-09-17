#' Function to create a bedPe file from HiC-Pro data
#' @param HiCPath directory where to find HiC-Pro matrices
#' @author sergio.espeso-gil
#' @export

importHiC_ICED_interactions<- function (HiCPath){

matrices<-list.files(HiCPath, pattern="\\.matrix$", full.names=TRUE, recursive = TRUE)

RAW<-matrices[grep("/raw/", matrices)]
ICED<-matrices[grep("/iced/", matrices)]


for (i in 1:length(ICED)) {
        file<-ICED[i]
        filename<-basename(file)
        head(filename)
        matrix <- read.table(file, sep = "\t" ,col.names = c("ID1", "ID2", "value"), stringsAsFactors = F)
        to_bed <-gsub('_iced.matrix','_abs.bed', filename)
        resolution<-basename(dirname(file))
        to_dir <-gsub(paste("_",resolution,'_iced.matrix', sep=''),'', filename)
        bedfile<-paste(HiCPath, to_dir, "raw", resolution, to_bed, sep="/")
        tmp <- read.table(bedfile, sep = "\t", col.names=c("chr","start","end","ID2"))
        tmp2 <- read.table(bedfile, sep = "\t", col.names=c("chr","start","end","ID1"))
        m1 = merge(matrix, tmp, by=c("ID2"))
        m2 = merge(m1, tmp2, by=c("ID1"))
        final_table<-m2[c('chr.y','start.y','end.y','chr.x','start.x','end.x','value')]
        colnames(final_table)<-c("A_chr","A_start","A_end","T_chr","T_start","T_end","value")
        trans<-subset(final_table, final_table$A_chr != final_table$T_chr)
        cis<-subset(final_table, final_table$A_chr == final_table$T_chr)
        dir.create(paste(HiCPath, "iced_interactions", sep="/"))
        out<-paste(HiCPath, "iced_interactions/", sep="/")
        write.table(final_table, paste(out, filename, ".txt", sep=""), quote = FALSE, row.names=F,col.names=T,sep="\t")
        write.table(trans, paste(out, filename, "TRANS.txt", sep=""), quote = FALSE, row.names=F,col.names=T,sep="\t")
        write.table(cis, paste(out, filename, "CIS.txt", sep=""), quote = FALSE, row.names=F,col.names=T,sep="\t")
        trans.hashed<-data.frame(ID= paste(trans$A_chr,trans$A_start,trans$A_end,trans$T_chr,trans$T_start,trans$T_end, sep="_"))
        trans.hashed$value<-trans$value
        cis.hashed<-data.frame(ID= paste(cis$A_chr,cis$A_start,cis$A_end,cis$T_chr,cis$T_start,cis$T_end, sep="_"))
        cis.hashed$value<-cis$value
        write.table(trans.hashed, paste(out, filename, "TRANS.hash.txt", sep=""), quote = FALSE, row.names=F,col.names=T,sep="\t")
        write.table(cis.hashed, paste(out, filename, "CIS.hash.txt", sep=""), quote = FALSE, row.names=F,col.names=T,sep="\t")
}

}
