#' Compartment computation and saddle plots visualization using GENOVA
#' ATACseq peaks consensus file computed using Genrich were used to discern active (A) versus B compartments. 

library(GENOVA)
library(stringr)

Dir<-"/hic_results/matrix/"
ATACseqConsensus<-"/ATAC-seq/0290821.ATACseq.Genrich.condition/consensusGenrich.narrowPeak"

# Processing data:
flxBIC_100kb <- load_contacts(
  signal_path = paste0(Dir, 'flx_BIC/iced/100000/flx_BIC_100000_iced.matrix'),
  indices_path = paste0(Dir, 'flx_BIC/raw/100000/flx_BIC_100000_abs.bed'),
  sample_name = "flxBIC",
  colour = "black"
)

# Loading data 
ckoBIC_100kb <- load_contacts(
  signal_path = paste0(Dir, 'cKO_BIC/iced/100000/cKO_BIC_100000_iced.matrix'),
  indices_path = paste0(Dir, 'cKO_BIC/raw/100000/cKO_BIC_100000_abs.bed'),
  sample_name = "ckoBIC",
  colour = "red"
)

ATACseqConsensusPeaks = read.delim(ATACseqConsensus, header = FALSE)

CS_out = compartment_score(list(flxBIC_100kb, ckoBIC_100kb),
                           bed = ATACseqConsensusPeaks)

saddle_out = saddle(list(flxBIC_100kb, ckoBIC_100kb),
                   CS_discovery = CS_out, bins = 50)

save(CS_out, 
     saddle_out, file=paste0(Dir, "SaddleCompartmentData.rda"))


visualise(saddle_out)


CSS <- quantify(saddle_out)
compared <- tidyr::spread(unique(CSS[,-c(3,4)]), key = 'exp', value = 'strength'
with(compared, plot(WT, WAPL, xlim = c(0,4), ylim = c(0,4), pch = 20))
abline(a = 0, b = 1, lty = 1)
