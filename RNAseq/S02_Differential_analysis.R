#DESeq2 analysis - three groups
##RUVseq: https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html#differential-expression-analysis-with-deseq2'
##DESeq2:https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variance-stabilizing-transformation

if (!requireNamespace("BiocManager"))
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn", "tidyverse", "RUVSeq", "DESeq2", "ggplot2"))


library(RUVSeq)
library(ggplot2)
library(tximport)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(xlsx)
library(pheatmap)
library(cluster)
library(PoiClaClu)


setwd("~/Rtest/SATB2 rescue")
seqdata <- read.delim("~/SATB2 rescue/RNAseqMouseRescueCounts - Kopie.txt", stringsAsFactors = FALSE)
seqdata
sampleinfo <- read.delim("~/SATB2 rescue/Phenotypes_Rescue.txt", stringsAsFactors = TRUE)
countdata <- seqdata[,-1]
rownames (countdata) <- seqdata [,1]
head(countdata)
colnames(countdata)


filter <- apply(countdata, 1, function(x) length(x[x>5])>=2)
filtered <- countdata[filter,]
dim(filtered)
genes<-rownames(filtered)
x <- as.factor(c(rep("cKO",4), rep("flx",4), rep("SATB2",4)))
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(x, row.names=colnames(filtered)))

## Boxplots of relative log expression (RLE = log-ratio of read count to median readcount across sample) and PCA plots of raw data

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=0.8)


## RUVg - Empirical control genes
## If no genes are known a priori not to be influenced by the covariates of interest, one can obtain a set of "in-silico empirical" negative controls, e.g., least significantly DE genes based on a first-pass DE analysis performed prior to RUVg normalization.
## k is the number of factors that affect variability

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set, empirical, k=5)
pData(set2)
set2
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=0.8)
head(normCounts(set2))
write.table(normCounts(set2), file = "Norm counts_RUVg_k5.txt", sep = "\t")

#The RUVg function returns two pieces of information: the estimated factors of unwanted variation (added as columns to the phenoData slot of set) and the normalized counts obtained by regressing the original counts on the unwanted factors. 
#The normalized values are stored in the normalizedCounts slot of set and can be accessed with the normCounts method. 
#These counts should be used only for exploration. 
#It is important that subsequent DE analysis be done on the original counts (accessible through the counts method), 
#as removing the unwanted factors from the counts can also remove part of a factor of interest (Gagnon-Bartsch, Jacob, and Speed 2013).




dds <- DESeqDataSetFromMatrix(countData = counts(set2),colData = pData(set2), design = ~ W_1 + W_2 + W_3 + W_4 + W_5 + x)

dds$x=relevel (dds$x, "flx")

dds
dds <- DESeq(dds)
resultsNames(dds)

res = results(dds, contrast=c("x","cKO","flx"))
res
sum(res$padj < 0.05, na.rm=TRUE)
resSig_cKOvsFlx <- subset(res, padj < 0.05)
resSig_cKOvsFlx
resLFC_cKOvsFlx <- subset(resSig_cKOvsFlx, abs(log2FoldChange) > 0.3)
nrow(resLFC_cKOvsFlx)
write.table(resLFC, file = "DEGs_cKOvsFlx_RUVg_k5.txt", sep = "\t")


res = results(dds, contrast=c("x","cKO", "SATB2"))
sum(res$padj < 0.05, na.rm=TRUE)
resSig_cKOvsRescue <- subset(res, padj < 0.05)
resLFC_cKOvsRescue <- subset(resSig_cKOvsRescue, abs(log2FoldChange) > 0.3)
nrow(resLFC_cKOvsRescue)
resLFC_cKOvsRescue
write.table(resLFC_cKOvsRescue, file = "DEGs_cKOvsRescue_RUVg_k5.txt", sep = "\t")

res = results(dds, contrast=c("x","SATB2", "flx"))
sum(res$padj < 0.05, na.rm=TRUE)
resSig_SATB2vsFlx <- subset(res, padj < 0.05)
resLFC_SATB2vsFlx <- subset(resSig_SATB2vsFlx, abs(log2FoldChange) > 0.3)
nrow(resLFC_SATB2vsFlx)
resLFC_SATB2vsFlx
write.table(resLFC_SATB2vsFlx, file = "DEGs_SATB2vsFlx_RUVg_k5.txt", sep = "\t")

##Gene clustering
#Since the clustering is only relevant for genes that actually carry a signal, one usually would only cluster a subset of the most highly variable genes. 
#Here, for demonstration, let us select the 50 genes with the highest variance across samples. We will work with the VST data (generated using RUV-seq-corrected counts).
#The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
#Hence, we center each genes’ values across samples, and plot a heatmap 

library("genefilter")
dds_r <- DESeqDataSetFromMatrix(countData = normCounts(set2),colData = pData(set2), design = ~ x)
vsd_r <- vst(dds_r, blind = FALSE)
topVarGenes <- head(order(rowVars(assay(vsd_r)), decreasing = TRUE), 50)

mat  <- assay(vsd_r)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd_r)[c("x")])
pheatmap(mat, annotation_col = anno,fontsize_row = 7)

topVarGenes <- head(order(rowVars(assay(vsd_r)), decreasing = TRUE), 100)

mat  <- assay(vsd_r)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd_r)[c("x")])
pheatmap(mat, annotation_col = anno,fontsize_row = 7)



## In order to test for differential expression, we operate on raw counts and use discrete distributions.
## However, for other downstream analyses - e.g. for visualization or clustering - it might be useful to work with transformed versions of the count data.
## Visualizations using raw counts (not RUVseq corrected!!!)

nrow(dds)

library(pheatmap)
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[c("x")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
sampleDists <- dist(t(assay(ntd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$x, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(ntd, intgroup=c("x"))


vsd <- vst(dds, blind=FALSE)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[c("x")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$x, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("x"))





rld <- rlog(dds, blind = FALSE)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$x, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(rld, intgroup=c("x"))

pcaData <- plotPCA(rld, intgroup=c("x"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=x, shape=x)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$x)
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$x <- dds$x
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = x, shape = x)) +
geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


##Heatmaps are a nice visualisation to examine hierarchical clustering of your samples. 
##We can do this using the heatmap.2 function from the gplots package. 
##In this example heatmap.2 calculates a matrix of euclidean distances from the logCPM (counts(set2) object) for the 500 most variable genes.

library(gplots)

var_genes <- apply(counts(set2), 1, var)
head(var_genes)
var_genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:2000]
head(select_var)
highly_variable_lcpm <- counts(set2)[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange", "green")[x]
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 2000 most variable genes across samples",
          ColSideColors=col.cell,scale="row")


library("pheatmap")
library("cluster")

rld <- rlog(dds, blind = FALSE)
allSig_genes <- rownames(resSig_cKOvsFlx)
allSig_genes
select <- order(rowVars(assay(vsd)[allSig_genes,]), decreasing = TRUE) [1:500]
select

df <- as.data.frame(colData(dds)[c("x")])
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)


rld <- rlog(dds, blind = FALSE)
allSig_genes <- rownames(resSig_cKOvsFlx)
allSig_genes
df <- as.data.frame(colData(dds)[c("x")])
pheatmap(assay(rld)[allSig_genes,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)


sampleDists <- dist(t(assay(rld)[allSig_genes,]))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$x, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA((rld)[allSig_genes,], intgroup=c("x"))




allSig_genes <- rownames(resSig_cKOvsFlx)
topVarGenes <- head(order(rowVars(assay(rld)[allSig_genes,]), decreasing = TRUE), 2819)
set.seed(1574)
k <-   pheatmap(assay(rld)[topVarGenes,], scale="row",kmeans_k = 2)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"

OrderByCluster <- assay(rld)[topVarGenes,][order(clusterDF$Cluster),]

pheatmap(OrderByCluster,
         scale="row",annotation_row = clusterDF,
         show_rownames = FALSE,cluster_rows = FALSE)




allSig_genes <- rownames(resSig_cKOvsFlx)
k <-   pheatmap(assay(rld)[allSig_genes,], scale="row",kmeans_k = 2)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"

OrderByCluster <- assay(rld)[allSig_genes,][order(clusterDF$Cluster),]
pheatmap(OrderByCluster,
         scale="row",annotation_row = clusterDF,
         show_rownames = FALSE,cluster_rows = FALSE, cluster_cols=FALSE)



pheatmap(OrderByCluster,
         scale="row",annotation_row = clusterDF,
         show_rownames = FALSE,cluster_rows = FALSE)

## In order to test for differential expression, we operate on raw counts and use discrete distributions.
## However for other downstream analyses - e.g. for visualization or clustering - it might be useful to work with transformed versions of the count data.
## Generating transformed RUVseq-corrected counts!!!:

## Useful first step in an RNA-seq analysis is often to assess overall similarity between samples: 
## Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment’s design?
## We use the R function dist to calculate the Euclidean distance between samples.
## To ensure we have a roughly equal contribution from all genes, we use it on the VST data (or rlog data). 
## We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

library("pheatmap")
dds_r <- DESeqDataSetFromMatrix(countData = normCounts(set2),colData = pData(set2), design = ~ x)
rld_r <- rlog(dds_r, blind = FALSE)
sampleDists <- dist(t(assay(rld_r)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld_r$x, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(rld_r, intgroup=c("x"))


dds_r <- DESeqDataSetFromMatrix(countData = normCounts(set2),colData = pData(set2), design = ~ x)
vsd_r <- vst(dds_r, blind = FALSE)
sampleDists <- dist(t(assay(vsd_r)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_r$x, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

##Heatmap of sample-to-sample distances using the variance stabilizing transformed values.
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

##Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. 
##This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples. 
##The PoissonDistance function takes the original count matrix (not normalized, RUVseq-corrected!!) with samples as rows instead of columns, so we need to transpose the counts in dds.

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds_r)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds_r$x)
colnames(samplePoisDistMatrix) <- NULL
## Heatmap of sample-to-sample distances using the Poisson Distance
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


#Another way to visualize sample-to-sample distances is a principal components analysis (PCA). 
#In this ordination method, the data points (here, the samples) are projected onto the 2D plane such that they spread out in the two directions that explain most of the differences (figure below).
#The x-axis is the direction that separates the data points the most. 
#The values of the samples in this direction are written PC1. 
#The y-axis is a direction (it must be orthogonal to the first direction) that separates the data the second most. 
#The values of the samples in this direction are written PC2. 
#The percent of the total variance that is contained in the direction is printed in the axis label. 
#Note that these percentages do not add to 100%, because there are more dimensions that contain the remaining variance (although each of these remaining dimensions will explain less than the two that we see).

#PCA plot using the VST data:

plotPCA(vsd_r, intgroup=c("x"))



pcaData <- plotPCA(rld_r, intgroup=c("x"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=x, shape=x)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#Another technique for performing dimension reduction on data that is not Normally distributed (e.g. over-dispersed count data) is
#generalized principal component analysis, or GLM-PCA, (Townes et al. 2019) as implemented in the CRAN package glmpca. 
#This package takes as input the count matrix, as well as the number of latent dimensions to fit (here, we specify 2). 
#As stated by Townes et al. (2019):"…we propose the use of GLM-PCA, a generalization of PCA to exponential family likelihoods. 
#GLM-PCA operates on raw counts, avoiding the pitfalls of normalization. 
#We also demonstrate that applying PCA to deviance or Pearson residuals provides a useful and fast approximation to GLM-PCA."

library("glmpca")
gpca <- glmpca(counts(dds_r), L=2)
gpca.dat <- gpca$factors
gpca.dat$x <- dds_r$x
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = x, shape = x)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


#Another plot, very similar to the PCA plot, can be made using the multidimensional scaling (MDS) function in base R. 
#This is useful when we don’t have a matrix of data, but only a matrix of distances. 
#Here we compute the MDS for the distances calculated from the VST data and plot these in a figure below.

mds <- as.data.frame(colData(vsd_r))  %>%
cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = x, shape = x)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")








##Heatmaps are a nice visualisation to examine hierarchical clustering of your samples. 
##We can do this using the heatmap.2 function from the gplots package. 
##In this example heatmap.2 calculates a matrix of euclidean distances from the logCPM normCounts(set2) for the 500 most variable genes.

library(gplots)

var_genes <- apply(normCounts(set2), 1, var)
head(var_genes)
var_genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:2000]
head(select_var)
highly_variable_lcpm <- (normCounts(set2))[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange", "green")[x]
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 2000 most variable genes across samples",
          ColSideColors=col.cell,scale="row")


## Clustering, PCA, and heatmap based on DEGs (cKO vs Flx)

library(pheatmap)
library(cluster)

allSig_genes <- rownames(resSig_cKOvsFlx)
allSig_genes
df <- as.data.frame(colData(dds_r)[c("x")])
pheatmap(assay(rld_r)[allSig_genes,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

sampleDists <- dist(t(assay(rld_r)[allSig_genes,]))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld_r$x, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA((rld_r)[allSig_genes,], intgroup=c("x"))


allSig_genes <- rownames(resSig_cKOvsFlx)
k <-   pheatmap(assay(rld_r)[allSig_genes,], scale="row",kmeans_k = 2)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- assay(rld)[allSig_genes,][order(clusterDF$Cluster),]
pheatmap(OrderByCluster,
         scale="row",annotation_row = clusterDF,
         show_rownames = FALSE,cluster_rows = FALSE)






