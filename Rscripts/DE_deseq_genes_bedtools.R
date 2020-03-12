#================================================================
#- IMPLEMENTATION
#-    version         0.0.1
#-    author          Richa Bharti
#-    copyright       Copyright (c) 2020
#-    license         GNU General Public License
#================================================================



#check and install libraries
if("DESeq2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install('DESeq2')}
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}
if("gplots" %in% rownames(installed.packages()) == FALSE) {install.packages("gplots")}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("argparse" %in% rownames(installed.packages()) == FALSE) {install.packages("argparse")}

library('argparse')

p <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
p$add_argument("--genecount", help="Enter the bedtools generated gene count file")
p$add_argument("--metadata", help="Enter the metadata file containing samples information")
p$add_argument("--condition",  help="Enter the column name in metadata file used for pairwise comparison")
p$add_argument("--outputpath", help="Enter the output folder path")

args <- p$parse_args()
#args <- commandArgs(trailingOnly = TRUE)
geneCountingFile <- args$genecount
metadataFile =  args$metadata
condition_Metadata = args$condition
outputFolder <- args$outputpath

#loading libraries
suppressPackageStartupMessages({
	library('DESeq2')
	library('RColorBrewer')
	library('gplots')
	library('dplyr')
})


#print (geneCountingFile)
#print (metadataFile)
#print (condition_Metadata)
#print (outputFolder)


rawCountTable <- read.table(geneCountingFile, sep='\t', header=TRUE,check.names = FALSE)
rownames(rawCountTable) = rawCountTable$Attributes
rawCountTable$Attributes = NULL
countTable <- round(rawCountTable[,9:length(names(rawCountTable))])

countTable <- as.matrix(countTable)
libs <- colnames(countTable)
metadata = read.table(metadataFile, sep='\t', header=TRUE,check.names = FALSE)
condition = metadata[condition_Metadata]
comp1 = (unique(condition))[1,1]
comp2 = (unique(condition))[2,1]
comp = paste(comp1,"_vs_",comp2,sep="")


dds <- DESeqDataSetFromMatrix(countData=countTable, colData=metadata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
pdf(paste(outputFolder,"qc-dispersions_bedtools.pdf",sep=""))#, w=1000, h=1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
hist(assay(rld))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
pdf(paste(outputFolder,"qc-heatmap-samples_bedtools.pdf",sep=""))#, w=1000, h=1000, pointsize=20)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- with(colData(dds), paste(libs, sep=' : '))
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(13, 13))
dev.off()

# Principal components analysis
pdf(paste(outputFolder,"PCA-samples_bedtools.pdf",sep=""))#, w=1000, h=1000, pointsize=20)
rld <- rlog(dds)
print(plotPCA(rld, intgroup=c('condition')))
dev.off()


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Attributes"
head(resdata)
## Write results
write.table(resdata,sep="\t",row.names=FALSE,quote=F, file=paste(outputFolder,comp,"_diffexpr-results_bedtools.txt",sep=""))

## Examine plot of p-values
#hist(res$pvalue, breaks=50, col="grey")

## Examine plot of p-values
pdf(paste(outputFolder,"pvalue_histogram_50_bedtools.pdf", sep=""))#, 1500, 1000, pointsize=20)
hist(res$pvalue, breaks=50, col="grey")
dev.off()

## Examine independent filtering
#attr(res, "filterThreshold")
#plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:

pdf(paste(outputFolder,"diffexpr-maplot_bedtools.pdf", sep=""))#, 1500, 1000, pointsize=20)
plotMA(res, main="MA Plot")
dev.off()
