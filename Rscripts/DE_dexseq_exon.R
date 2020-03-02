#================================================================
#- IMPLEMENTATION
#-    version         0.0.1
#-    author          Richa Bharti
#-    copyright       Copyright (c) 2020
#-    license         GNU General Public License
#================================================================



#check and install libraries
if("DEXSeq" %in% rownames(installed.packages()) == FALSE) {BiocManager::install('DEXSeq')}
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}
if("gplots" %in% rownames(installed.packages()) == FALSE) {install.packages("gplots")}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("argparse" %in% rownames(installed.packages()) == FALSE) {install.packages("argparse")}

library('argparse')

p <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
p$add_argument("--exoncountpath", help="Enter the dexseq generated exon count file path")
p$add_argument("--gffFilePath", help="Enter the dexseq generated gff file path")
p$add_argument("--metadata", help="Enter the metadata file containing samples information")
p$add_argument("--outputpath", help="Enter the output folder path")

args <- p$parse_args()
#args <- commandArgs(trailingOnly = TRUE)
countFiles<- list.files(args$exoncountpath, pattern=".csv$", full.names=TRUE) 
metadataFile =  args$metadata
outputFolder <- args$outputpath
flattenedFile = list.files(args$gffFile, pattern="gff$", full.names=TRUE)
#plotGene

#loading libraries
suppressPackageStartupMessages({
	library('DEXSeq')
	library('RColorBrewer')
	library('gplots')
	library('dplyr')
})


metadata = read.table(metadataFile, sep='\t', header=TRUE,check.names = FALSE)
row.names(metadata) = metadata$Sample 
metadata$Sample = NULL

dxd = DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData=metadata,
    design= ~ sample + exon + condition:exon,
    flattenedfile=flattenedFile)

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )

dxd = testForDEU( dxd )

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")


dxr1 = DEXSeqResults( dxd )

## Write results
write.csv(dxr1, file=paste(outputFolder,"diffexpr-results_dexseq.csv",sep=""))



# Plot MA 
pdf(paste(outputFolder,"ma_plot_dexseq.pdf",sep=""))#, w=1000, h=1000, pointsize=20)
plotMA( dxr1, cex=0.8, main="MA plot" )
dev.off()

# Plot dispersions
pdf(paste(outputFolder,"dispersion_plot_dexseq.pdf",sep=""))#, w=1000, h=1000, pointsize=20)
plotDispEsts( dxd, main="Dispersion plot" )
dev.off()

# Plot dispersions
pdf(paste(outputFolder,"fitted_expression_plot_dexseq.pdf",sep=""))#, w=1000, h=1000, pointsize=20)
plotDEXSeq( dxr1, "ENSMUSG00000038738.15", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()
#quit("no")
