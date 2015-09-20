#!/usr/bin/Rscript

# usage: $ Rscript DE_analysis_ENSG.R read_counts_matrix

args <- commandArgs(trailingOnly = TRUE)

# Read data file
dataRNAseq = read.table(args[1], header = TRUE, row.names = 1, sep="\t")

pdf("DE_analysis_graphics.pdf")

# Calculate logFC values using read counts

# mean values for melanocytes and cancerous cells
meanMcounts = apply(dataRNAseq[,1:2],1,mean)
meanCcounts = apply(dataRNAseq[,3:4],1,mean)

# logFC on raw data
logFC = log2((meanCcounts + 1)/(meanMcounts + 1))

# distribution of logFC on raw data
hist(logFC, nclass = 100, main = "logFC(cancerous/melanocytes) \n(raw data)", 
xlab = "log(cancerous/melanocytes) value")
abline(v = 0, col = "red")


# DESeq package

library(DESeq2)

# Loading data for the experiment
# M = "normal" melanocyte
# C = cancerous cell
# design.txt = text file with 2 columns, first experiment and second condition (M/C)

colData = read.table("design.txt", row.names = 1, header = TRUE)

# DESeqDataSet object creation
dds = DESeqDataSetFromMatrix(countData = dataRNAseq[,1:4], colData = colData, design = ~condition)
#nrow(dds)
#60234

# Pre-filtering the data set (removing rows with no counts or a single count)
dds = dds[rowSums(counts(dds))>1,]
#nrow(dds)
#47451

# calculation of sizeFactors
dds = estimateSizeFactors(dds)
sizeFactors(dds)

# Visual exploring of the data 

# rlog transformation (regularized log transforlation, stabilize variance across the mean)
# for fully unsupervised transformation, set blind=TRUE
rld = rlog(dds,blind=TRUE)

# Effect of the rlog transformation, first two samples
par(mfrow=c(1,2))
dds=estimateSizeFactors(dds)
plot(log2(counts(dds,normalized=TRUE)[,1:2]+1),pch=16,cex=0.3,
main="Before rlog transformation")
plot(assay(rld)[,1:2],pch=16,cex=0.3,
main="After rlog transformation")

# PCA plot
par(mfrow=c(1,1))
p_rld = plotPCA(rld,intgroup=c("condition"))
p_rld = update(p_rld, panel = function(x, y, ...) {lattice::panel.xyplot(x, y, ...);
lattice::ltext(x=x, y=y, labels=rownames(colData(rld)), pos=1, offset=1, cex=0.5)})
print(p_rld)

# Sample distances
sampleDists = dist(t(assay(rld)))

# Heatmaps distances
library("RColorBrewer")
library("pheatmap")

sampleDistMatrix = as.matrix(sampleDists)
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists, col = colors,
main="Heatmat of sample distances")


# Normalization of the data

# get normalized count values
cdsNorm = counts(dds, normalized = TRUE)

# mean values
meanMcountsNorm = apply(cdsNorm[,1:2], 1, mean)
meanCcountsNorm = apply(cdsNorm[,3:4], 1, mean)
# sd values for log(H/N) replicates
sdMcountsNorm   = apply(cdsNorm[,1:2], 1, sd)
sdCcountsNorm   = apply(cdsNorm[,3:4], 1, sd)

# logFC (after normalization)
logFCNorm = log2((meanCcountsNorm + 1)/(meanMcountsNorm + 1))

hist(logFCNorm, nclass = 100, main = "logFC (C/M) distribution \n(normalized data)",
xlab = "log(C/M) value")
abline(v = 0, col = "red")

# thresold can be chosen (here the values are 2 and 5) to select up and down regulated genes  
abline(v = 2, col = "red", lty = "dashed")
abline(v = -2, col = "green", lty = "dashed")
abline(v = 5, col = "red", lty = "dashed")
abline(v = -5, col = "green", lty = "dashed")

upGenes2 = names(logFCNorm[logFCNorm > 2])
downGenes2 = names(logFCNorm[logFCNorm < -2])
upGenes5 = names(logFCNorm[logFCNorm > 5])
downGenes5 = names(logFCNorm[logFCNorm < -5])

# evaluate expression level of genes
exprLevel = apply(cdsNorm, 1, mean)

# logFC versus the level of gene expression
plot(log(exprLevel), logFCNorm, pch = 20,
	xlab = "Gene expression level (log scale)", ylab = "logFC",
	main = "RNAseq data")
abline(h = 2, col = "green", lty = "dashed")
abline(h = -2, col = "red", lty = "dashed")
points(log(exprLevel[upGenes2]), logFCNorm[upGenes2], pch = 20,
	col = "green")
points(log(exprLevel[downGenes2]), logFCNorm[downGenes2], pch = 20,
	col = "red")
	
plot(log(exprLevel), logFCNorm, pch = 20,
	xlab = "Gene expression level (log scale)", ylab = "logFC",
	main = "RNAseq data")
abline(h = 5, col = "green", lty = "dashed")
abline(h = -5, col = "red", lty = "dashed")
points(log(exprLevel[upGenes5]), logFCNorm[upGenes5], pch = 20,
	col = "green")
points(log(exprLevel[downGenes5]), logFCNorm[downGenes5], pch = 20,
	col = "red")


######
# Perform the DE analysis with DESeq
######

## Differential analysis
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
res = results(dds)
mcols(res,use.names=TRUE)
	
hist(-log(res$padj), breaks = 20, col = "black", border="white", 
	xlab = "-log(p-value)", 
	main = "Distribution of -log(adjusted pvalues)")

# writing of the results
write.table(res, "DESeq2_statistics.txt", row.name=T, quote=F, sep='\t')
write.table(upGenes2, "up_genes_2.txt", row.name=F, col.name=F, quote=F)
write.table(downGenes2, "down_genes_2.txt", row.name=F, col.name=F, quote=F)
write.table(upGenes5, "up_genes_5.txt", row.name=F, col.name=F, quote=F)
write.table(downGenes5, "down_genes_5.txt", row.name=F, col.name=F, quote=F)

dev.off()

topGenes = head(order(res$padj),100)
write.table(res[topGenes,],"results_DESeq_100topGenes.txt",sep="\t",quote=F,row.name=T)

#table(res$padj<0.05)
#FALSE  TRUE 
#17019 23315 

# lower FDR threshold to 5% (default: 10%)
res.05 = results(dds,alpha=0.05)
#table(res.05$padj<0.05)
# idem
res.01 = results(dds,alpha=0.01)
#table(res.01$padj<0.05)
#FALSE  TRUE 
#17019 23315 

# raise logFC threshold
res.FC2 = results(dds,lfcThreshold=2)

res.FC5 = results(dds,lfcThreshold=5)


# plotMA topGene in graphics
pdf("plotMA_resFC2_topGene.pdf")
plotMA(res.FC2,ylim=c(-15,15))
topGene = rownames(res.FC2)[which.min(res.FC2$padj)]
with(res[topGene,], {
	points(baseMean,log2FoldChange,col="black",cex=2,lwd=2)
	text(baseMean,log2FoldChange,topGene,pos=2,col="black")
	})
dev.off()

pdf("plotMA_resFC5_topGene.pdf")
plotMA(res.FC5,ylim=c(-15,15))
topGene_LC1 = rownames(res.FC5)[which.min(res.FC5$padj)]
with(res[topGene,], {
	points(baseMean,log2FoldChange,col="black",cex=2,lwd=2)
	text(baseMean,log2FoldChange,topGene,pos=2,col="black")
	})
dev.off()
