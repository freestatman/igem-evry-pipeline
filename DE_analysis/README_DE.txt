##########################################################
#####		         RNAseq data    	     #####
#####	     Differential expression analysis        #####
##########################################################

1) USAGE

Open a terminal window and run:

$ Rscript DE_analysis_ENSG.R read_counts_matrix


2) INPUT

read_counts_matrix should be a tab-separated file containing the read counts for each gene and for each sample.
The current version of this script works with a file containing the information for 4 samples:
	- the first 2 columns represent the 2 samples for normal cells
	- the last 2 columns represent the 2 samples for tumoral cells
The script can be easily changed to take larger files as input (change the dataRNAseq[,] column selection)


3) OUTPUTS

The script automatically generates the following files in the user's working directory:
	- DE_analysis_graphics.pdf: pdf file containing graphs for visual inspection and interpretation of results
	- DESeq2_statistics.txt: tab separated file with analysis statistics for each gene / transcript
	- up_genes_2.txt: list of upregulated genes (logFC > 2)
	- up_genes_5.txt: list of upregulated genes (logFC > 5)
	- down_genes_2.txt: list of downregulated genes (logFC < -2)
	- down_genes_5.txt: list of downregulated genes (logFC < -5)
	- results_DESeq_100topGenes.txt: tab separated file with analysis statistics for top 100 genes / transcripts

	
4) SOFTWARE / DEPENDENCIES NEEDED

You should install the following:

- R version > 3.1.1
- R package DESeq2
- R package RColorBrewer
- R package pheatmap

