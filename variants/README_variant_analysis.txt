##########################################################
#####			      Variant calls 				 #####
##### 		   Analysis - Association tests		     #####
##########################################################

1) USAGE

Open a terminal window and run (script should be made executable):

$ ./VAT_asso.sh

or simply run each of the commands appearing in the script to perform the process step by step and check results (recommended)

!!! We defined common variants as having a Minor Allele Frequency (MAF) > 25% as we had few samples,
this can be easy changed back to MAF > 5% by uncommenting lines 44 and 45.


2) INPUT

The working directory should contain the VCF file produced after quality control (afterQC_variants.vcf) ;
as well as a "pheno.txt" tab separated file with the following format (e.g.):

sample_name	phenotype
SRR1993908_bowtie_mapping.sorted.bam	1
SRR1993909_bowtie_mapping.sorted.bam	1
SRR1993910_bowtie_mapping.sorted.bam	2
SRR1993911_bowtie_mapping.sorted.bam	2

(in here "1" -> normal cells and "2" -> tumoral cells)


3) OUTPUTS

The script sets-up a vtools project (SQL-like project, organised in tables) ; the different tables
can be easily accessed from a command line:

$ vtools show tables

or, for a particular table:

$ vtools show <table_name>

The scripts also automatically generates a tab-separated file with results (for each variant) of Fisher's exact test
(13 columns: chromosome ; position ; reference base ; alternative base ; CDS starting position ; CDS end position ;
mutation type ; region type ; number of variant alleles in phenotype 2 samples ; number of variant alleles in phenotype 1 samples ;
number of heterozygous variants ; number of homozygous variants ; p-value Fisher's exact)

	
4) SOFTWARE / DEPENDENCIES NEEDED

You should install the following (latest versions):

- Variant Tools: http://sourceforge.net/projects/varianttools/
- ANNOVAR: http://annovar.openbioinformatics.org/en/latest/user-guide/download/
- Python should be also installed