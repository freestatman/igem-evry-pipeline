##########################################################
#####			      Variant calls 				 #####
##### 		 Quality control - filtering steps	     #####
##########################################################

1) USAGE

Open a terminal window and run (script should be made executable):

$ ./filtering_variants_stringent.sh

or simply run each of the commands appearing in the script to perform the filtering step by step

The user can change thresholds for the HWE p-value test, call rate minimum, Fisher Strand values and Qual by Depth, 
to make them less strict.


2) INPUT

The working directory should contain the VCF file produced after variant calling ;
(melanocytes_melanomes_var.vcf is the default name, it should be change in line 13 for other files)


3) OUTPUTS

The script generates VCF files for each of the filtering steps:
	- N_C_diallelic.recode.vcf: VCF file, diallelic variants only
	- N_C_diallelic_annot.vcf: VCF file, annotation BEFORE quality control
	- N_C_HF_hwe.recode.vcf: VCF file, variants passing HWE p-value filter
	- N_C_CR98.recode.vcf: VCF file, variants having CR>98%
	- afterQC_variants.vcf: VCF file, variants passing QC steps
	- afterQC_variants_annot.vcf: VCF file, annotation AFTER quality control


	
4) SOFTWARE / DEPENDENCIES NEEDED

You should install the following (latest versions):

- PLINK/SEQ: https://atgu.mgh.harvard.edu/plinkseq2/download.shtml
- VCFTOOLS: http://vcftools.sourceforge.net/downloads.html
- SnpEff: http://snpeff.sourceforge.net/download.html
- SnpSift: http://snpeff.sourceforge.net/SnpSift.html
- GATK: https://www.broadinstitute.org/gatk/download/