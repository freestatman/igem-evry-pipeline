######################################################
### Variant Filtering                              ###
### Hard filters -> optimize both high sensitivity ### 
### and specificity together ;                     ###
### !!! some real sites will get filtered out !!!  ###
######################################################

# Compute basic statistics (num variants, quality, ...) with plinkseq
#pseq  melanocytes_melanomes_var.vcf v-stats > stats_raw.out
#pseq  melanocytes_melanomes_var.vcf i-stats > stats_raw_individual.out

#1) Keep diallelic variants only
vcftools --vcf melanocytes_melanomes_var.vcf --min-alleles 2 --max-alleles 2 --recode --out N_C_diallelic

#2) Annotation, beforeQC
java -Xmx32g -jar snpEff.jar -v GRCh37.75 N_C_diallelic.recode.vcf > N_C_diallelic_annot.vcf

# annotate unknown variants only (unknown as not reported in dbSNP)
java -jar /SnpSift.jar annotate -dbsnp N_C_diallelic.recode.vcf > bQC_dbsnp.vcf
java -Xmx4g -jar snpEff.jar eff -v GRCh37.75 bQC_dbsnp.vcf > bQC_eff.vcf
java -jar SnpSift.jar filter -f bQC_eff.vcf "! exists ID" > bQC_eff_not_in_dbSnp.vcf
java -Xmx32g -jar snpEff.jar eff -v GRCh37.75 bQC_eff_not_in_dbSnp.vcf > bQC_not_in_db_annot.vcf

#pseq bQC_not_in_db_annot.vcf v-stats > stats_novel.out

#3) High Quality variants (CR>98% and HWE p > 10-7)
vcftools --vcf N_C_diallelic.recode.vcf --hwe 0.0000001 --recode --out N_C_HF_hwe
vcftools --vcf N_C_HF_hwe.recode.vcf --max-missing 0.98 --recode --out N_C_CR98

#pseq N_C_CR98.recode.vcf v-stats  > stats_hard_filt_CR98.out

#4) GATK filters (as in BEST PRACTICES for RNAseq data and variant calling)

# Filtering based on Fisher Strand values and Qual by Depth 
# Filter out clusters of at least 3 SNPs in a window of 35 bases between them
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V N_C_CR98.recode.vcf \
 -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD \
  -filter "QD < 2.0" -o afterQC_variants.vcf

#5) Annotation, afterQC
java -Xmx32g -jar snpEff.jar -v GRCh37.75 afterQC_variants.vcf > afterQC_variants_annot.vcf

#pseq afterQC_variants_annot.vcf v-stats  > afterQC_variants.out
