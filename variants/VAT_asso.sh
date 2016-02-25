#####################################
##### PIPELINE VARIANT ANALYSIS #####
##### python - variant tools    #####
#####################################

## vtools project set-up ##

# initialize project and import vcf file
vtools init proj
vtools import afterQC_variants.vcf --build hg19 --var_info AA AC AN DP --geno_info DP_geno

# import phenotypes
# pheno.txt is a tab separated file: column 1 = sample_name ; column_2 =
# phenotype N (controls) or C (cancerous cells)
vtools phenotype --from_file pheno.txt

# ANNOVAR annotations
# if necessary, download database
#/path/to/annovar/annotate_variation.pl --downdb refGene /path/to/annovar/humandb/ -build hg19

vtools export variant --output ANNOVAR.input --format ANNOVAR

perl annotate_variation.pl -geneanno ANNOVAR.input -buildver hg19 humandb/

vtools update variant --format ANNOVAR_exonic_variant_function --from_file ANNOVAR.input.exonic_variant_function --var_info mut_type function genename

vtools update variant --format ANNOVAR_variant_function --from_file ANNOVAR.input.variant_function --var_info region_type region_name

# annotation: refGene, dbSNP, dbNSFP and 1K annotations
vtools use refGene
vtools use dbSNP
vtools use dbNSFP

# alternative allele frequency calculations
vtools update variant --from_stat 'total_ie=#(GT)' 'num_ie=#(alt)' 'het_ie=#(het)' 'hom_ie=#(hom)' 'other_ie=#(other)' 'num_var=#(mutGT)'
vtools update variant --set 'af_ie=num_ie/(total_ie * 2.0)'

### fisher's exact
# all variants
vtools update variant --from_stat 'num_gt_case=#(GT)' 'num_var_alleles_case=#(alt)' --samples "phenotype = 2 "
vtools update variant --from_stat 'num_gt_ctrl=#(GT)' 'num_var_alleles_ctrl=#(alt)' --samples "phenotype = 1 "
vtools update variant --set "prop_pval=Fisher_exact(num_var_alleles_case, num_var_alleles_ctrl, 2*num_gt_case, 2*num_gt_ctrl)"

### creating variant subsets
vtools select variant "af_ie > 0.005" -t variants "variant table (MAF>0.5%)"
#vtools select variants "af_ie<=0.05 AND af_ie > 0.005" -t rare_var "rare variants defined as having a MAF≤5%"
#vtools select variants "af_ie > 0.05" -t com_var "common variants defined as having a MAF>5%"
vtools select variants "af_ie<=0.25 AND af_ie > 0.005" -t rare_var "rare variants defined as having a MAF≤25%"
vtools select variants "af_ie > 0.25" -t com_var "common variants defined as having a MAF>25%"

# non-synonymous variants only
vtools select variants "mut_type like 'nonsynonymous%' OR mut_type like 'stoploss%' OR mut_type like 'stopgain%' OR mut_type like 'splicing%' OR mut_type like 'frameshift%' OR mut_type like 'nonframeshift%'" -t fvar
vtools select rare_var "mut_type like 'nonsynonymous%' OR mut_type like 'stoploss%' OR mut_type like 'stopgain%' OR mut_type like 'splicing%' OR mut_type like 'frameshift%' OR mut_type like 'nonframeshift%'" -t rare_fvar "nonsynonymous, stoploss, stopgain, splicing and indel variants selected from table rare_var"
vtools select com_var "mut_type like 'nonsynonymous%' OR mut_type like 'stoploss%' OR mut_type like 'stopgain%' OR mut_type like 'splicing%' OR mut_type like 'frameshift%' OR mut_type like 'nonframeshift%'" -t com_fvar "nonsynonymous, stoploss, stopgain, splicing and indel variants selected from table com_var"

# exonic variants only
vtools select variants "region_type = 'exonic' OR region_type = 'exonic;splicing' OR region_type = 'ncRNA_exonic'" -t exo_var "exonic variants from table variant"
vtools select rare_var "region_type = 'exonic' OR region_type = 'exonic;splicing' OR region_type = 'ncRNA_exonic'" -t exo_RV "exonic variants from table rare_var"
vtools select com_var "region_type = 'exonic' OR region_type = 'exonic;splicing' OR region_type = 'ncRNA_exonic'" -t exo_CV "exonic variants from table comm_var"

vtools select fvar "region_type = 'exonic' OR region_type = 'exonic;splicing' OR region_type = 'ncRNA_exonic'" -t exo_fvar "exonic variants from table fvar"
vtools select rare_fvar "region_type = 'exonic' OR region_type = 'exonic;splicing' OR region_type = 'ncRNA_exonic'" -t exo_fRV "exonic variants from table rare_fvar"
vtools select com_fvar "region_type = 'exonic' OR region_type = 'exonic;splicing' OR region_type = 'ncRNA_exonic'" -t exo_fCV "exonic variants from table com_fvar"



########## Association testing ##########

## COMMON variants

# Fisher's exact test
vtools output com_var \
chr pos ref alt refGene.name2 refGene.cdsStart refGene.cdsEnd refGene.strand \
mut_type region_type num_var_alleles_case num_var_alleles_ctrl het_ie hom_ie prop_pval \
--header CHR POS REF ALT GENE CDS_START CDS_END STRAND \
MUT_TYPE REGION NUM_VAR_ALLELES_C NUM_VAR_ALLELES_N NUM_HTZ NUM_HMZ PVAL_FISHER > pval_CV_fisher.txt

# Logistic regression
vtools associate com_var phenotype \
	       --discard_variants "%(NA)>0.1" \
	       --method "LogitRegBurden --name logReg --alternative 2" \
	       --group_by refGene.name2 \
	       --to_db logReg_CV \
	       -j8 > logReg_CV.txt


## RARE variants

# Burden test on exonic, non-synonymous variants
vtools associate exo_fRV phenotype \
	       --discard_variants "%(NA)>0.1" \
	       --method "BurdenBt --name BurdenTest --alternative 2" \
	       --group_by refGene.name2 \
	       --to_db BurdenTest \
	       -j8 > BurdenTest_exo_nsyn.txt


## COMMON and RARE variants

# SKAT on exonic variants
vtools associate exo_var pheno --discard_samples "%(NA)>0.1" --discard_variants "%(NA)>0.1" -m "SKAT --name skat disease" "SKAT --name skato disease --p_method optimal.adj" "SKAT --name skatLog disease --logistic_weights 0.07 150" "SKAT --name skatoLog disease --p_method optimal.adj --logistic_weights 0.07 150" --group_by refGene.name2 --to_db skat_exo -j8 > skat_exo.txt

# SKAT on exonic non-synonymous variants
vtools associate exo_fvar pheno --discard_samples "%(NA)>0.1" --discard_variants "%(NA)>0.1" -m "SKAT --name skat disease" "SKAT --name skato disease --p_method optimal.adj" "SKAT --name skatLog disease --logistic_weights 0.07 150" "SKAT --name skatoLog disease --p_method optimal.adj --logistic_weights 0.07 150" --group_by refGene.name2 --to_db skat_exof -j8 > skat_exof.txt
