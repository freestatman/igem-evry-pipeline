##########################################################
## MHC-I affinity predictions and final antigen scoring ##
##########################################################

1) USAGE

Open a terminal window and run (script should be made executable):

$ ./mhcI_affinity.R -lt 500 -th 50 -l 9 input_file.fasta output_file.tsv

The user can change thresholds for low-binding (-lt) and high-binding (-th) scores, antigen length (-l). 

Help about command-line arguments can be printed by running: $ ./mhcI_affinity.R --help


2) INPUT

The input file must be a output file of 'proteasome_cleavage.R' (predicted antigens table)


3) OUTPUT

The script generates a table containing the results of the predictions, structured as follow: 
- 'antigene': this column contains the potential antigen sequence
- 'insertGene': this column contans the sequence to insert in our system (it must be longer because it needs to be cleaved by the proteasome) 
- 'netchopScore': Netchop score associated with the antigen. 
- 'nmpScore': predicted MHC-I affinity score by NetMHCpan. The predicted output is given in units of IC50nM, therefore a lower number indicates higher affinity. It is generally assumed that peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity, hence the default parameters values.
- 'hlaAllele': the HLA allele of CMH-I for which the affinity socre was the best
- 'combinedScore': the final score associated to the potential antigen (sum of Netchop socre and normalized NetMHCpan score)

	
4) SOFTWARE / DEPENDENCIES NEEDED

You should install the following (latest versions):

- R version > 3.1.1
- R package 'argparse'
- NetMHCpan (free for academic use - http://www.cbs.dtu.dk/services/NetMHCpan/)
