##########################################################
###        Proteasome cleavage sites prediction        ###
##########################################################

1) USAGE

Open a terminal window and run (script should be made executable):

$ ./preoteasome_cleavage.R -nt 0.5 input_file.fasta output_file.tsv

The user can change threshold for the proteasome cleavage 

Help about command-line arguments can be printed by running: $ ./proteasome_cleavage.R --help


2) INPUT

The input file must be a fasta containing the sequence of one or more proteins.


3) OUTPUT

The script generates a table containing the results of the predictions, structured as follow: 
- 'antigene': this column contains the potential antigen sequence
- 'insertGene': this column contans the sequence to insert in our system (it must be longer because it needs to be cleaved by the proteasome) 
- 'netchopScore': Netchop score associated with the antigen. 

	
4) SOFTWARE / DEPENDENCIES NEEDED

You should install the following (latest versions):

- R version > 3.1.1
- R package 'argparse'
- Netchop (free for academic use - http://www.cbs.dtu.dk/services/NetChop/)
