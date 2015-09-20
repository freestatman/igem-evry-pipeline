# A simple command (one-line) to get all protein sequences associated to a list of identifiers in a single .fasta file (here 'output.fa')
# This script uses the RESTful Ensembl API
# Use it with a file containing Ensembl gene identifiers (here 'up_genes.txt')

while IFS='' read -r line || [[ -n "$line" ]]; do curl 'http://rest.ensembl.org/sequence/id/'"$line"'?multiple_sequences=1;type=protein' -H 'Content-type:text/x-fasta'; done < up_genes.txt > output.fa
