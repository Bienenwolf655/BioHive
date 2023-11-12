## Description:

Extracts the promotor sequence (default defined as 1000 bp upstream sequence, can be changed via promotorregion option) for the given genes from the genome fasta file.
It uses in the background the bedtools command line tool. An example output can be found under the /example_data/Promotor.fasta.

# Usage:

python BEDT.py bedfile genelist genomefasta intermediatefolder shscript --promotorregion

## with example data (from the Magnaporthe oryzae MG8):

1. cd /example_data
2. python ../BEDT.py MO.bed ./genelist_subject_x.txt ./MO.fasta ./test42 ../bedtools.sh
