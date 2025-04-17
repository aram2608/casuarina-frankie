#!/bin/bash

input_seq=$1 # full path to sequences 
$data_base_dir=$2 #full path to database, example = ~/path/to/swissprot ---- no extension needed

# diamond params

./diamond blastp \
    -d swissprot \
    -q $input_seq \
    --outfmt 6 qseqid sseqid qlen slen qstart qend sstart send length pident positive mismatch gaps evalue bitscore qcovs \
    -p 12

# downloading the tool
# wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz
# tar xzf diamond-linux64.tar.gz

# downloading and using a BLAST database
#update_blastdb.pl --decompress swissprot
# ./diamond prepdb -d swissprot

# creating a diamond-formatted database file
# ./diamond makedb --in reference.fasta -d reference