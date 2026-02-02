#!/bin/bash

input_seq=$1 # full path to sequences 
data_base=$2 #full path to database, example = ~/path/to/swissprot ---- no extension needed
diamond_path=$3 # path to diamond binaries example = ~/path/to/diamond
output_path=$4 # path to desired output diretory

# safety check for arguments
if [ ! -f $input_seq ] || [ -z $data_base_dir ] || [ -z $diamond_path ] || [ -z $output_path ]; then
    echo "Usage: <input_to_sequences> <database_path> <path_to_diamond_binaries> <output_path>"
    exit 1
fi

# makes output just in case
mkdir -p $output_path

# diamond params for a protein blast search
$diamond_path/./diamond blastp \
    -d "$data_base" \
    -q "$input_seq" \
    --evalue 1e-5 \
    --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
    -p 15 >> "$output_path"/diamond_output.txt

# downloading the tool
# wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz
# tar xzf diamond-linux64.tar.gz

# downloading and using a BLAST database
#update_blastdb.pl --decompress swissprot
# ./diamond prepdb -d swissprot

# creating a diamond-formatted database file
# ./diamond makedb --in reference.fasta -d reference

# ron path
# /home/share/databases/ncbi_nr/nr