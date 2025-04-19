#!/bin/bash

input_seqs=$1 # full path to input seqs example = ~/path/to/sequences/fasta.fasta
interpro_path=$2 # path to interproscan binaries
output_dir=$3 # path to output directory

if [ -z "$input_seqs" ] || [ "$interpro_path" ]; then
    echo "Usage: <input_sequences> <path_to_interpro_binaries>"

# interpro search for pathways and goterms
$interpro_path/./interproscan.sh -i "$input_seqs" \
    -goterms \
    -pa \
    -cpu 15

echo "Finished annotating $input_seqs"