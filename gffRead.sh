#!/bin/bash

braker_input=$1 # input to braker gtf ----- /full/path/to/file
genome_input=$2 # input to clean genome ------- /full/path/to/file
protein_output=$(basename $braker_input .gtf).faa # extracts name for protein fasta
transcript_output=$(basename $braker_input .gtf).fna # extracts name for transcript fasta

# safety check to ensure all args are met
if [ ! -f "$braker_input" ] || [ ! -f "$genome_input" ]; then
    echo "Usage: <braker_input> <genome_input>"
    exit 1
fi

# gffread input to extract fastas from the genome and braker output
gffread "$braker_input" \
    -g $genome_input \
    -y $protein_output \
    -w $transcript_output

echo "Finished extracting sequences for $braker_input"