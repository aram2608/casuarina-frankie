#!/bin/bash

date  #script start time

protein_fasta="$1"      #full path to merged BAM file
genome="$2"             #path to protein FASTA
output_directory="$3"   #output directory

#argument check
if [ -z "$protein_fasta" ] || [ -z "$genome" ] || [ -z "$output_directory" ]; then
    echo "Usage: ./braker3.sh <protein_fasta> <genome> <output_directory>"
    exit 1
fi

#bam file validation
if [ ! -f "$protein_fasta" ]; then
    echo "Merged BAM file not found: $protein_fasta"
    exit 1
fi

#run BRAKER2
braker.pl \
    --species=casuarina_glauca_prot \
    --genome="$genome" \
    --prot_seq="$protein_fasta" \
    --softmasking \
    --workingdir="$output_directory" \
    --gff3 \
    --cores=26 \
    --UTR=on

echo "Finished annotating Casuarina glauca"
date  #time to finish script

#added gff3 for evm