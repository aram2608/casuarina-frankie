#!/bin/bash

date  #script start time

merged_bam="$1"         #full path to merged BAM file
genome="$2"             #path to genome FASTA
output_directory="$3"   #output directory

#argument check
if [ -z "$merged_bam" ] || [ -z "$genome" ] || [ -z "$output_directory" ]; then
    echo "Usage: ./braker3.sh <merged_bam> <genome> <output_directory>"
    exit 1
fi

#bam file validation
if [ ! -f "$merged_bam" ]; then
    echo "Merged BAM file not found: $merged_bam"
    exit 1
fi

#bam index validation
if [ ! -f "${merged_bam}.bai" ]; then
    echo "BAM index not found: ${merged_bam}.bai"
    echo "Generating index..."
    samtools index "$merged_bam"
fi

#run BRAKER2
braker.pl \
    --species=casuarina_glauca_rna \
    --genome="$genome" \
    --bam="$merged_bam" \
    --stranded=- \
    --softmasking \
    --workingdir="$output_directory" \
    --gff3 \
    --cores=26 \
    --UTR=on

echo "Finished annotating Casuarina glauca"
date  #time to finish script

#added strandness flag
#added gff3 for evm
#--stranded=- for reverse