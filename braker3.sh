#!/bin/bash

date  #script start time

merged_bam="$1"         #full path to merged BAM file
genome="$2"             #path to genome FASTA
stringtie_dir=$3        #merged gtf file directory
output_directory="$4"   #output directory

#argument check
if [ -z "$merged_bam" ] || [ -z "$genome" ] || [ -z $stringtie_dir ] || [ -z "$output_directory" ]; then
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

#genemark environment check
export GENEMARK_PATH="/home/users/ja1473/gmes_linux_64_4"
export GM_KEY="$GENEMARK_PATH/gm_key"

if [ ! -f "$GM_KEY" ]; then
    echo "Missing GeneMark key file"
    exit 1
fi

#run BRAKER2
braker.pl \
    --species=casuarina_glauca \
    --etmode \
    --genome="$genome" \
    --bam="$merged_bam" \
    --softmasking \
    --workingdir="$output_directory" \
    --cores=26 \
    --UTR=on

echo "Finished annotating Casuarina glauca"
date  #time to finish script
