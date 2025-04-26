#!/bin/bash

date #prints start time
bam_dir=$1 #input directory for bams
output_dir=$2 #output directory for stringtie

#a check to ensure args are met
if [ -z "$bam_dir" ] || [ -z "$output_dir" ]; then
    echo "Usage: ./stringtie.sh <bam_dir> <output_dir>"
    exit 1
fi

#makes directory just in case
mkdir -p "$output_dir"

#loop for bam files
for bam in "$bam_dir"/*.bam; do
    #file name is crazy so chat said to do this
    sample=$(basename "$bam" | sed 's/\.fastq.*//; s/\.fq.*//; s/\.bam//')
    echo Processing "$sample"

    #stringtie params
    stringtie \
        $bam \
        -p 20 \
        --rf \
        -o "$output_dir/${sample}.gtf" \
        -l "$sample"
    echo "finished processing $bam"
done

date #print end time of script

#added --rf flag for reverse strandness of my reads