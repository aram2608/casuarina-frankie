#!/bin/bash

date #prints start time
bam_dir=$1 #input directory for bams
output_dir=$2 #output directory for stringtie

#makes directory just in case
mkdir -p $output_dir

#a check to ensure args are met
if [ -z $bam_dir ] || [ -z $output_dir ]; then
    echo "Usage: ./stringtie.sh <bam_dir> <output_dir>"
    exit 1
fi

#loop for bam files
for bam in $bam_dir/*.bam; do
    $sample=$(basename $bam .bam)
    echo Processing $sample

    #stringtie params
    stringtie \
        -p 20 \
        -o $output_dir/${sample}.gtf \
        -l ${$sample} 
done