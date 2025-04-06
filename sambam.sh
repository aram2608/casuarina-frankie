#!/bin/bash

date #prints start time
sam_file_directory=$1
genome_directory=$2
masked_repeats_directory=$3
species='casuarina_glauca'

#creates repeats for use in BRAKER2

#tests for arguments
if [ -z $sam_file_directory ] || [ -z $genome_directory ] || [ -z $masked_repeats_directory ];then
    echo "Usage: ./sambambraker2 <sam_file_directory> <masked_repeats_directory>"
    exit 1 #safely exists program
fi

#loops through sam directory to perform conversion and braker2 annotation
for sam_file in $sam_file_directory/*.sam; do
    echo Converting SAM to BAM for $sam_file

    base_name=$(basename $sam_file)
    output_bam=$sam_file_directory/${base_name%.sam}.bam #replaces extension
    sorted_bam=$sam_file_directory/sorted_${base_name%.sam}.bam #creates output file name for sorted BAMs
 
    #samtool conversion
    gunzip -c $sam_file | samtools view -@ 20 -Sb $output_bam $sam_file | samtools sort -O bam -o $sorted_bam $output_bam

    if [ $PIPESTATUS[0] -ne 0 ] || [ $PIPESTATUS[1] -ne 0 ] || [ $PIPESTATUS[3] ]; then #checks for command success and logs failures
        echo "Warning: Failed to convert SAM file" >> failurelog.txt
        continue
    fi

    rm $sam_file
    echo $sam_file processed and replaced with sorted BAM
done