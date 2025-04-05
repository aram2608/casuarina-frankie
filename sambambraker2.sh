#!/bin/bash

date #prints start time
sam_file_directory=$1

#tests for arguments
if [ -z $sam_file_directory ];then
    echo "Usage: ./sambambraker2 <sam_file_directory>"
    exit 1 #safely exists program
fi
#loops through sam directory to perform conversion and braker2 annotation
for sam_file in $sam_file_directory/*.sam; do
    echo Converting SAM to BAM
    base_name=$(basename $sam_file)
    output_bam=$sam_file_directory/${base_name%.sam}.bam #replaces extension
    sorted_bam=sorted_$output_bam
 
    #samtool conversion
    samtools view -@ 20 -Sb $output_bam $sam_file | samtools sort -O bam -o $sorted_bam $output_bam
    rm $sam_file
    echo "$sam_file processed and replaced with sorted BAM"

    if [ -? -ne 0 ]; then #checks for command success and logs failures
        echo "Warning: Failed to convert SAM file" >> failurelog.txt
    fi
done