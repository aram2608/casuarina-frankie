#!/bin/bash

date #prints start time
sam_file_directory=$1

#tests for arguments
if [ -z $sam_file_directory ]; then
    echo "Usage: ./sambambraker2 <sam_file_directory>"
    exit 1 #safely exists program
fi

#loops through sam directory to perform conversion and braker2 annotation
for sam_file in $sam_file_directory/*.sam.gz; do
    echo Converting SAM to BAM for $sam_file

    base_name=$(basename $sam_file .sam.gz) #allows removal of suffixes
    sorted_bam=$sam_file_directory/sorted_${base_name}.bam #creates output file name for sorted BAMs
 
    #samtool conversion and failure checking
    if gunzip -c $sam_file | samtools view -@ 20 -Sb | samtools sort -O bam -o $sorted_bam; then
        echo $sam_file processed and replaced with sorted BAM
        rm $sam_file
    else
        echo failure processing $sam_file >> failurelog.txt
    fi
done

date #tells you the end time