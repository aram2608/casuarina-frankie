#!/bin/bash

date #script start time

merged_sorted_bam=$1 #input dir
trinity_output=$2 #output dir

#check for args and bam files
if [ ! -f "$merged_sorted_bam" ] || [ -z "$trinity_output" ]; then
    echo "Usage: ./trinity.sh <bam_dir> <trinity_output>"
    exit 1
fi

#make output just in case
mkdir -p "$trinity_output"

#trinity params
Trinity \
    --max_memory 50g \
    --genome_guided_bam "$merged_sorted_bam" \
    --output "$trinity_output" \
    --SS_lib_type R \
    --CPU 26 \
    --genome_guided_max_intron 10000 

echo "Finished processing of files"

date #end time of script

#trinity params
#max memory manages resources, edit if fails
#genome_guided uses rna seq that has been aligned
#output is for output directory
#SS_lib details strand direction
#cpu is for resource usage
#genome_guided_max this is to tell max intron, according to literature 10k is fine for plants