#!/bin/bash

date #print start time

bam_directory=$1 #directory for bam files
genome=$2 #path to genome
gene_mark_key_dir=$3 #where the gene_mark key is
output_directory=$4 #output directory

#check to ensure all args met
if [ -z $bam_directory ] || [ -z $genome ] || [ -z $gene_mark_key_dir ] || [ -z $output_directory ]; then
    echo "Usage: ./braker2.sh <bam_directory> <gene_mark_key_dir>"
    exit 1
fi

#look through bam files
for bam_file in $bam_directory/*.bam; do 
    braker.pl --species=casuarina_glauca \
    --genome=$genome \
    --bam=$bam_file \
    --softmasking \
    --workingdir=$output_directory
    echo Finished annotating $bam_file
done

date #print the end time of script
