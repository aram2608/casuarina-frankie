#!/bin/bash

date #print start time

bam_directory=$1 #directory for bam files
genome=$2 #path to genome
gene_mark_key_dir=$3 #where the gene_mark key is
output_directory=$4 #output directory

#create a list of bam files to train all at once
bam_list=$(ls "$bam_directory"/*.bam | paste -sd, -)

#check to see if all bams included
echo BAM list: $bam_list

#check to ensure all args met
if [ -z $bam_directory ] || [ -z $genome ] || [ -z $gene_mark_key_dir ] || [ -z $output_directory ]; then
    echo "Usage: ./braker2.sh <bam_directory> <genome> <gene_mark_key_dir> <output_directory>"
    exit 1
fi

export GENEMARK_PATH=$gene_mark_key_dir
export GM_KEY="$gene_mark_key_dir/gm_key"

#test to ensure genemark info is set
if [ ! -f GENEMARK_PATH ] || [ ! -f $gm_key ]; then
    echo Missing GeneMark parameters
    exit 1
fi

#running braker2
braker.pl \
    --species=casuarina_glauca \
    --genome=$genome \
    --bam=$bam_list \
    --softmasking \
    --workingdir=$output_directory \
    --cores=26 \
    --gff3

echo Finished annotating $sample_name


date #print the end time of script
