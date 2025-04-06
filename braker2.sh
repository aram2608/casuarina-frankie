#!/bin/bash

date #print start time

bam_directory=$1 #directory for bam files
genome=$2 #path to genome
output_directory=$3 #output directory

#create a list of bam files to train all at once
bam_list=$(ls "$bam_directory"/*.bam | paste -sd, -)

#check to see if all bams included
echo BAM list: $bam_list

#check to ensure all args met
if [ -z $bam_directory ] || [ -z $genome ] || [ -z $output_directory ]; then
    echo "Usage: ./braker2.sh <bam_directory> <genome> <output_directory>"
    exit 1
fi

export GENEMARK_PATH="/home/users/ja1473/gmes_linux_64_4"
export GM_KEY="/home/users/ja1473/gmes_linux_64_4/gm_key"

#test to ensure genemark info is set
if [ ! -f $GM_KEY ]; then
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
    --gff3 \
    --UTR=on

echo Finished annotating Casuarina glauca

date #print the end time of script
