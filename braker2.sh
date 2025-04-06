#!/bin/bash

date #print start time

bam_directory=$1 #directory for bam files
genome=$2 #path to genome
output_directory=$3 #output directory

#check to ensure all args met
if [ -z $bam_directory ] || [ -z $genome ] || [ -z $output_directory ]; then
    echo "Usage: ./braker2.sh <bam_directory> <genome> <output_directory>"
    exit 1
fi

#create a list of bam files separated by commas
bam_list=$(ls "$bam_directory"/*.bam 2>/dev/null | paste -sd, -) #redirects file descriptor 2 to the unix trash can /dev/null

#check for empty bam list
if [ -z $bam_list ]; then
    echo No BAM files found in $bam_directory
    exit 1
fi

#check to see if all bams included
echo BAM list: $bam_list

#set genemark variable
export GENEMARK_PATH="/home/users/ja1473/gmes_linux_64_4"
export GM_KEY="$GENEMARK_PATH/gm_key"

#test to ensure genemark info is set
if [ ! -f $GM_KEY ]; then
    echo Missing GeneMark parameters
    exit 1
fi

#run braker2
braker.pl \
    --species=casuarina_glauca_v2 \
    --genome=$genome \
    --bam=$bam_list \
    --workingdir=$output_directory \
    --cores=26 \
    --UTR=on

echo Finished annotating Casuarina glauca

date #print the end time of script

#run braker2
#braker.pl \
    #--species=casuarina_glauca \
    #--genome=$genome \
    #--bam=$bam_list \
    #--softmasking \
    #--workingdir=$output_directory \
    #--cores=26 \
    #--UTR=on

#i tried these parameters last time and it failed, gonna retry with diff ones