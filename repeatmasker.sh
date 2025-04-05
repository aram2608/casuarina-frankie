#!/bin/bash

date #print start time of script
genome_directory=$1 #input to genome directory
masked_repeats_directory=$2 #output directory

mkdir -p $output_directory #makes output directroy just in case

#check to ensure all args are used
if [ -z $genome_directory ] || [ -z $output_directory ]; then 
    echo "Usage: ./repeatmasker.sh <genome_directory> <output_directory>"
    exit 1
fi

for genome in $genome_directory/*.fna; do #loop through directory for genome
    #check to ensure genome is in directory
    if [ -z $genome ];then
        echo "No genome in directory"
        exit 1
    fi
    #Chunk of script to mask repeats
    BuildDatabase -name casuarina $genome
    RepeatModeler -casuarina -pa 20 -LTRStruct >> mask.log
    RepeatMasker -pa 36 -gff -lib consensi.fa.classified -dir $masked_repeats_directory -xsmall $genome
done

date #lets you know the end time of script