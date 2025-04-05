#!/bin/bash

date #print start time of script
genome_directory=$1 #input to genome directory
masked_repeats_directory=$2 #output directory

mkdir -p $masked_repeats_directory #makes output directory just in case

#check to ensure all args are used
if [ -z $genome_directory ] || [ -z $masked_repeats_directory ]; then 
    echo "Usage: ./repeatmasker.sh <genome_directory> <output_directory>"
    exit 1
fi

for genome in $genome_directory/*.fna; do #loop through directory for genome
    #check to ensure genome is in directory
    if [ ! -f $genome ];then
        echo "No genome in $genome_directory"
        exit 1
    fi

    #make a database name
    database="casuarina"

    #make a working directory
    mkdir -p $masked_repeats_directory/$database
    cd $masked_repeats_directory/$database || exit 1

    echo "Processing genome: $genome"

    #Build repeat modeler database
    BuildDatabase -name $database $genome
    
    #run repeatmodeler (output is consensi.fa.classified)
    if [ ! -f "consensi.fa.classified" ]; then #a check to make sure output does not exist already
        RepeatModeler -database $database -pa 20 -LTRStruct >> mask.log
    else
        echo "Repeat library already exist, skipping RepeatModeler"
    fi

    #run repeatmasker with repeat library
    RepeatMasker -pa 20 -gff -xsmall -lib consensi.fa.classified -dir . $genome

    #extract masked fasta for braker2
    masked_file=$(basename $genome).masked
    cp $masked_file $masked_repeats_directory/${database}.softmasked

    echo "Soft-masked genome saved to $masked_repeats_directory/${database}.softmasked.fasta"

    cd - >/dev/null #go back to previous directory
done

date #lets you know the end time of script