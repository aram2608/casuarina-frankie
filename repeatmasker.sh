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
    if [ -! -f $genome ];then
        echo "No genome in $genome_directory"
        exit 1
    fi

    #make a database name
    database="casuarina"

    #make a working directory
    mkdir -p $masked_repeats_directory/$base_file
    cd $masked_repeats_directory/$base_file || exit 1

    #Build repeat modeler database
    BuildDatabase -name $database $genome
    
    #run repeatmodeler (output is consensi.fa.classified)
    RepeatModeler -database $database -pa 20 -LTRStruct >> mask.log

    #run repeatmasker with repeat library
    RepeatMasker -pa 20 -gff -xsmall -lib consensi.fa.classified -dir . $genome

    cd - >/dev/null #go back to previous directory
done

date #lets you know the end time of script