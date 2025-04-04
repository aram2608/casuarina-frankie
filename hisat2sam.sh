#!/bin/bash

date #print start time
trimmed_fastq_directory=$1
input_index_hisat2=$2
output_SAM_directory=$3

#safety feature to ensure all inputs are included
if [ -z $trimmed_fastq_directory ] || [ -z $input_index_hisat2 ] || [ -z $output_SAM_directory ]; then
    echo "Usage: ./gene_model_prep.sh <trimmed_fastq_directory> <input_index_directory>"
    exit 1
fi

mkdir -p $output_SAM_directory #makes directory just in case

#using hisat2 to create SAM files
for trimmed_fastq in $input_fastq_directory/*.fastq.gz; do #loop directory and search for zipped fastq files
    echo Starting alignment for $trimmed_fastq
    base_name=$(basename $trimmed_fastq) #extracts basename from input
    output_SAM=$output_SAM_directory/$base_name_name.sam #creates SAM file name from fastq file name

    hisat2 --phred33 --dta -x $input_index_hisat2 -U $trimmed_fastq -S $output_sum

    #check if hisat2 command failed

    if [ $? -ne 0]; then
        echo "Failure for $trimmed_fastq alignment"
    fi
done

date #tells you how longs the program has run

#hisat2 parameters
#--phred33 quality score
#--dta reports in a form useful for downstream transcriptome analysis
#-x is for index directory
#-S is for output sam file