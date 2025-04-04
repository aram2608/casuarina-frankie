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

mkdir -p $output_BAM_directory #makes directory just in case

#using hisat2 to create SAM files
for trimmed_fastq in $input_fastq_directory/*; do #loop directory
    echo Starting alignment for $trimmed_fastq
    if [[ $trimmed_fastq == *.fastq.gz ]]; then #tests if file is a fastq
        file_name=$(basename $trimmed_fastq) #extracts basename from input
        output_SAM=$output_BAM_directory/$file_name.sam #creates

        for index in $input_index_hisat2/*; do #loop through index directory
        
            hisat2 --phred33 --dta -x $index -U $trimmed_fastq -S $output_SAM

            #safety check to ensure previous command was successful
            if [ $? -ne 0 ]; then
                echo Failure in processing $trimmed_fastq
            fi
        done
    fi
done

date #tells you how longs the program has run

#hisat2 parameters
#--phred33 quality score
#--dta reports in a form useful for downstream transcriptome analysis
#-x is for index directory
#-S is for output sam file