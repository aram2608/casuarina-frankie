#!/bin/bash

#a script using trimmmomatic to trim single-end illumina reads
#loaded from the RON server

source /home/ja1473/anaconda3/etc/profile.d/conda.sh
conda activate genomics

date #prints start time of script
input_directory=$1 
output_diretory=$2
adapters=

echo $input_directory will now be processed

#loop to trim files
for fastq_file in $input_directory/*;do #loop through files of input directory
    output_prefix=${trimmed_} #in order to attach trimmed prefix
    output_file=${output_prefix}${fastq_file} #new output file
    
    if [[ $fastq_file == *fastq.gz ]] #tests for fastq file
        then
        echo now trimming $fastq_file -> $output_file; #prints file to be trimmed

        trimmomatic SE -phred33 $fastq_file $output_file \
            ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    fi
done
date 

#following parameters applied
#TrueSeq3-SE adapters removed
#2 mismatches allowed
#30 palindrome clip threshold
#10 simple clip threshold
#leading:3 for the removal of low-qual base pairs from the beginning (below 3)
#trailing:3 for removing low-qual bases from the end (below 3)
#slidingwindow:4:15, scans with a 4 base window and cuts after drop of 15 in qual
#minlen:36 discards reads shorter than 36 bases