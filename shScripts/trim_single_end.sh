#!/bin/bash

#a script using trimmmomatic to trim single-end illumina reads
#loaded from the RON server

#run the following command to download the adapter file for Single End Illumina
#wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa

#whatever directory you ran the command in will be the input for the adapter
#directory

date #prints start time of script
input_directory=$1 #provide input directory
output_directory=$2 #provide output directory
adapters=$3 #provide path to adapter file

#added safety test to ensure all directories are provided
if [ -z $input_directory ] || [ -z $output_directory ] || [ -z $adapters ]; then
    echo "Usage: $0 <input_directory> <output_directory> </path/to/adapters>"
    exit 1
fi

echo $input_directory will now be processed
mkdir -p $output_directory #makes output directory just in case

#loop to trim files
for fastq_file in $input_directory/*;do #loop through files of input directory
    if [[ $fastq_file == *fastq.gz ]]; then #tests for fastq file
        file_name=$(basename $fastq_file) #creates basename for files
        output_file=$output_directory/trimmed_$file_name

        echo now trimming $fastq_file -> $output_file; #prints file to be trimmed

        trimmomatic SE -phred33 $fastq_file $output_file \
            ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
            
        #an error handling step in case a file is not processed for some reason
        if  [ $? -ne 0 ]; then # $? is a special variable that stores the status of the previous command, -ne 0 means a non zero error was thrown
            echo "Error Trimming $fastq_file"
        fi
    fi
done
date #logging of script run time

#following parameters applied
#TrueSeq3-SE adapters removed
#2 mismatches allowed
#30 palindrome clip threshold
#10 simple clip threshold
#leading:3 for the removal of low-qual base pairs from the beginning (below 3)
#trailing:3 for removing low-qual bases from the end (below 3)
#slidingwindow:4:15, scans with a 4 base window and cuts after drop of 15 in qual
#minlen:36 discards reads shorter than 36 bases