#!/bin/bash

#to display the time at the start of the script
date
input_directory=$1 #provide path to input directory
output_directory=$2 #provide path to output_directory

#an added safety test to make sure both inputs provided
if [ -z "$input_directory" ] || [ -z "$output_directory" ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

echo "$input_directory will now be processed" #echos the directory to be processed
mkdir -p $output_directory #makes output directory just incase

#FastQC analysis for the files in the provided directory
for fastq_file in $input_directory/*;do #loops through directory
    echo Starting $fastq_file
    if [[ $fastq_file == *.fastq.gz ]];then #test if files are fastq files, must be zipped files
        fastqc -o $output_directory $fastq_file #performs fastQC analysis for all files in directory
    else
    echo $fastq_file is not a fastq file #prints which files can not be processed
    fi
done
#lets you now how long the script has run
date  