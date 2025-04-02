#!/bin/bash

#to display the time at the start of the script
date
directory=$1 #stores a directory from a command line argument

echo "$directory will now be processed" #echos the directory to be processed

#FastQC analysis for the files in the provided directory
for fastq_file in $directory/*;do #loops through directory
    echo Starting $fastq_file
    if [[ $fastq_file == *fastq.gz ]];then #test if files are fastq files, must be zipped files
        fastqc -o . $fastq_file #perfomrs fastQC analysis for all files in directory
    else
    echo $fastq_file is not a fastq file #prints which files can not be processed
    fi
done
#lets you now how long the script has run
date 
