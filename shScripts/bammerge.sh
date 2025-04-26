#!/bin/bash

date #print start time

bam_dir=$1 #input directory
output_dir=$2 #output directory
merged_name="merged_all.bam" #outputname

mkdir -p $output_dir

#safety check for args
if [ -z "$bam_dir" ] || [ -z "$output_dir" ]; then
    echo "Usage: ./bammerge.sh <bam_dir> <out_put>"
    exit 1
fi

#merge all bam files in directory
samtools merge -@ 8 "$output_dir/$merged_name" "$bam_dir"/*.bam

#index new file
samtools index "$output_dir/$merged_name"

#sucess message
echo "Merged BAM created at: $output_dir/$merged_name"
echo "Index created at: $output_dir/${merged_name}.bai"

date #script end time