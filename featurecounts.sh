#!/bin/bash

input_bams_file=$1 # file with names of bam files
annotation_file=$2 # file with annotations in GTF format
output_dir=$3 # output directory for count matrix

# safety check for input file
if [ ! -f "$input_bams_file" ] || [ ! -f "$annotation_file" ]; then
    echo "Usage: <input_bam_txt_file>"
    exit 1
fi

# makes the output dir just in case
mkdir -p $output_dir

# creates a base name for the output file
base_name=$(basename $input_bam_txt_file .txt)

# params for featurecounts run
featurecounts -s 2 \
    exon \
    -g gene_id \
    -a "$annotation_file" \
    -o "$output_dir"/${base_name}_counts.txt \
    $input_bam_txt_file

echo "Finished calculating counts"

# params for feature counts
# exon counts only coding sequences
# -s tells you strandness of rna, 2 is reverse
# -a is for annotation files, GTF
# -o is for output
# the input bam file is postional

# !!! Make sure this script is run in the same directory as the bam files