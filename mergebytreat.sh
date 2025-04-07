#!/bin/bash

# Create output directory
mkdir -p merged_bams

# Define treatment groups and time points
declare -a treatments=("BR" "SNR" "R")
declare -a timepoints=("0" "24" "48")

# Loop over each treatment and timepoint
for treatment in "${treatments[@]}"; do
  for time in "${timepoints[@]}"; do
    # Create a subdirectory for each treatment group
    outdir="merged_bams/${treatment}_${time}"
    mkdir -p "$outdir"
    
    # Find relevant BAM files
    files=$(ls trimmed_${time}${treatment}*sorted.bam 2>/dev/null)
    
    # Check if files exist for the group
    if [ -z "$files" ]; then
      echo "No files found for ${treatment} at ${time}h"
      continue
    fi

    # Output filename
    merged_bam="${outdir}/${treatment}_${time}_merged.bam"
    sorted_bam="${outdir}/${treatment}_${time}_merged_sorted.bam"

    echo "Merging files for ${treatment} at ${time}h..."
    # Merge the BAM files
    samtools merge -f "$merged_bam" $files

    # Sort the merged BAM file
    echo "Sorting merged file..."
    samtools sort -o "$sorted_bam" "$merged_bam"

    # Index the sorted BAM
    echo "Indexing sorted file..."
    samtools index "$sorted_bam"

    # Remove the unsorted merged BAM
    rm "$merged_bam"

    echo "Finished processing ${treatment} at ${time}h"
    echo "Output: $sorted_bam"
  done
done
