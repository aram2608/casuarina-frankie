#!/bin/bash

date  # Print start time of script

genome_directory=$1               # Input to genome directory
masked_repeats_directory=$2       # Output directory

# Check to ensure all args are used
if [ -z "$genome_directory" ] || [ -z "$masked_repeats_directory" ]; then 
    echo "Usage: ./repeatmasker.sh <genome_directory> <output_directory>"
    exit 1
fi

mkdir -p "$masked_repeats_directory"  # Make output directory if needed

# Loop through genomes
for genome in "$genome_directory"/*.fna; do
    if [ ! -f "$genome" ]; then
        echo "No genome found in $genome_directory"
        exit 1
    fi

    # Set hardcoded database name
    database="casuarina"

    # Make a working directory
    mkdir -p "$masked_repeats_directory/$database"
    cd "$masked_repeats_directory/$database" || exit 1

    echo "Processing genome: $genome"

    # Build RepeatModeler database
    BuildDatabase -name "$database" "$genome"
    
    # Run RepeatModeler (without -LTRStruct)
    if [ ! -f "consensi.fa.classified" ]; then
        RepeatModeler -database "$database" -pa 20 >> mask.log
    else
        echo "Repeat library already exists, skipping RepeatModeler"
    fi

    # Check that RepeatModeler succeeded
    if [ ! -s "consensi.fa.classified" ]; then
        echo "ERROR: RepeatModeler failed or produced empty output."
        exit 1
    fi

    # Run RepeatMasker
    RepeatMasker -pa 20 -gff -xsmall -lib consensi.fa.classified -dir . "$genome"

    # Extract soft-masked FASTA for BRAKER2
    masked_file="$(basename "$genome").masked"
    if [ -f "$masked_file" ]; then
        cp "$masked_file" "$masked_repeats_directory/${database}.softmasked.fasta"
        echo "Soft-masked genome saved to: $masked_repeats_directory/${database}.softmasked.fasta"
    else
        echo "ERROR: RepeatMasker did not produce the masked FASTA."
        exit 1
    fi

    cd - > /dev/null  # Return to previous directory
done

date  # Print end time of script
