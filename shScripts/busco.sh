#!/bin/bash

protein_dir=$1 # input to directory

# safety check to make sure args are met
if [ -z "$protein_dir" ]; then
    echo "Usage: <protein_directory>"
    exit
fi

# looping through protein directory
for prot_GTF in "$protein_dir"/*.faa; do
    # creating output directories
    base_name=$(basename "$prot_GTF" .faa)

    output_dir_1="${base_name}_busco_fabales"
    output_dir_2="${base_name}_busco_eudicots"
    output_dir_3="${base_name}_busco_embryophyta"

    # busco run parameters
    busco -i "$prot_GTF" -l fabales_odb10 -o "$output_dir_1" -m protein
    busco -i "$prot_GTF" -l eudicots_odb10 -o "$output_dir_2" -m protein
    busco -i "$prot_GTF" -l embryophyta_odb10 -o "$output_dir_3" -m protein

    echo "Finished checking quality of $prot_GTF"
done

# important params to know for working with casuarina
# -l fabales_odb10    # most biologically relevant, legumes like medicago
# -l eudicots_odb10   # eudicots a lot broader and more general
# -l embryophyta_odb10 # all land plants, a lot more general