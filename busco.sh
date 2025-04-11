#!/bin/bash

input_GTF_protein_dir=$1 # input file to GTF file ---- /full/path/to/file
output_dir_1=$(basename "$input_GTF_protein_dir" .gtf)_busco_fabales # output dir
output_dir_2=$(basename "$input_GTF_protein_dir" .gtf)_busco_eudicots # output dir

# safety check to make sure args are met
if [ -z "$input_GTF_protein_dir" ]; then
    echo "Usage: <input_GTF>"
    exit
fi

for prot_GTF in $input_GTF_protein_dir/*.faa; do
    # busco run parameters
    busco -i "$prot_GTF" \
        -l fabales_odb10 \
        -o "$output_dir_1" \
        -m protein

    busco -i "$prot_GTF" \
        -l eudicots_odb10 \
        -o "$output_dir_2" \
        -m protein

echo "Finished checking quality of $prot_GTF"

# important params to know for working with casuarina
# -l fabales_odb10 # most biologically relevant, legumes like medicago
# -l eudicots_odb10 # eudicots a lot braoder and more general
# -l embryophyta_odb10 # all land plants, a lot more general