#!/bin/bash

input_GTF=$1 # input file to GTF file ---- /full/path/to/file
output_file=$(basename $input_GTF .gtf)_busco_fabales # output file name extraction from basename

# safety check to make sure args are met
if [ ! -f "$input_GTF" ]; then
    echo "Usage: <input_GTF>"
    exit
fi

# busco run parameters
busco -i braker_proteins.faa \
    -l fabales_odb10 \
    -o "$output_file" \
    -m protein

echo "Finished checking quality of $input_GTF"

# important params to know for working with casuarina
# -l fabales_odb10 # most biologically relevant, legumes like medicago
# -l eudicots_odb10 # eudicots a lot braoder and more general
# -l embryophyta_odb10 # all land plants, a lot more general