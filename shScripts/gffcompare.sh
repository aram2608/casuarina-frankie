#!/bin/bash

stringtie_output=$1 # input to stringtie output, use full path /path/to/file.gtf
braker_annotation=$2 # input to braker annotation, use full path /path/to/file.gtf

# safety check for arguments

if [ ! -f "$stringtie_output" ] || [ ! -f "$braker_annotation" ]; then
    echo "Usage: input_to_stringtie_gtf input_to_annotation>"
    exit 1
fi

# gffcompare params

gffcompare -r "$braker_annotation" \
    -o combined \
    "$stringtie_output"

# this script merges the stringtie assemblies with the braker run to improve annotation