#!/bin/bash

# this line is for installing all eggnog databases
# the default is the current directory so comment out
# if you already have the databases downloaded
# download_eggnog_data.py --data_dir .

date # print script start time
protein_seqs=$1
annotations=$2

# safety check for arguments
if [ -z "$protein_seqs" ] || [ -z $annotations ]; then
    echo "Usage: <protein_sequences> <ouput_file_name>"
    exit 1
fi

# base eggnogmapper params
emapper.py -i $protein_seqs \
    -o $annotations \
    --cpu 10 \
    --data_dir . \
    --tax_scope 33090 \
    --target_orthologs all \
    --evalue 1e-3 \
    --override

echo "Finished annotating $protein_seqs"
date # script end time

# -i is the default input flag, takes proteins
# -o is the output flag
# --itype if using another type of input specifiy
# -m using mmseq search, needs prior setup before use
# --cpu is the # of threads
# --data_dir sets the path to the eggnog database
# --list_taxa can be used to see which taxa are in the databse
# --tax_scope sets a boundary for which taxa annotations can come from
# --target_orthologs all, adding this too, should increase our hits hopefully
# --evalue setting a cutoff as well to get only relevant matches

# yup gotta apply taxa scope, pretty sure my plant is not sniffing things or developing cancer
# 33090 is for viridaeplanta so all green plants