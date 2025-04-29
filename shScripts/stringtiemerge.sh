#!/bin/bash

stringtie_list=$1 # list of string tie assemblies to merge, full path

if [ ! -f $stringtie_list ]; then
    echo "Usage: <path_to_list_of_assemblies>"

# stringtie merge params
stringtie --merge \
    -o merged.gtf \
    -m 200 \
    -c 1.5 \
    -F 0.5 \
    -f 0.05 \
    -p 8 \
    $stringtie_list

# the following params were applied
# -m minimum input transcript length 200
# -c min coverage length 1.5
# -F min fpkm 0.5
# -f minmum isofrom fraction 0.05
# -p 8 threads