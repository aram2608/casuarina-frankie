#!/bin/bash

date #print script start time

trinity_file=$1 #path to trinity file
$genome=$2 #path to genome, use cleaned verison!

#arg check to ensure files exist
if [ ! -f "$trinity_file" ] || [ ! -f "$genome" ]; then
    echo "Usage: ./pasa.sh <trinity_file> <genome>"
    exit 1
fi

#run pasa
Launch_PASA_pipeline.pl \
    -c alignAssembly.config \
    -C \
    -R \
    -g "$genome" \
    --ALIGNERS blat,gmap\
    -t "$trinity_file" \
    --transcribed_is_aligned_orient

echo "Finished procesing '$trinity_file'"
date #script end time
