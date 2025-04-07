#!/bin/bash

date #print start time

pasa_alignments=$1 #full path to pasa_gff3
pasa_orfs=$2 #full path to pasa_orfs
braker_hints=$3 #full path to braker hints
genome=$4 #full path to clean!! genome

#argument checks
if [ ! -f $pasa_alignments ] || [ ! -f $pasa_orfs ] || [ ! -f $braker_hints ] || [ ! -f $genome ]; then
    echo "Usage: ./evidencemodeler.sh <pasa_alignments> <pasa_orfs> <braker_hints> <genome>"
    exit 1
fi

#params for evidence modeler
EVidenceModeler \
    --sample_id casuarina_glauca \
    --genome $genome \
    --gene_predictions $braker_hints \
    --transcript_alignments $pasa_alignments \
    --OTHER_EVIDENCE $pasa_orfs \
    --segmentSize 100000 \
    --overlapSize 10000

echo "Finished creating gene model"
date #print end time

#all files should be in gff3 format