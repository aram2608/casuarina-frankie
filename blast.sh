#!/bin/bash

sequences=$1 # full path to protein seqs
output=$2 # full path to output dir

# safety feature for arguments
if [ -f $sequences ] || [ -z $output ]; then
    echo "Usage: <path_to_proteins> <output_dir_path>"
    exit 1
fi

echo "Beginning blast search for $sequences"
# blastp params
blastp -query $sequences \
    -db nr \
    -outfmt '6 qseqid sseqid qlen slen qstart qend sstart send length pident positive mismatch gaps evalue bitscore qcovs' \
    -num_alignments 10 
    | > blastp_10hits.out

echo -e 'qseqid\tsseqid\tqlen\tslen\tqstart\tqend\tsstart\tsend\tlength\tpident\tpositive\tmismatch\tgaps\tevalue\tbitscore\tqcovs' > header1.tab
cat header1.tab blastp_10hits.out > blastp_10hits.tsv

echo "Finished processesing $sequences: new file = blastp_10hits.tsv"

# blastp params
# -query = input sequences
# -db is the database to be used
# -outfmt = produces tabular format
# -num_alignments = sets the number of alignments for hits