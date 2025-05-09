#!/bin/bash

sequences=$1 # full path to protein seqs
output=$2 # full path to output dir

# safety feature for arguments
if [ ! -f "$sequences" ] || [ -z "$output" ]; then
    echo "Usage: <path_to_proteins> <output_dir_path>"
    exit 1
fi

# make output just in case
mkdir -p "$output"

echo "Beginning blast search for '$sequences'..."

# blastp params
blastp -query $sequences \
    -db nr \
    -outfmt '6 qseqid sseqid qlen slen qstart qend sstart send length pident positive mismatch gaps evalue bitscore qcovs' \
    -num_alignments 10 \
    -num_threads 12 \
    -out "$output/blastp_10hits.out"

echo -e 'qseqid\tsseqid\tqlen\tslen\tqstart\tqend\tsstart\tsend\tlength\tpident\tpositive\tmismatch\tgaps\tevalue\tbitscore\tqcovs' > "$output/header1.tab"
cat "$output/header1.tab" "$output/blastp_10hits.out" > "$output/blastp_10hits.tsv"

echo "Finished processesing $sequences: new file = blastp_10hits.tsv"

# blastp params
# -query = input sequences
# -db is the database to be used
# -outfmt = produces tabular format
# -num_alignments = sets the number of alignments for hits

# use the following command to download a database
# update_blastdb.pl --decompress nr [*]
# in this case it will be the non-redundant database

# ron path
# /home/share/databases/ncbi_nr/nr