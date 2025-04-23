import pandas as pd
import argparse

def main(blast_file, gene_file, blast_index, gene_index, output_file):
    """
    A script for finding matching genes between a blast output and differentially
    expressed genes from and RNA-seq study. 

    The goal is to add functional annotations/descriptions to the diff expressed genes.

    This script is ideal/meant for organisms that are non models and more difficult to work with.

    For proper usage make sure you know which indexes to use for gene matching.

    Usage:
    python3 geneMatcher.py input_file input_file2 index_for_file1 index_for_file2
    """

    # read in files for matching
    df = pd.read_csv(blast_file, sep="\t", header=None, index_col=False)
    df2 = pd.read_csv(gene_file, sep="\t", header=None, index_col=False)
    gene_column_1 = int(blast_index) # sets gene_column for blast file with
    gene_column_2 = int(gene_index) # sets gene_column for gene file
    col1 = df.columns[gene_column_1]
    col2 = df2.columns[gene_column_2]

    # finds overlapping data between both dataframe in the specified indexes
    overlap_data = pd.merge(df, df2[[col2]], left_on=col1, right_on=col2)

    # drop dupes
    overlap_data = overlap_data.drop_duplicates()

    # write to new text file
    overlap_data.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adds blast terms to a file with diff expressed genes.")
    parser.add_argument("blast_file", help="Input file with blast terms")
    parser.add_argument("gene_file", help="Input file with diff expressed genes")
    parser.add_argument("output_file", help="Output file with blast terms matched to diff expressed genes")
    parser.add_argument("--ID1", dest="blast_index", required=True, help="Index containing gene IDs")
    parser.add_argument("--ID2", dest="gene_index", required=True, help="Index containing gene IDs")
    args = parser.parse_args()

    main(args.blast_file, args.gene_file, args.blast_index, args.gene_index, args.output_file)