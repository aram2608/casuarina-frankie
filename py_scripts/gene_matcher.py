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
    blast_df = pd.read_csv(blast_file, sep="\t", header=None, index_col=False, low_memory=False)
    gene_df = pd.read_csv(gene_file, sep="\t", header=None, index_col=False)
    gene_column = int(gene_index)
    blast_column = int(blast_index)

    # extract and cleanup upregulated genes
    upregulated_genes = set(gene_df[gene_column].astype(str).str.strip())

    # find matching genes and save to new frame
    matched = blast_df[blast_df[blast_column].astype(str).str.strip().isin(upregulated_genes)]

    # write to new text file
    matched.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adds blast terms to a file with diff expressed genes.")
    parser.add_argument("blast_file", help="Input file with blast terms")
    parser.add_argument("gene_file", help="Input file with diff expressed genes")
    parser.add_argument("output_file", help="Output file with blast terms matched to diff expressed genes")
    parser.add_argument("--ID1", dest="blast_index", required=True, help="Index containing gene IDs")
    parser.add_argument("--ID2", dest="gene_index", required=True, help="Index containing gene IDs")
    args = parser.parse_args()

    main(args.blast_file, args.gene_file, args.blast_index, args.gene_index, args.output_file)
