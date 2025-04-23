import pandas as pd
import argparse

def main(input_file, input_file2, index_1, index_2, output_txt):
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
    df = pd.read_csv(input_file, sep="\t", header=None, index_col=False)
    df2 = pd.read_csv(input_file2, sep="\t", header=None, index_col=False)
    gene_column_1 = int(index_1)
    gene_column_2 = int(index_2)

    # finds overlapping data between both dataframe in the specified indexes
    overlap_data = pd.merge(df[[gene_column_1]], df2[[gene_column_2]]).drop_duplicates()

    # write to new text file
    overlap_data.to_csv(output_txt, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adds blast terms to a file with diff expressed genes.")
    parser.add_argument("input_txt", help="Input file with blast terms")
    parser.add_argument("input_txt2", help="Input file with diff expressed genes")
    parser.add_argument("output_txt", help="Output file with blast terms matched to diff expressed genes")
    parser.add_argument("--ID1", dest="index_1", required=True, help="Index containing gene IDs")
    parser.add_argument("--ID2", dest="index_2", required=True, help="Index containing gene IDs")
    args = parser.parse_args()

    main(args.input_txt,args.input_txt2 ,args.output_txt)