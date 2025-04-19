import pandas as pd
import argparse

# a function to fetch GO terms from an interproscan run
def main(input_tsv, output_file):

    # load TSV as a dataframe
    df = pd.read_csv(input_tsv, sep="\t", header=None)

    # creates a new empty data frame
    new_data_frame = pd.DataFrame()

    # insert geneids and goterms into new data frame
    new_data_frame.insert(0, "gene_id", df[0])
    new_data_frame.insert(1, "go_terms", df[13])
    
    # writes to new csv file
    new_data_frame.to_csv(output_file, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove isoforms from an eggnogmapper annotation file.")
    parser.add_argument("input_tsv", help="Input annotation file (TSV format)")
    parser.add_argument("output_file", help="Name of output file")
    args = parser.parse_args()

    main(args.input_tsv, args.output_file)
