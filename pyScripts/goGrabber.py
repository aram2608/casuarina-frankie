import pandas as pd
import argparse

def main(input_tsv, output_file):

    """
    Extracts gene IDs and GO terms from an InterProScan TSV file and writes to a new TSV file.
    To prepare the InterProScan file for extraction, use the following command on your output...

    grep "GO" example_interproscan.tsv > my_new_file.tsv

    This will ensure you have only pulled genes/proteins that are guaranteed to have GO terms.
    
    Next ensure that the GO terms are in the 14th column, in pandas this is represented as index[13]

    Example Usage: python3 goGrabber.py /path/to/my_new_file.tsv /path/to/workdir/goterms.tsv
    """

    # load TSV as a dataframe
    df = pd.read_csv(input_tsv, sep="\t", header=None)

    # creates a new empty data frame
    new_data_frame = pd.DataFrame()

    # insert geneids and goterms into new data frame
    new_data_frame.insert(0, "gene_id", df[0]) # inserts into first index
    new_data_frame.insert(1, "go_terms", df[13]) # inserts into second index
    
    # writes to new csv file
    new_data_frame.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Gene IDs and GO terms from InterProScan tsv.")
    parser.add_argument("input_tsv", help="Input annotation file (TSV format)")
    parser.add_argument("output_file", help="Name of output file")
    args = parser.parse_args()

    main(args.input_tsv, args.output_file)
