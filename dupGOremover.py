import argparse
import pandas as pd

def main(input_txt, output_txt,):

    # import file with go terms as a data frame
    df = pd.read_csv(input_txt, sep="\t", index_col=False)
    
    # creates new empty data frame
    new_frame = pd.DataFrame()

    # creates an empty list to store go terms

    for item in df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Removes duplicate GO terms for each gene ID.")
    parser.add_argument("input_txt", help="Input go term file (TSV format)")
    parser.add_argument("output_txt", help="Output file with duplicates removed")
    #parser.add_argument("--ID", dest="geneID", required=True, help="Column name or index containing gene IDs (e.g., 'query_name' or 0)")
    args = parser.parse_args()

    main(args.input_txt, args.output_txt)#, args.geneID)