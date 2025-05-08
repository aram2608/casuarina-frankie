import argparse
import pandas as pd

def main(input_txt, output_txt):

    """
    This is a function to remove redundant KEGG map terms. The function assumes a csv
    file with quotes. The function assumes the first index contains the query ID and
    the second index contains the KEGG_Pathway.

    Example Usage:
    python3 dupKEGGremover.py input_csv output_csv
    """

    # import file with kegg terms as a data frame
    df = pd.read_csv(input_txt)
    df.columns = ['query', 'KEGG_Pathway']
    print(df.head())

    # initialize an empty dictionary to store KEGG_terms
    kegg_terms = {}

     # iterates over each row and stores as a key if a gene_id and value if a GO: term list
    for _, row in df.iterrows(): # iterates over rows, the _, tells python to ignore the index on purpose
        query = row['query']
        terms = row['KEGG_Pathway'].split(',')
        unique_KEGG = sorted(set(terms)) # removes duplicates
        kegg_terms[query] = unique_KEGG # set query as key and unique_KEGG as values

    # write output file
    with open(output_txt, "w") as out:
        for query, terms in kegg_terms.items():
            out.write(f"{query}\t{','.join(terms)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Removes duplicate GO terms for each gene ID.")
    parser.add_argument("input_txt", help="Input go term file (TSV format)")
    parser.add_argument("output_txt", help="Output file with duplicates removed")
    args = parser.parse_args()

    main(args.input_txt, args.output_txt)

# I don't know how to add column names so run the following in your terminal

# echo -e 'query\tKEGG_Pathway' >> header.txt
# cat header.txt output_file_name >> new_file_name