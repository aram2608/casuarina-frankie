import argparse
import pandas as pd

def main(input_txt, output_txt,):

    """
    This is a function to keep unique goterms from an input TSV file. The script assumes
    there are no headers and that index [0] are the gene_ids and index[1] are the goterms.

    The purpose of this script is to remove redundancy in goterms for downstream GO
    enrichment analysis from custom GO annotations. 

    Ideally, the goterms were annotated using InterProScan where redundance is introduced
    from the multiple search tools use such as Panther, Pfam, InterPro, etc..

    Example Usage:
    python3 dupGOremover.py input_tsv output_tsv
    """

    # import file with go terms as a data frame
    df = pd.read_csv(input_txt, sep="\t", header=None, index_col=False)
    df.columns = ['gene_ids', 'go_terms']

    # creates an empty dictionary to store our gene_ids and goterms
    gene_GO = {}

    # iterates over each row and stores as a key if a gene_id and value if a GO: term list
    for id, row in df.iterrows(): # iterates over rows
        gene = row['gene_ids']
        terms = row['go_terms'].split(',')
        unique_gos = sorted(set(terms)) # removes duplicates
        gene_GO[gene] = unique_gos

    # write output file
    with open(output_txt, "w") as out:
        for gene, terms in gene_GO.items():
            out.write(f"{gene}\t{','.join(terms)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Removes duplicate GO terms for each gene ID.")
    parser.add_argument("input_txt", help="Input go term file (TSV format)")
    parser.add_argument("output_txt", help="Output file with duplicates removed")
    #parser.add_argument("--ID", dest="geneID", required=True, help="Column name or index containing gene IDs (e.g., 'query_name' or 0)")
    args = parser.parse_args()

    main(args.input_txt, args.output_txt)