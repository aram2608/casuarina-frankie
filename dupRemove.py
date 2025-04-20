import argparse
import pandas as pd

def remove_duplicate_fasta(input_txt, output_txt, geneID):
    
    """
    Removes isoform duplicates based on gene IDs from a tab-delimited file.
    You can provide the gene column as a name or index (0-based).

    Example:
    python3 script.py input.tsv output.tsv --ID 0
    python3 script.py input.tsv output.tsv --ID gene_name
    """

    df = pd.read_csv(input_txt, sep='\t')
    original_count = len(df)

    # Determine if geneID is a digit (index) or column name
    try:
        geneID_index = int(geneID)
        gene_col = df.columns[geneID_index]
    except ValueError:
        gene_col = geneID
        if gene_col not in df.columns:
            raise ValueError(f"Column '{gene_col}' not found. Available columns: {list(df.columns)}")
    except IndexError:
        raise IndexError(f"Column index {geneID_index} is out of bounds. File has {len(df.columns)} columns.")

    # Create base ID by removing isoform suffix
    df['base_id'] = df[gene_col].apply(lambda x: x.rsplit('.', 1)[0])

    # Drop duplicates based on the base_id
    df_unique = df.drop_duplicates(subset='base_id', keep='first')

    # Replace the original ID column with the base_id
    df_unique[gene_col] = df_unique['base_id']
    df_unique = df_unique.drop(columns=['base_id'])

    new_count = len(df_unique)
    print(f"Original entries: {original_count}")
    print(f"Unique genes kept: {new_count}")
    print(f"Isoforms removed: {original_count - new_count}")

    df_unique.to_csv(output_txt, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove isoforms from an annotation file.")
    parser.add_argument("input_txt", help="Input annotation file (TSV format)")
    parser.add_argument("output_txt", help="Output file with duplicates removed")
    parser.add_argument("--ID", dest="geneID", required=True,
                        help="Column name or index containing gene IDs (e.g., 'query_name' or 0)")
    args = parser.parse_args()

    remove_duplicate_fasta(args.input_txt, args.output_txt, args.geneID)
