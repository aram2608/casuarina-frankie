import argparse
import pandas as pd

def remove_duplicate_fasta(input_txt, output_txt):
    df = pd.read_csv(input_txt, sep='\t')
    original_count = len(df)

    # Create base ID by removing isoform suffix
    df['base_id'] = df['qseqid'].apply(lambda x: x.rsplit('.', 1)[0])

    # Keep only one entry per base gene
    df_unique = df.drop_duplicates(subset='base_id', keep='first')

    # Replace the query field with the cleaned ID
    df_unique.loc[:, 'qseqid'] = df_unique['base_id']
    df_unique = df_unique.drop(columns=['base_id'])

    new_count = len(df_unique)
    print(f"Original entries: {original_count}")
    print(f"Unique genes kept: {new_count}")
    print(f"Isoforms removed: {original_count - new_count}")

    df_unique.to_csv(output_txt, sep='\t', index=False)

# main arguments to run the script
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove isoforms from an eggnogmapper annotation file.")
    parser.add_argument("input_txt", help="Input annotation file (TSV format)")
    parser.add_argument("output_txt", help="Output file with duplicates removed")
    args = parser.parse_args()

    remove_duplicate_fasta(args.input_txt, args.output_txt)
