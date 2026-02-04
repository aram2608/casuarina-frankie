import pandas as pd
from typing import List
from pathlib import Path
import sys
from argparse import ArgumentParser, ArgumentTypeError


def get_files(directory: Path, suffix: str) -> List[Path]:
    print("Searching for files...\n")
    data_dir = Path(directory)
    glob_pattern = "*" + suffix
    enriched_files = list(data_dir.glob(glob_pattern))
    print(f"Retrieved {len(enriched_files)} files from {directory}")
    return enriched_files


def process_files(enriched_files: List[Path], id: str) -> List[str]:
    go_terms = []
    for file in enriched_files:
        # TODO: This expects a tsv, allow csvs as well
        df = pd.read_csv(file, sep="\t")
        if id not in df.keys():
            print(f"The following column ID was not found in {file}.\nID: {id}")
            print("Aborting...")
            sys.exit(1)
        sub = df[id]
        for i in sub:
            go_terms.append(i)
    return go_terms

def main(directory: Path, id: str, suffix: str, output: str):
    ok_suffix: List[str] = ["csv", "tsv"]
    if suffix not in ok_suffix:
        print(f"Unexpected file type: {suffix}")
        return

    enriched_files: List[Path] = get_files(directory, suffix)
    go_terms: List[str] = process_files(enriched_files, id)

    go_terms: List[str] = list(set(go_terms))
    new_df: pd.DataFrame = pd.DataFrame(go_terms, columns=[id])
    print(f"Creating {output}")
    new_df.to_csv(output, sep="\t")


def validate_file(filepath: str):
    path: Path = Path(filepath)
    if not path.exists():
        raise ArgumentTypeError(f"File {filepath} does not exist")
    return path


if __name__ == "__main__":
    parser: ArgumentParser = ArgumentParser(
        description="Concat go terms from separate files."
    )
    parser.add_argument(
        "-p",
        "--path",
        dest="path",
        type=validate_file,
        required=True,
        help="Path to directory.",
    )
    parser.add_argument(
        "-ID",
        dest="ID",
        required=True,
        type=str,
        help="The column name containing the gos terms.",
    )
    parser.add_argument(
        "-t",
        "--type",
        dest="file_type",
        required=True,
        help="The type of file to be parsed. Either csv or tsv",
    )
    parser.add_argument(
        "-o", "--output", dest="output", required=True, help="The output file path."
    )

    args = parser.parse_args()

    main(args.path, args.ID, args.file_type, args.output)
