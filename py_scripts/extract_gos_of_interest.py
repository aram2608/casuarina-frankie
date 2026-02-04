import pandas as pd
from typing import Dict, List
import json
from argparse import ArgumentParser, ArgumentTypeError
from pathlib import Path

gos_df = pd.read_csv("data/gos_of_interest.tsv", sep="\t")
query = pd.read_csv("data/gos_to_annotate.tsv", sep="\t")

go_series = gos_df["GO.ID"]

go_series = go_series.to_list()

terms: Dict[str, List[str]] = {}

it = 0
for gene in query["query"]:
    matched_term = []
    row = query.iloc[it, 1]
    row = row.split(",")
    for go in row:
        if go in go_series:
            matched_term.append(go)
    it += 1
    terms[gene] = matched_term

out = json.dumps(terms, indent=4)

with open("gos_of_interest.json", "w") as file:
    file.write(out)


counts: Dict[str, List[str]] = {}
for k, v in terms.items():
    for c in v:
        if c not in counts:
            counts[c] = 1
        else:
            counts[c] += 1

counts = dict(sorted(counts.items(), key=lambda item: item[1], reverse=True))

out = json.dumps(counts, indent=4)

with open("go_counts.json", "w") as file:
    file.write(out)

def validate_file(filepath: str):
    path: Path = Path(filepath)
    if not path.exists():
        raise ArgumentTypeError(f"File {filepath} does not exist")
    return path

if __name__ == "__main__":
    parser:ArgumentParser = ArgumentParser(description="Extracts go terms of interest from a tsv.")
    parser.add_argument()
    # TODO: Create a cli interface to take files