import pandas as pd

go_df = pd.read_csv("go_terms.tsv", sep="\t")

gen = (go_list.split(",") for go_list in go_df["GOs"])

for go in next(gen):
    print(go)