import sys
import pandas as pd
from matplotlib_venn import venn2
import numpy as np
import csv
import matplotlib.pyplot as plt


def main():

    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in proteoform single files")
    
    r1 = pd.read_csv(args[0], sep="\t")
    r2 = pd.read_csv(args[1], sep="\t")

    r1["Origin"] = 1
    r2["Origin"] = 2

    combined = pd.concat([r1, r2], ignore_index=True)

    combined.sort_values(by="E-value")

    combined["1"] = combined["Origin"] == 1
    combined["2"] = combined["Origin"] == 2

    resultlist = pd.DataFrame(columns=combined.columns).astype(combined.dtypes)

    for i, row in combined.iterrows():
        is_duplicate = False
        for j, result in resultlist.iterrows():
            if row['Protein accession'] == result['Protein accession'] and abs(row['Precursor mass'] - result['Precursor mass']) < 1.2:
                if row["Origin"] == 1:
                    resultlist.loc[j, "1"] = True
                elif row["Origin"] == 2:
                    resultlist.loc[j, "2"] = True
                is_duplicate = True
                break
        
        if not is_duplicate:
            resultlist.loc[0 if pd.isnull(resultlist.index.max()) else resultlist.index.max() + 1] = row

    # Calculate the sizes for the Venn diagram
    only_set1 = resultlist[(resultlist["1"] == True) & (resultlist["2"] == False)].shape[0]
    only_set2 = resultlist[(resultlist["1"] == False) & (resultlist["2"] == True)].shape[0]
    intersection = resultlist[(resultlist["1"] == True) & (resultlist["2"] == True)].shape[0]

    # resultlist[(resultlist["1"] == True) & (resultlist["2"] == True)].to_csv("overlapR1R2.tsv", sep="\t")

    # Plot the Venn diagram
    venn2(subsets=(only_set1, only_set2, intersection), set_labels=('R1', 'R2'))
    plt.savefig(args[0].rsplit("/", maxsplit=1)[0] + "/proteoform.png", dpi=200)

if __name__ == "__main__":
    main()