import sys
import pandas as pd
from matplotlib_venn import venn3
import numpy as np
import csv
import matplotlib.pyplot as plt


def main():

    args = sys.argv[1:]
    if (len(args) != 3):
        raise Exception(
            "Please pass in proteoform single files")
    
    r1 = pd.read_csv(args[0], sep="\t", skiprows=29)
    r2 = pd.read_csv(args[1], sep="\t", skiprows=29)
    r3 = pd.read_csv(args[2], sep="\t", skiprows=29)

    r1 = r1[~r1['Protein accession'].str.contains('DECOY')]
    r2 = r2[~r2['Protein accession'].str.contains('DECOY')]
    r3 = r3[~r3['Protein accession'].str.contains('DECOY')]

    r1["Origin"] = 1
    r2["Origin"] = 2
    r3["Origin"] = 3

    combined = pd.concat([r1, r2, r3], ignore_index=True)

    combined.sort_values(by="E-value")

    combined["1"] = combined["Origin"] == 1
    combined["2"] = combined["Origin"] == 2
    combined["3"] = combined["Origin"] == 3

    resultlist = pd.DataFrame(columns=combined.columns).astype(combined.dtypes)

    for i, row in combined.iterrows():
        is_duplicate = False
        for j, result in resultlist.iterrows():
            if row['Protein accession'] == result['Protein accession'] and abs(row['Precursor mass'] - result['Precursor mass']) < 1.2:
                if row["Origin"] == 1:
                    resultlist.loc[j, "1"] = True
                elif row["Origin"] == 2:
                    resultlist.loc[j, "2"] = True
                elif row["Origin"] == 3:
                    resultlist.loc[j, "3"] = True
                is_duplicate = True
                break
        
        if not is_duplicate:
            resultlist.loc[0 if pd.isnull(resultlist.index.max()) else resultlist.index.max() + 1] = row

    # Calculate the sizes for the Venn diagram
    only_set1 = resultlist[(resultlist["1"] == True) & (resultlist["2"] == False) & (resultlist["3"] == False)].shape[0]
    only_set2 = resultlist[(resultlist["1"] == False) & (resultlist["2"] == True) & (resultlist["3"] == False)].shape[0]
    only_set3 = resultlist[(resultlist["1"] == False) & (resultlist["2"] == False) & (resultlist["3"] == True)].shape[0]
    intersection12 = resultlist[(resultlist["1"] == True) & (resultlist["2"] == True) & (resultlist["3"] == False)].shape[0]
    intersection13 = resultlist[(resultlist["1"] == True) & (resultlist["2"] == False) & (resultlist["3"] == True)].shape[0]
    intersection23 = resultlist[(resultlist["1"] == False) & (resultlist["2"] == True) & (resultlist["3"] == True)].shape[0]
    intersection123 = resultlist[(resultlist["1"] == True) & (resultlist["2"] == True) & (resultlist["3"] == True)].shape[0]

    # Plot the Venn diagram
    venn3(subsets=(only_set1, only_set2, intersection12, only_set3, intersection13, intersection23, intersection123),
        set_labels=('Replicate 1', 'Replicate 2', 'Replicate 3'))

    plt.savefig("proteoformoverlap.png", dpi=200)

    plt.show()

if __name__ == "__main__":
    main()