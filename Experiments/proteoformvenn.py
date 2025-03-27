import sys
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import util

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in proteoform single files")
    
    r1 = util.read_tsv(args[0])
    r2 = util.read_tsv(args[1])

    # r1 = r1[~(r1["Protein accession"].str.contains("DECOY"))]
    # r1 = r1[r1["E-value"] < 0.01]
    # print(r1.shape[0])

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
    only_set1 = resultlist[(resultlist["1"] == True) & (resultlist["2"] == False)]
    only_set2 = resultlist[(resultlist["1"] == False) & (resultlist["2"] == True)]
    intersection = resultlist[(resultlist["1"] == True) & (resultlist["2"] == True)]

    only_set1.to_csv("Lost.tsv", sep="\t", index=False)
    only_set2.to_csv("New.tsv", sep="\t", index=False)

    # interestset = pd.concat([only_set1, intersection], ignore_index=True)

    # interestset.to_csv("largeoverlap.tsv", sep="\t", index=False)

    # Calculate the sizes for the Venn diagram
    only_set1 = resultlist[(resultlist["1"] == True) & (resultlist["2"] == False)].shape[0]
    only_set2 = resultlist[(resultlist["1"] == False) & (resultlist["2"] == True)].shape[0]
    intersection = resultlist[(resultlist["1"] == True) & (resultlist["2"] == True)].shape[0]

    # resultlist[(resultlist["1"] == False) & (resultlist["2"] == True)].to_csv("newproteoforms.tsv", sep="\t")

    plt.rcParams.update({'font.size': 18})

    venn = venn2(subsets=(only_set1, only_set2, intersection), set_labels=('No Charge Condition', 'With Charge Condition'))

    # # Get the set labels (titles) and manually align them
    # set_labels = venn.set_labels

    # if set_labels[0]:  # Set A label
    #     set_labels[0].set_position((-0.3, -0.6))  # Adjust manually for alignment

    # if set_labels[1]:  # Set B label
    #     set_labels[1].set_position((0.3, -0.6)) 

    plt.tight_layout()
    plt.savefig("./venndiagram.svg", format='svg', dpi=800)

if __name__ == "__main__":
    main()