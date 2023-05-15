import sys
import pandas as pd
import numpy as np

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in prsm.tsv files")
    


    base = pd.read_csv(args[0], sep="\t")
    #peaksunRemovefd = pd.read_csv(args[1], sep="\t")
    peaksremoved = pd.read_csv(args[1], sep="\t")


    #result1 = pd.concat([base, peaksunRemovefd], ignore_index=True)
    result2 = pd.concat([base, peaksremoved], ignore_index=True)

    #result1["Scan(s)"].replace('', np.nan, inplace=True)
    #result1["Scan(s)"].fillna(method="pad", inplace=True)

    result2["Scan(s)"].replace('', np.nan, inplace=True)
    result2["Scan(s)"].fillna(method="pad", inplace=True)


    #result1 = result1[result1.duplicated(subset=["Scan(s)", "Protein accession"], keep=False)]
    result2 = result2[result2.duplicated(subset=["Scan(s)", "Protein accession"], keep=False)]

    #result1 = result1.sort_values(by=["Scan(s)"])
    result2 = result2.sort_values(by=["Scan(s)"])

    #result1.to_csv("result1.csv", sep="\t")
    result2.to_csv("result2.csv", sep="\t")
        


if __name__ == "__main__":
    main()