import sys
import pandas as pd
import numpy as np

def main():
    args = sys.argv[1:]
    if (len(args) != 4):
        raise Exception(
            "Please pass in prsm.tsv files")
    


    r1 = pd.read_csv(args[0], sep="\t")
    r2 = pd.read_csv(args[1], sep="\t")
    r1New = pd.read_csv(args[2], sep="\t")
    r2New = pd.read_csv(args[3], sep="\t")
    
    combined = pd.concat([r1, r2], ignore_index=True)
    combined.drop_duplicates(subset=["Scan(s)", "Proteoform"], inplace=True)
    combined.sort_values(by=["Scan(s)"], inplace=True)
    combined.to_csv("combined.csv", sep="\t")

    combinedNew = pd.concat([r1New, r2New], ignore_index=True)
    combinedNew.drop_duplicates(subset=["Scan(s)", "Proteoform"], inplace=True)
    combinedNew.sort_values(by=["Scan(s)"], inplace=True)
    combinedNew.to_csv("combinedNew.csv", sep="\t")
        


if __name__ == "__main__":
    main()