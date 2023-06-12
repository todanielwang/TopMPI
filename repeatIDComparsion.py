import sys
import pandas as pd
import numpy as np

def main():
    args = sys.argv[1:]
    if (len(args) != 12):
        raise Exception(
            "Please pass in prsm.tsv files")
    


    r1 = pd.read_csv(args[0], sep="\t")
    r2 = pd.read_csv(args[1], sep="\t")
    r3 = pd.read_csv(args[2], sep="\t")
    r4 = pd.read_csv(args[3], sep="\t")
    r5 = pd.read_csv(args[4], sep="\t")
    r6 = pd.read_csv(args[5], sep="\t")
    r7 = pd.read_csv(args[6], sep="\t")
    r8 = pd.read_csv(args[7], sep="\t")
    r9 = pd.read_csv(args[8], sep="\t")
    r10 = pd.read_csv(args[9], sep="\t")
    r11 = pd.read_csv(args[10], sep="\t")
    r12 = pd.read_csv(args[11], sep="\t")
    
    
    combined = pd.concat([r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12], ignore_index=True)
    combined.drop_duplicates(subset=["Protein accession"], inplace=True)
    combined.to_csv("combined.csv", sep="\t")
        


if __name__ == "__main__":
    main()