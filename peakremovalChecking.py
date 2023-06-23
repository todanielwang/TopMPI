import sys
import pandas as pd
import numpy as np

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in prsm.tsv files")
    


    r1 = pd.read_csv(args[0], sep="\t")
    r2 = pd.read_csv(args[1], sep=",")
    
    new = r2[r2.Proteoform.isin(r1.Proteoform) == False]
    new.sort_values("E-value", inplace=True)

    new.to_csv("stuff.csv")   


if __name__ == "__main__":
    main()