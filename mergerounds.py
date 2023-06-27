import sys
import pandas as pd
import numpy as np

def main():
    args = sys.argv[1:]
    


    r1 = pd.read_csv(args[0], sep="\t")
    r2 = pd.read_csv(args[1], sep="\t")
    
    
    combined = pd.concat([r1, r2], ignore_index=True)
    Proteoform = combined.drop_duplicates(subset=["Proteoform"])
    Protein = combined.drop_duplicates(subset=["Protein accession"])
    Proteoform.to_csv("Proteoform.csv", sep="\t")
    Protein.to_csv("Protein.csv", sep="\t")
        


if __name__ == "__main__":
    main()