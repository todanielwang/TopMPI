#Map the scans with first feature
import numpy as np
import sys
import pandas as pd
from matplotlib import pyplot as plt

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception()
    
    df = pd.read_csv(args[0], delimiter="\t")
    with open(args[1]) as f:
        scanList = [np.int64(x) for x in filter(None, f.read().split("\n"))]

    after = df[df["Scan(s)"].isin(scanList)]

    after["Precursor mass"].plot.hist(alpha=.5, legend=True)

    plt.savefig('Mass.png')

    plt.close()

    after["Charge"].plot.hist(alpha=.5, legend=True)

    plt.savefig('Charge.png')

    #after.to_csv("MissingPrSMs.csv", sep="\t", index=False)

    

if __name__ == "__main__":
    main()