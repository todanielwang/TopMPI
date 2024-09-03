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

    r1 = r1[~r1['Protein accession'].str.contains('DECOY')]
    r2 = r2[~r2['Protein accession'].str.contains('DECOY')]

    # Sort both DataFrames by E-value
    df1_sorted = r1.sort_values(by='E-value')
    df2_sorted = r2.sort_values(by='E-value')

    # Function to identify duplicates based on custom condition
    def identify_duplicates(df1, df2, threshold):
        set1 = set()
        set2 = set()
        overlap = set()

        # Compare each entry in df1 with df2
        for i, row1 in df1.iterrows():
            is_duplicate = False
            for j, row2 in df2.iterrows():
                if row1['Protein accession'] == row2['Protein accession'] and abs(row1['Precursor mass'] - row2['Precursor mass']) < threshold:
                    overlap.add((row1['Protein accession'], row1['Precursor mass']))
                    is_duplicate = True
                    break
            if not is_duplicate:
                set1.add((row1['Protein accession'], row1['Precursor mass']))

        # Compare each entry in df2 with df1 (not required to find overlaps again)
        for i, row2 in df2.iterrows():
            is_duplicate = False
            for j, row1 in df1.iterrows():
                if row2['Protein accession'] == row1['Protein accession'] and abs(row2['Precursor mass'] - row1['Precursor mass']) < threshold:
                    is_duplicate = True
                    break
            if not is_duplicate:
                set2.add((row2['Protein accession'], row2['Precursor mass']))
        
        return set1, set2, overlap

    # Get the sets for the Venn diagram
    set1, set2, overlap = identify_duplicates(df1_sorted, df2_sorted, 1.2)

    # Calculate the sizes for the Venn diagram
    only_set1 = len(set1)
    only_set2 = len(set2)
    intersection = len(overlap)

    # Plot the Venn diagram
    venn2(subsets=(only_set1, only_set2, intersection), set_labels=('Vanilla TopPIC', 'MSDeplex'))
    plt.savefig("proteoform.png", dpi=200)

if __name__ == "__main__":
    main()