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
    
    r1 = pd.read_csv(args[0])
    r2 = pd.read_csv(args[1])
    r3 = pd.read_csv(args[2])

    r1 = r1[~r1['Protein accession'].str.contains('DECOY')]
    r2 = r2[~r2['Protein accession'].str.contains('DECOY')]
    r3 = r3[~r3['Protein accession'].str.contains('DECOY')]

    df1_sorted = r1.sort_values(by='E-value')
    df2_sorted = r2.sort_values(by='E-value')
    df3_sorted = r3.sort_values(by='E-value')

    # Function to identify duplicates based on custom condition
    def identify_duplicates(df1, df2, df3, threshold):
        set1 = set()
        set2 = set()
        set3 = set()
        overlap12 = set()
        overlap13 = set()
        overlap23 = set()
        overlap_all = set()

        # Compare each entry in df1 with df2 and df3
        for i, row1 in df1.iterrows():
            is_duplicate12 = False
            is_duplicate13 = False

            for j, row2 in df2.iterrows():
                if row1['Protein accession'] == row2['Protein accession'] and abs(row1['Precursor mass'] - row2['Precursor mass']) < threshold:
                    overlap12.add((row1['Protein accession'], row1['Precursor mass']))
                    is_duplicate12 = True
                    break

            for k, row3 in df3.iterrows():
                if row1['Protein accession'] == row3['Protein accession'] and abs(row1['Precursor mass'] - row3['Precursor mass']) < threshold:
                    overlap13.add((row1['Protein accession'], row1['Precursor mass']))
                    is_duplicate13 = True
                    break

            if is_duplicate12 and is_duplicate13:
                overlap_all.add((row1['Protein accession'], row1['Precursor mass']))
            elif not is_duplicate12 and not is_duplicate13:
                set1.add((row1['Protein accession'], row1['Precursor mass']))

        # Compare each entry in df2 with df1 and df3 using the duplicate condition
        for i, row2 in df2.iterrows():
            is_duplicate_with_df1 = any(
                row2['Protein accession'] == row1['Protein accession'] and abs(row2['Precursor mass'] - row1['Precursor mass']) < threshold
                for j, row1 in df1.iterrows()
            )

            if is_duplicate_with_df1:
                if (row2['Protein accession'], row2['Precursor mass']) not in overlap_all:
                    continue

            is_duplicate23 = False

            for j, row3 in df3.iterrows():
                if row2['Protein accession'] == row3['Protein accession'] and abs(row2['Precursor mass'] - row3['Precursor mass']) < threshold:
                    overlap23.add((row2['Protein accession'], row2['Precursor mass']))
                    is_duplicate23 = True
                    break

            if not is_duplicate23 and not is_duplicate_with_df1:
                set2.add((row2['Protein accession'], row2['Precursor mass']))

        # Compare each entry in df3 with df1 and df2 using the duplicate condition
        for i, row3 in df3.iterrows():
            is_duplicate_with_df1 = any(
                row3['Protein accession'] == row1['Protein accession'] and abs(row3['Precursor mass'] - row1['Precursor mass']) < threshold
                for j, row1 in df1.iterrows()
            )

            is_duplicate_with_df2 = any(
                row3['Protein accession'] == row2['Protein accession'] and abs(row3['Precursor mass'] - row2['Precursor mass']) < threshold
                for j, row2 in df2.iterrows()
            )

            if not is_duplicate_with_df1 and not is_duplicate_with_df2:
                set3.add((row3['Protein accession'], row3['Precursor mass']))

        return set1, set2, set3, overlap12, overlap13, overlap23, overlap_all
    
    set1, set2, set3, overlap12, overlap13, overlap23, overlap_all = identify_duplicates(df1_sorted, df2_sorted, df3_sorted, 1.2)

    # Calculate the sizes for the Venn diagram
    only_set1 = len(set1)
    only_set2 = len(set2)
    only_set3 = len(set3)
    intersection12 = len(overlap12) - len(overlap_all)
    intersection13 = len(overlap13) - len(overlap_all)
    intersection23 = len(overlap23)
    intersection123 = len(overlap_all)

    # Plot the Venn diagram
    venn3(subsets=(only_set1, only_set2, intersection12, only_set3, intersection13, intersection23, intersection123),
        set_labels=('Replicate 1', 'Replicate 2', 'Replicate 3'))

    plt.savefig("proteoformoverlap.png", dpi=200)

    plt.show()

if __name__ == "__main__":
    main()