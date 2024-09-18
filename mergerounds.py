import sys
import pandas as pd
import numpy as np
import csv

def main():
    args = sys.argv[1:]
    

    r1 = pd.read_csv(args[0], sep="\t")
    r2 = pd.read_csv(args[1], sep="\t")

    # Concatenate two proteoform files
    combined_df = pd.concat([r1, r2], ignore_index=True)

    combined_df = combined_df[~combined_df['Protein accession'].str.contains('DECOY')]

    # Define the threshold for the absolute difference in mass
    threshold = 1.2

    # Drop duplicates using feature IDs and keeping the one with the lowest E-value
    if args[2] == "Drop":
        combined_df = combined_df.sort_values(by='E-value').drop_duplicates(subset='Feature ID', keep='first')

    # Function to find duplicates based on the condition
    def drop_custom_duplicates(group):
        # Sort the group by E-value to prioritize rows with the lowest value in E-value
        group = group.sort_values(by='E-value')
        
        # Initialize a list to store indices of rows to keep
        keep_indices = []

        # Iterate through the sorted group
        for index, row in group.iterrows():
            # Check if this row is a duplicate of any previously kept row
            is_duplicate = False
            for keep_index in keep_indices:
                if abs(row['Precursor mass'] - group.loc[keep_index, 'Precursor mass']) < threshold:
                    is_duplicate = True
                    break
            # If not a duplicate, add it to the list of indices to keep
            if not is_duplicate:
                keep_indices.append(index)
        
        # Return only the rows to keep
        return group.loc[keep_indices]

    # Apply the function to groups defined by 'ColumnA'
    result_df = combined_df.groupby('Protein accession', group_keys=False).apply(drop_custom_duplicates)

    result_df.to_csv(args[0].rsplit("/", maxsplit=1)[0] + '/merged_results.tsv', sep='\t', index=False)
        


if __name__ == "__main__":
    main()