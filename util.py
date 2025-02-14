import pandas as pd
import json
import copy
import re
import os
import argparse


def calculate_q_values(df):
    """
    Calculates q-values for a target-decoy search.
    
    Parameters:
    df (pd.DataFrame): Input DataFrame containing protein accession and score columns.
    protein_column (str): The name of the column containing protein accession data (default: 'Protein accession').
    score_column (str): The name of the column containing identification scores (default: 'Score').
    decoy_identifier (str): The string that identifies decoy entries in the protein accession column (default: 'DECOY').
    
    Returns:
    pd.DataFrame: DataFrame with additional columns for cumulative decoy/target counts, FDR, and q-values.
    """
    # Copy the input DataFrame to avoid modifying the original data
    df = df.copy()

    # Add a column to indicate if the protein is a decoy
    df['IsDecoy'] = df["Protein accession"].str.contains("DECOY")

    # Sort by score (assuming higher score means better identification)
    df = df.sort_values(by="E-value")

    # Initialize counters for decoy and target counts
    df['Cumulative_Decoy'] = df['IsDecoy'].cumsum()
    df['Cumulative_Target'] = (~df['IsDecoy']).cumsum()

    # Calculate FDR: FDR = (# decoys / # total)
    df['FDR'] = df['Cumulative_Decoy'] / (df['Cumulative_Decoy'] + df['Cumulative_Target'])

    # Calculate q-value: the minimum FDR at or above this score
    df['q-value'] = df['FDR'][::-1].cummin()[::-1]  # Reverse cummin to get the minimum FDR for each score

    df.sort_index()

    return df["q-value"]

def read_tsv(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        
        if first_line == "********************** Parameters **********************":
            # Find the next occurrence of the delimiter line
            skip_rows = 1
            for line in file:
                if line.strip() == "********************** Parameters **********************":
                    skip_rows += 1
                    break
                skip_rows += 1
        else:
            skip_rows = 0  # No lines to skip
    
    # Read the TSV file with the appropriate skiprows value
    return pd.read_csv(file_path, sep='\t', skiprows=skip_rows)

def getMatchedPeaks(prsmID, dir, spec):
    with open(os.path.join(dir, "prsm" + str(prsmID) + ".js")) as file:
        file.readline()
        toppic = json.loads(file.read())
        peak_list = toppic["prsm"]["ms"]["peaks"]["peak"]
        matched_list = []
        nonmatched_list = []
        if len(spec.peak_list) == 1:
            matched_list.append(copy.deepcopy(spec.peak_list[0]))
        else:
            for idx in range(0, len(peak_list)):
                if "matched_ions" in peak_list[idx]:
                    matched_list.append(copy.deepcopy(spec.peak_list[idx]))
                else:
                    nonmatched_list.append(copy.deepcopy(spec.peak_list[idx]))
        return matched_list, nonmatched_list
    
def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts



def _gene_theo_ions(prot_seq):
    left_ions = []
    right_ions = []

    prot_seq = prot_seq.split(".", 1)[1]
    prot_seq = prot_seq.rsplit(".", 1)[0]
  
    acetylation = False

    if "[Acetyl]-" in prot_seq:
        acetylation = True
        prot_seq = prot_seq.replace('[Acetyl]-', '')

    if "(C)[Carbamidomethylation]" in prot_seq:
        prot_seq = prot_seq.replace("(C)[Carbamidomethylation]", "X")


    acetylation_weight = 42.0106
    # c57_weight = 57.021464

    weights = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841, 'G': 57.02146, 
            'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
            'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203, 'T': 101.04768, 'V': 99.06841,
            'W': 186.07931, 'Y': 163.06333, "X": 160.030654}


    if acetylation: 
        left_ions.append(acetylation_weight) 

    idx = 0
    previousweight = 0
    while idx < len(prot_seq):
        if (prot_seq[idx] == "("):
            startidx = idx
            endidx = prot_seq[idx:].find(")") + idx
            modS = endidx + 1
            modE = prot_seq[modS:].find("]") + modS
            weight = float(prot_seq[modS + 1:modE])
            frags = prot_seq[startidx + 1:endidx]
            for fragIdx in range(0, len(frags)):
                baseweight = previousweight
                left_ions.append(baseweight + weights[frags[fragIdx]] + weight)
                if fragIdx != len(frags) - 1:
                    left_ions.append(baseweight + weights[frags[fragIdx]])
                previousweight = left_ions[-1]
            idx = modE + 1
        else:
            left_ions.append(previousweight + weights[prot_seq[idx]])
            previousweight = left_ions[-1]
            idx += 1

    idx = len(prot_seq) - 1
    previousweight = 0
    while idx >= 0:
        if (prot_seq[idx] == "]"):
            startidx = prot_seq[:idx].rfind("(")
            endidx = prot_seq[:idx].rfind(")")
            modS = endidx + 1
            modE = idx
            weight = float(prot_seq[modS + 1:modE])
            frags = prot_seq[startidx + 1:endidx]
            for fragIdx in range(0, len(frags)):
                baseweight = previousweight
                right_ions.append(baseweight + weights[frags[len(frags) - fragIdx - 1]] + weight)
                if fragIdx != len(frags) - 1:
                    right_ions.append(baseweight + weights[frags[len(frags) - fragIdx - 1]])
                previousweight = right_ions[-1]
            idx = startidx - 1
        else:
            right_ions.append(previousweight + weights[prot_seq[idx]])
            previousweight = right_ions[-1]
            idx -= 1

    return left_ions, right_ions

def _get_modified_fragments(mass_list, shift):
    mod_mass_list = [x + shift for x in mass_list]
    return mod_mass_list

def _remove(peak_list, mass_list, shift):
    frags = _get_modified_fragments(mass_list, shift)
    count = 0
    for idx, peak in reversed(list(enumerate(peak_list))):
        peak_mass = peak.mass
        for j in range(len(frags)):
            frag_mass = frags[j]
            tol = (10 * frag_mass) / 1e6
            if (tol < 0.01):
                tol = 0.01
            if (abs(peak_mass - frag_mass) <= tol):
                del peak_list[idx]
                count += 1
                break
    return count


def removePeaks(peak_list, prot_sequence):
    bions, yions = _gene_theo_ions(prot_sequence)
    Proton = 1.00727647
    H = 1.007825035
    O = 15.99491463
    CO = 12.0000 + O
    NH3 = 14.003074 + H + H + H
    H2O = H + H + O
    
    count = 0
    # b -ion
    count += _remove(peak_list, bions, 0.0)
    # y -ion
    count += _remove(peak_list, yions, 19.0184-Proton)
    # # b - 1
    # count += remove(peak_list, bions, -Proton)
    # # y - 1
    # count += remove(peak_list, yions, 19.0184-Proton-Proton)
    # # b + 1
    # count += remove(peak_list, bions, Proton)
    # # y + 1
    # count += remove(peak_list, yions, 19.0184)

    # if (count > 0):
    #     print(str(count) + " peaks was removed from this spectra")

    return peak_list

def getProteoforms(inputdf, threshold=1.2, filterbyfeature = True):
    combined_df = inputdf.copy(deep=True)
    # Drop duplicates using feature IDs and keeping the one with the lowest E-value
    if filterbyfeature:
      combined_df = combined_df.sort_values(by='E-value').drop_duplicates(subset='Feature ID', keep='first').reset_index(drop=True)

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

    result_df = result_df.sort_values(by="Scan(s)").reset_index(drop=True)

    return result_df

def str_to_bool(value):
    """Convert a string to a boolean value."""
    if isinstance(value, bool):
        return value
    if value.lower() in ('true', '1', 'yes', 'y', 't'):
        return True
    elif value.lower() in ('false', '0', 'no', 'n', 'f'):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected (true/false, 1/0, yes/no, y/n).")