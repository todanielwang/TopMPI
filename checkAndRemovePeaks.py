import pandas as pd
import read_msalign
import json
import copy
import argparse

def getMatchedPeaks(prsmID, dir, spec):
    with open(dir + "prsm" + str(prsmID) + ".js") as file:
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

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="The main body of TopMPI, use -h to see manual")

    # Add the positional argument
    parser.add_argument("directory", help="The directory to the TopMPI folder")

    # Add the optional flags
    parser.add_argument("-a", "--alpha", type=float, default=0.2, help="The intensity ratio between the first and second precursor required for a spectrum to treated as multiplexed")
    parser.add_argument("-b", "--beta", type=float, default=0.9, help="The percentage of shared matched peaks required for the identifications of the 2 precursors to be treated as inconsistent")
    parser.add_argument("-d", "--delta", type=int, default=5, help="The offset to the # of matched peaks condition based on number of unexpected mass shifts")
    parser.add_argument("-g", "--gamma", type=int, default=4, help="The number of matched peaks difference required to switch to the second precursor")

    # Parse the arguments
    args = parser.parse_args()

    A = pd.read_csv(args.directory + "/A_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=26)
    B = pd.read_csv(args.directory + "/B_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=26)

    spectra = read_msalign.read_spec_file(args.directory + "/A_ms2.msalign")

    spec_dict = {}
    for spec in spectra:
        spec_dict[str(spec.header.spec_scan)] = spec

    multiplexedspectra = [int(spec.header.spec_scan) for spec in spectra if (not (spec.header.pre_inte_list[0] == '' or float(spec.header.pre_inte_list[0]) == float(0) or len(spec.header.pre_inte_list) == 1 or float(spec.header.pre_inte_list[1]) == float(0))) and float(spec.header.pre_inte_list[1]) / float(spec.header.pre_inte_list[0]) > args.alpha]

    print("We have " + str(len(multiplexedspectra)) + " multiplexed spectra out of total of " + str(len(spectra)) + " spectra")

    dirA = args.directory + "/A_html/toppic_prsm_cutoff/data_js/prsms/"
    dirB = args.directory + "/B_html/toppic_prsm_cutoff/data_js/prsms/"

    multiplexA = A[A["Scan(s)"].isin(multiplexedspectra)]
    multiplexB = B[B["Scan(s)"].isin(multiplexedspectra)]

    merge = multiplexA.merge(multiplexB, on="Scan(s)", how="inner")
    bothscanlist = merge["Scan(s)"].tolist()

    sameproteinscans = merge[merge["Protein accession_x"] == merge["Protein accession_y"]]["Scan(s)"].tolist()

    print("We have " + str(len(sameproteinscans)) + " scans where both precursor reported the same protein")

    for spec in spectra:
        scan = int(spec.header.spec_scan)
        if scan not in sameproteinscans and scan in bothscanlist:
            matchedListA, nonMatchedList = getMatchedPeaks(merge[merge["Scan(s)"] == int(scan)].iloc[0]["Prsm ID_x"], dirA, spec_dict[str(scan)])
            matchedListB, nonMatchedList = getMatchedPeaks(merge[merge["Scan(s)"] == int(scan)].iloc[0]["Prsm ID_y"], dirB, spec_dict[str(scan)])
            setA = set(matchedListA)
            setB = set(matchedListB)
            if (len(setA.intersection(setB)) / min(int(len(setA)), int(len(setB))) > args.beta):
                sameproteinscans.append(scan)

    print("We have total of " + str(len(sameproteinscans)) + " candidates that are incosistent")

    filteredmerge = merge[merge["Scan(s)"].isin(sameproteinscans)].copy().reset_index(drop=True)

    filteredmerge["difference"] = filteredmerge["#matched peaks_y"] - filteredmerge["#matched peaks_x"] - int(args.gamma)

    filteredmerge["balanceddifference"] = filteredmerge["difference"]

    weight = int(args.delta)

    filteredmerge.loc[filteredmerge["#unexpected modifications_x"] > filteredmerge["#unexpected modifications_y"], "balanceddifference"] += weight

    filteredmerge.loc[filteredmerge["#unexpected modifications_x"] < filteredmerge["#unexpected modifications_y"], "balanceddifference"] -= weight

    correctedscanlist = filteredmerge[filteredmerge["balanceddifference"] >= 0]["Scan(s)"].tolist()

    print("We have " + str(len(correctedscanlist)) + " scans based on matched peaks comparsions")

    onlyB = set(multiplexB["Scan(s)"].tolist()) - set(multiplexA['Scan(s)'].tolist())

    print("We also have " + str(len(onlyB)) + " spectra who only reported under second precuror")

    switchlist = onlyB.union(set(correctedscanlist))

    print("We will be switching precursor for a total of " + str(len(switchlist)) + " scans")

    for spec in spectra:
        if int(spec.header.spec_scan) in switchlist:
            pre_mz_list = spec.header.pre_mz_list
            pre_charge_list = spec.header.pre_charge_list
            pre_mass_list = spec.header.pre_mass_list
            pre_inte_list = spec.header.pre_inte_list
            pre_id_list = spec.header.pre_id_list

            pre_mz_list[0], pre_mz_list[1] = pre_mz_list[1], pre_mz_list[0]
            pre_charge_list[0], pre_charge_list[1] = pre_charge_list[1], pre_charge_list[0]
            pre_mass_list[0], pre_mass_list[1] = pre_mass_list[1], pre_mass_list[0]
            pre_inte_list[0], pre_inte_list[1] = pre_inte_list[1], pre_inte_list[0]
            pre_id_list[0], pre_id_list[1] = pre_id_list[1], pre_id_list[0]

            spec.header.pre_mz_list = pre_mz_list
            spec.header.pre_charge_list = pre_charge_list
            spec.header.pre_mass_list = pre_mass_list
            spec.header.pre_inte_list = pre_inte_list
            spec.header.pre_id_list = pre_id_list

    read_msalign.write_spec_file(args.directory + "Primary_ms2.msalign", spectra)

    prsm = A
    prsmother = B

    filteredB = prsmother[prsmother["Scan(s)"].isin(switchlist)]

    filteredA = prsm[~prsm["Scan(s)"].isin(switchlist)]

    result = pd.concat([filteredA, filteredB], ignore_index=True).sort_values(by="Scan(s)")

    result["Spectrum-level Q-value"] = calculate_q_values(result)

    outputresult = result[result["Spectrum-level Q-value"] < 0.01]

    outputresult = outputresult[~outputresult["Protein accession"].str.contains("DECOY")]

    outputresult.to_csv(args.directory + "Primary_ms2_toppic_prsm_single.tsv", sep="\t", index=False)

    proteoformoutput = result

    # Define the threshold for the absolute difference in mass
    threshold = 1.2

    # Drop duplicates using feature IDs and keeping the one with the lowest E-value
    proteoformoutput = proteoformoutput.sort_values(by='E-value').drop_duplicates(subset='Feature ID', keep='first')

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
    proteoformresult = proteoformoutput.groupby('Protein accession', group_keys=False).apply(drop_custom_duplicates)

    proteoformresult["Proteoform-level Q-value"] = calculate_q_values(proteoformresult)

    proteoformresult = proteoformresult[~proteoformresult['Protein accession'].str.contains('DECOY')].reset_index(drop=True)

    proteoformresult[proteoformresult["Proteoform-level Q-value"] < 0.01].to_csv(args.directory + "Primary_ms2_toppic_proteoform_single.tsv", sep="\t", index=False)

    for spec in spectra:
        spec.header.pre_mz_list = spec.header.pre_mz_list[1:]
        spec.header.pre_charge_list = spec.header.pre_charge_list[1:]
        spec.header.pre_mass_list = spec.header.pre_mass_list[1:]
        spec.header.pre_inte_list = spec.header.pre_inte_list[1:]
        spec.header.pre_id_list = spec.header.pre_id_list[1:]

    dirend = "_html/toppic_prsm_cutoff/data_js/prsms/"
    count = 0
    for index, row in result.iterrows():
        dirmid = row["Data file name"].rsplit('/', 1)[-1].split('_')[0]
        prsmid = row["Prsm ID"]
        filename = args.directory + "/" + str(dirmid) + dirend + "prsm" + str(prsmid) + ".js"
        with open(filename) as file:
            file.readline()
            toppic = json.loads(file.read())
            scan = str(toppic["prsm"]["ms"]["ms_header"]["scans"])
            peak_list = toppic["prsm"]["ms"]["peaks"]["peak"]
            for idx in range(len(peak_list) -1, -1, -1):
                if "matched_ions" in peak_list[idx]:
                    del spec_dict[scan].peak_list[idx]
            count += 1
    
    print("Number of scans with peaks removed is {}".format(count))


    read_msalign.write_spec_file(args.directory + "Secondary_ms2.msalign", spectra)

if __name__ == "__main__":
    main()