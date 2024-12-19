import sys
import pandas as pd
import read_msalign
import json
import copy
import numpy as np

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
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the directory")
    A = pd.read_csv(args[0] + "/A_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)
    B = pd.read_csv(args[0] + "/B_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)

    spectra = read_msalign.read_spec_file(args[0] + "/A_ms2.msalign")

    spec_dict = {}
    for spec in spectra:
        spec_dict[str(spec.header.spec_scan)] = spec

    multiplexedspectra = [int(spec.header.spec_scan) for spec in spectra if (not (spec.header.pre_inte_list[0] == '' or float(spec.header.pre_inte_list[0]) == float(0) or len(spec.header.pre_inte_list) == 1 or float(spec.header.pre_inte_list[1]) == float(0))) and float(spec.header.pre_inte_list[1]) / float(spec.header.pre_inte_list[0]) > 0.2]

    print("We have " + str(len(multiplexedspectra)) + " multiplexed spectra out of total of " + str(len(spectra)) + " spectra")

    dirA = args[0] + "/A_html/toppic_prsm_cutoff/data_js/prsms/"
    dirB = args[0] + "/B_html/toppic_prsm_cutoff/data_js/prsms/"

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
            if (len(setA.intersection(setB)) / min(int(len(setA)), int(len(setB))) > 0.9):
                sameproteinscans.append(scan)

    print("We have total of " + str(len(sameproteinscans)) + " candidates that are incosistent")

    filteredmerge = merge[merge["Scan(s)"].isin(sameproteinscans)].copy()
    completemerge = filteredmerge[filteredmerge["#unexpected modifications_x"] - filteredmerge["#unexpected modifications_y"] >= 0].copy()
    completemerge.loc[:, "difference # ions"] = completemerge["#matched fragment ions_y"] - completemerge["#matched fragment ions_x"]
    completemerge.loc[:, "difference logged E-value"] = np.log10(completemerge["E-value_y"]) - np.log10(completemerge["E-value_x"])

    correctedscanlist = completemerge[completemerge["difference # ions"] >= 3]["Scan(s)"].tolist()

    print("We will switch precursor for " + str(len(correctedscanlist)) + " scans that are inconsistent")

    consistentscans = set(bothscanlist) - set(sameproteinscans)

    print("We have a total of " + str(len(consistentscans)) + " of scans that are consistent")

    consistentrows = merge[merge["Scan(s)"].isin(list(consistentscans))].copy()
    consistentswitchscans = consistentrows[consistentrows["#matched peaks_y"] > consistentrows["#matched peaks_x"]]["Scan(s)"].tolist()

    print("We will switch precursor for " + str(len(consistentswitchscans)) + " scans that are consistent")

    onlyB = set(multiplexB["Scan(s)"].tolist()) - set(multiplexA['Scan(s)'].tolist())

    print("We also have " + str(len(onlyB)) + " spectra who only reported under second precuror")

    switchlist = onlyB.union(set(correctedscanlist)).union(set(consistentswitchscans))

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

    read_msalign.write_spec_file(args[0] + "FirstPrSM_ms2.msalign", spectra)

    prsm = A
    prsmother = B

    # Filter df2 for rows where "Scan(s)" is in the scan_list
    filtered_other = prsmother[prsmother["Scan(s)"].isin(switchlist)]

    # Remove rows from df1 where "Scan(s)" is in the scan_list
    remaining_prsm = prsm[~prsm["Scan(s)"].isin(switchlist)]

    # Merge the dataframes, prioritizing df2 values (even if empty)
    result = pd.concat([remaining_prsm, filtered_other], ignore_index=True).sort_values(by="Scan(s)")

    #Change to E-value?
    result["Spectrum-level Q-value"] = calculate_q_values(result)

    result = result[result["Spectrum-level Q-value"] < 0.01]

    result.to_csv(args[0] + "FirstPrSM_toppic_prsm_single.tsv", sep="\t", index=False)

    splitspectra = [spec for spec in spectra if (len(spec.header.pre_mz_list) > 1 and int(spec.header.spec_scan) in multiplexedspectra)]

    for spec in splitspectra:
        spec.header.pre_mz_list = spec.header.pre_mz_list[1:]
        spec.header.pre_charge_list = spec.header.pre_charge_list[1:]
        spec.header.pre_mass_list = spec.header.pre_mass_list[1:]
        spec.header.pre_inte_list = spec.header.pre_inte_list[1:]
        spec.header.pre_id_list = spec.header.pre_id_list[1:]

    splitspecdict = {}
    for spec in splitspectra:
        splitspecdict[str(spec.header.spec_scan)] = spec

    dirend = "_html/toppic_prsm_cutoff/data_js/prsms/"
    count = 0
    for index, row in result.iterrows():
        dirmid = row["Data file name"].rsplit('/', 1)[-1].split('_')[0]
        prsmid = row["Prsm ID"]
        filename = args[0] + "/" + str(dirmid) + dirend + "prsm" + str(prsmid) + ".js"
        with open(filename) as file:
            file.readline()
            toppic = json.loads(file.read())
            scan = str(toppic["prsm"]["ms"]["ms_header"]["scans"])
            if scan in splitspecdict:
                peak_list = toppic["prsm"]["ms"]["peaks"]["peak"]
                for idx in range(len(peak_list) -1, -1, -1):
                    if "matched_ions" in peak_list[idx]:
                        del splitspecdict[scan].peak_list[idx]
                count += 1
    
    print("Number of scans with peaks removed is {}".format(count))


    read_msalign.write_spec_file(args[0] + "SecondPrSM_ms2.msalign", splitspectra)

if __name__ == "__main__":
    main()