import sys
import pandas as pd
import json
import read_msalign

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

    merge = A.merge(B, on="Scan(s)", how="inner")

    switchlist = merge[merge["#matched peaks_y"] > merge["#matched peaks_x"]]["Scan(s)"].tolist()

    print("We have " + str(len(switchlist)) + " scans based on matched peaks comparsions")

    onlyB = set(B["Scan(s)"].tolist()) - set(A['Scan(s)'].tolist())

    print("We also have " + str(len(onlyB)) + " spectra who only reported under second precuror")

    switchlist = onlyB.union(set(switchlist))

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

    filteredB = prsmother[prsmother["Scan(s)"].isin(switchlist)]

    filteredA = prsm[~prsm["Scan(s)"].isin(switchlist)]

    result = pd.concat([filteredA, filteredB], ignore_index=True).sort_values(by="Scan(s)")

    result["Spectrum-level Q-value"] = calculate_q_values(result)

    outputresult = result[result["Spectrum-level Q-value"] < 0.01]

    outputresult.to_csv(args[0] + "FirstPrSM_toppic_prsm_single.tsv", sep="\t", index=False)

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
        filename = args[0] + "/" + str(dirmid) + dirend + "prsm" + str(prsmid) + ".js"
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


    read_msalign.write_spec_file(args[0] + "SecondPrSM_ms2.msalign", spectra)

if __name__ == "__main__":
    main()