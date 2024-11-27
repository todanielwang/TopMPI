import sys
import pandas as pd
import numpy as np
import read_msalign
import copy
import json

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

# Function to find duplicates based on the condition
def drop_custom_duplicates(group):
    threshold = 1.2

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


def main():
    #Please pass in directory with option flag of cutofftype and cutoffvalue
    args = sys.argv[1:]
    
    cutofftype = "FDR"
    cutoffvalue = 0.01

    if len(args) > 1:
        cutofftype = args[1]
        cutoffvalue = float(args[2])
    
    inputdf = pd.read_csv(args[0] + "Result_final.tsv", delimiter="\t")

    a_spec_list = read_msalign.read_spec_file(args[0] + "A_ms2.msalign")

    result_a = pd.read_csv(args[0] + "A_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

    result_ab = pd.read_csv(args[0] + "AB_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

    result_b = pd.read_csv(args[0] + "B_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

    result_ba = pd.read_csv(args[0] + "BA_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

    # ratio = 10
    # speclist = [spec for spec in a_spec_list if (len(spec.header.pre_mz_list) > 1) and (float(spec.header.pre_inte_list[1]) > 0) and (float(spec.header.pre_inte_list[0]) / float(spec.header.pre_inte_list[1]) < ratio)]

    # coverage = 0.7
    # speclist = [spec for spec in a_spec_list if (len(spec.header.pre_mz_list) > 1) and (float(spec.header.pre_inte_list[1]) > 0) and (float(spec.header.pre_inte_list[0]) + float(spec.header.pre_inte_list[1]) > coverage * sum(map(float, spec.header.pre_inte_list)))]

    # intensity = 5
    # speclist = [spec for spec in a_spec_list if (len(spec.header.pre_mz_list) > 1) and (float(spec.header.pre_inte_list[1]) > 0) and (float(spec.header.pre_inte_list[1]) > intensity * 10000)]

    # scanlist = [int(spec.header.spec_scan) for spec in speclist]

    scanlist = []

    output1 = pd.DataFrame(columns=result_a.columns).astype(result_a.dtypes)

    # inputdf["choice"] = np.where((inputdf["F1 Con"] == "True") & (inputdf["F2 Con"] == "True") & (inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] == inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]), "-", inputdf["choice"])

    for index, row in inputdf.iterrows():
        scan = row["Scan"]
        if int(scan) not in scanlist:
            o1 = pd.Series()

            A1 = result_a[result_a["Scan(s)"] == scan]
            B1 = result_b[result_b["Scan(s)"] == scan]
            
            if (row["choice"] == "-"):
                o1 = A1.iloc[0]
            elif ((row["choice"] == "A")):
                if A1.shape[0] > 0:
                    o1 = A1.iloc[0]
            elif (row["choice"] == "B"):
                if B1.shape[0] > 0:
                    o1 = B1.iloc[0]

            if not o1.empty:
                output1.loc[0 if pd.isnull(output1.index.max()) else output1.index.max() + 1] = o1
        # else:
        #     A1 = result_a[result_a["Scan(s)"] == scan]
        #     if A1.shape[0] > 0:
        #         output1.loc[0 if pd.isnull(output1.index.max()) else output1.index.max() + 1] = A1.iloc[0]

    output1["Data file name"] = "First Identification"

    output3 = output1.sort_values(by='E-value').drop_duplicates(subset='Feature ID', keep='first').groupby('Protein accession', group_keys=False).apply(drop_custom_duplicates)

    if cutofftype == "E-value":
        prsm1 = output1[output1["E-value"] < cutoffvalue]

        proteoform1 = output3[output3["E-value"] < cutoffvalue]
    else:
        output1.drop(["Spectrum-level Q-value", "Proteoform-level Q-value"], axis=1, inplace=True)
        output3.drop(["Spectrum-level Q-value", "Proteoform-level Q-value"], axis=1, inplace=True)

        output1["Spectrum-level Q-value"] = calculate_q_values(output1)
        output3["Proteoform-level Q-value"] = calculate_q_values(output3)

        #  & ~(output1["Protein accession"].str.contains("DECOY"))
        prsm1 = output1[(output1["Spectrum-level Q-value"] < cutoffvalue)].sort_values(by="Scan(s)")

        #  & ~(output3["Protein accession"].str.contains("DECOY"))
        proteoform1 = output3[(output3["Proteoform-level Q-value"] < cutoffvalue)].sort_values(by="Scan(s)")
    
    prsm1scanlist = prsm1["Scan(s)"].tolist()

    output2 = pd.DataFrame(columns=result_a.columns).astype(result_a.dtypes)
    for index, row in inputdf.iterrows():
        scan = row["Scan"]
        if int(scan) not in scanlist:
            o2 = pd.Series()

            A1 = result_a[result_a["Scan(s)"] == scan]
            A2 = result_ab[result_ab["Scan(s)"] == scan]
            B1 = result_b[result_b["Scan(s)"] == scan]
            B2 = result_ba[result_ba["Scan(s)"] == scan]
            
            if (row["choice"] == "-"):
                if (A2.iloc[0]["Feature score"] > 0.8):
                    if scan in prsm1scanlist:
                        o2 = A2.iloc[0]
                    else:
                        o2 = B1.iloc[0]
            elif (row["choice"] == "A"):
                if scan in prsm1scanlist:
                    if (A2.shape[0] > 0) & (A2.iloc[0]["Protein accession"] != A1.iloc[0]["Protein accession"]):
                        o2 = A2.iloc[0]
                else:
                    if B1.shape[0] > 0:
                        o2 = B1.iloc[0]
            elif (row["choice"] == "B"):
                if scan in prsm1scanlist:
                    if (B2.shape[0] > 0) & (B2.iloc[0]["Protein accession"] != B1.iloc[0]["Protein accession"]):
                        o2 = B2.iloc[0]
                else:
                    if A1.shape[0] > 0:
                        o2 = A1.iloc[0]

            # # Filter on R1R2 + R3R4
            # if (not o1.empty) and (not o2.empty) and (o1["Protein accession"] == o2["Protein accession"]):
            #     if o1["E-value"] <= o2["E-value"]:
            #         o2 = pd.Series()
            #     else:
            #         o1 = pd.Series()

            if not o2.empty:
                output2.loc[0 if pd.isnull(output2.index.max()) else output2.index.max() + 1] = o2
        # else:
        #     A2 = result_ab[result_ab["Scan(s)"] == scan]
        #     if A2.shape[0] > 0:
        #         output2.loc[0 if pd.isnull(output2.index.max()) else output2.index.max() + 1] = A2.iloc[0]

    output2["Data file name"] = "Second Identification"

    output4 = output2.sort_values(by='E-value').drop_duplicates(subset='Feature ID', keep='first').groupby('Protein accession', group_keys=False).apply(drop_custom_duplicates)

    if cutofftype == "E-value":
        prsm2 = output2[output2["E-value"] < cutoffvalue]

        proteoform2 = output4[output4["E-value"] < cutoffvalue]
    else:
        output2.drop(["Spectrum-level Q-value", "Proteoform-level Q-value"], axis=1, inplace=True)
        output4.drop(["Spectrum-level Q-value", "Proteoform-level Q-value"], axis=1, inplace=True)

        output2["Spectrum-level Q-value"] = calculate_q_values(output2)
        output4["Proteoform-level Q-value"] = calculate_q_values(output4)
        
        # & ~(output2["Protein accession"].str.contains("DECOY"))
        prsm2 = output2[(output2["Spectrum-level Q-value"] < cutoffvalue)].sort_values(by="Scan(s)")

        # & ~(output4["Protein accession"].str.contains("DECOY"))
        proteoform2 = output4[(output4["Proteoform-level Q-value"] < cutoffvalue)].sort_values(by="Scan(s)")

    prsm1 = prsm1[~(prsm1["Protein accession"].str.contains("DECOY"))]
    prsm2 = prsm2[~(prsm2["Protein accession"].str.contains("DECOY"))]
    proteoform1 = proteoform1[~(proteoform1["Protein accession"].str.contains("DECOY"))]
    proteoform2 = proteoform2[~(proteoform2["Protein accession"].str.contains("DECOY"))]

    prsm1.to_csv(args[0] + "tests/R1R2+R3R4/prsm1.tsv", sep="\t", index=False)
    prsm2.to_csv(args[0] + "tests/R1R2+R3R4/prsm2.tsv", sep="\t", index=False)
    proteoform1.to_csv(args[0] + "tests/R1R2+R3R4/proteoform1.tsv", sep="\t", index=False)
    proteoform2.to_csv(args[0] + "tests/R1R2+R3R4/proteoform2.tsv", sep="\t", index=False)        

if __name__ == "__main__":
    main()