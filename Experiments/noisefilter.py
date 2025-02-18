import sys
import pandas as pd
import json
import read_msalign
import util

def main():
    args = sys.argv[1:]
    if (len(args) != 3):
        raise Exception(
            "Please pass in the directory, delta, and gamma")
    A = pd.read_csv(args[0] + "/A_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)
    B = pd.read_csv(args[0] + "/B_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)

    spectra = read_msalign.read_spec_file(args[0] + "/A_ms2.msalign")

    spec_dict = {}
    for spec in spectra:
        spec_dict[str(spec.header.spec_scan)] = spec

    dirA = args[0] + "/A_html/toppic_prsm_cutoff/data_js/prsms/"
    dirB = args[0] + "/B_html/toppic_prsm_cutoff/data_js/prsms/"

    merge = A.merge(B, on="Scan(s)", how="inner")

    bothscanlist = merge["Scan(s)"].tolist()

    sameproteinscans = merge[merge["Protein accession_x"] == merge["Protein accession_y"]]["Scan(s)"].tolist()

    print("We have " + str(len(sameproteinscans)) + " scans where both precursor reported the same protein")

    for spec in spectra:
        scan = int(spec.header.spec_scan)
        if scan not in sameproteinscans and scan in bothscanlist:
            matchedListA, nonMatchedList = util.getMatchedPeaks(merge[merge["Scan(s)"] == int(scan)].iloc[0]["Prsm ID_x"], dirA, spec_dict[str(scan)])
            matchedListB, nonMatchedList = util.getMatchedPeaks(merge[merge["Scan(s)"] == int(scan)].iloc[0]["Prsm ID_y"], dirB, spec_dict[str(scan)])
            setA = set(matchedListA)
            setB = set(matchedListB)
            if (len(setA.intersection(setB)) / min(int(len(setA)), int(len(setB))) > 0.9):
                sameproteinscans.append(scan)

    print("We have total of " + str(len(sameproteinscans)) + " candidates that are incosistent")

    filteredmerge = merge[merge["Scan(s)"].isin(sameproteinscans)].copy().reset_index(drop=True)

    filteredmerge["difference"] = filteredmerge["#matched peaks_y"] - filteredmerge["#matched peaks_x"] - int(args[2])

    filteredmerge["balanceddifference"] = filteredmerge["difference"]

    weight = int(args[1])

    filteredmerge.loc[filteredmerge["#unexpected modifications_x"] > filteredmerge["#unexpected modifications_y"], "balanceddifference"] += weight

    filteredmerge.loc[filteredmerge["#unexpected modifications_x"] < filteredmerge["#unexpected modifications_y"], "balanceddifference"] -= weight

    switchlist = filteredmerge[filteredmerge["balanceddifference"] >= 0]["Scan(s)"].tolist()

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

    result["Spectrum-level Q-value"] = util.calculate_q_values(result)

    outputresult = result[result["Spectrum-level Q-value"] < 0.01]

    # outputresult = outputresult[~outputresult["Protein accession"].str.contains("DECOY")]

    outputresult.to_csv(args[0] + "FirstPrSM_toppic_prsm_single.tsv", sep="\t", index=False)

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

    proteoformresult["Proteoform-level Q-value"] = util.calculate_q_values(proteoformresult)

    # proteoformresult = proteoformresult[~proteoformresult['Protein accession'].str.contains('DECOY')].reset_index(drop=True)

    proteoformresult[proteoformresult["Proteoform-level Q-value"] < 0.01].to_csv(args[0] + "FirstPrSM_toppic_proteoform_single.tsv", sep="\t", index=False)

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