import pandas as pd
import read_msalign
import json
import argparse
import util
import os

def main(args_list=None):
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
    args = parser.parse_args(args_list)

    print("Starting TopMPI precursor selection process with alpha " + str(args.alpha) + ", beta " + str(args.beta) + ", delta " + str(args.delta) + ", and gamma " + str(args.gamma))

    A = util.read_tsv(os.path.join(args.directory, "First_ms2_toppic_prsm_single.tsv"))
    B = util.read_tsv(os.path.join(args.directory, "Second_ms2_toppic_prsm_single.tsv"))

    spectra = read_msalign.read_spec_file(os.path.join(args.directory, "First_ms2.msalign"))

    spec_dict = {}
    for spec in spectra:
        spec_dict[str(spec.header.spec_scan)] = spec

    multiplexedspectra = [int(spec.header.spec_scan) for spec in spectra if (not (spec.header.pre_inte_list[0] == '' or float(spec.header.pre_inte_list[0]) == float(0) or len(spec.header.pre_inte_list) == 1 or float(spec.header.pre_inte_list[1]) == float(0))) and float(spec.header.pre_inte_list[1]) / float(spec.header.pre_inte_list[0]) > args.alpha]

    print("We have " + str(len(multiplexedspectra)) + " multiplexed spectra out of total of " + str(len(spectra)) + " spectra")


    dirA = os.path.join(args.directory, "First_html", "toppic_prsm_cutoff", "data_js", "prsms")
    dirB = os.path.join(args.directory, "Second_html", "toppic_prsm_cutoff", "data_js", "prsms")

    multiplexA = A[A["Scan(s)"].isin(multiplexedspectra)]
    multiplexB = B[B["Scan(s)"].isin(multiplexedspectra)]

    merge = multiplexA.merge(multiplexB, on="Scan(s)", how="inner")
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

    read_msalign.write_spec_file(os.path.join(args.directory, "Primary_ms2.msalign"), spectra)

    prsm = A
    prsmother = B

    filteredB = prsmother[prsmother["Scan(s)"].isin(switchlist)]

    filteredA = prsm[~prsm["Scan(s)"].isin(switchlist)]

    result = pd.concat([filteredA, filteredB], ignore_index=True).sort_values(by="Scan(s)")

    result.to_csv(os.path.join(args.directory, "Primary_ms2_temp_prsm.tsv"), sep="\t", index=False)

    for spec in spectra:
        spec.header.pre_mz_list = spec.header.pre_mz_list[1:]
        spec.header.pre_charge_list = spec.header.pre_charge_list[1:]
        spec.header.pre_mass_list = spec.header.pre_mass_list[1:]
        spec.header.pre_inte_list = spec.header.pre_inte_list[1:]
        spec.header.pre_id_list = spec.header.pre_id_list[1:]


    dirend = os.path.join("_html", "toppic_prsm_cutoff", "data_js", "prsms")
    count = 0
    for _, row in result.iterrows():
        dirmid = os.path.basename(row["Data file name"]).split('_')[0]
        prsmid = row["Prsm ID"]
        filename = os.path.join(args.directory, str(dirmid) + dirend, f"prsm{prsmid}.js")
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

    read_msalign.write_spec_file(os.path.join(args.directory, "Secondary_ms2.msalign"), spectra)

    print("Finished generating Primary_ms2.msalign and Secondary_ms2.msalign files")

if __name__ == "__main__":
    main()