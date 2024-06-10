import sys
import pandas as pd
import random
import json
import copy
import read_msalign
import math

def getMatchedPeaks(prsmID, dir, spec):
    with open(dir + "prsm" + str(prsmID) + ".js") as file:
        file.readline()
        toppic = json.loads(file.read())
        peak_list = toppic["prsm"]["ms"]["peaks"]["peak"]
        matched_list = []
        nonmatched_list = []
        for idx in range(0, len(peak_list)):
            if "matched_ions" in peak_list[idx]:
                matched_list.append(copy.deepcopy(spec.peak_list[idx]))
            else:
                nonmatched_list.append(copy.deepcopy(spec.peak_list[idx]))
        return matched_list, nonmatched_list

def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the directory")
    
    a_dir = args[0] + "A_html/toppic_prsm_cutoff/data_js/prsms/"

    ab_dir = args[0] + "AB_html/toppic_prsm_cutoff/data_js/prsms/"

    b_dir = args[0] + "B_html/toppic_prsm_cutoff/data_js/prsms/"

    ba_dir = args[0] + "BA_html/toppic_prsm_cutoff/data_js/prsms/"

    a_spec_list = read_msalign.read_spec_file(args[0] + "A_ms2.msalign")

    ab_spec_list = read_msalign.read_spec_file(args[0] + "AB_ms2.msalign")

    b_spec_list = read_msalign.read_spec_file(args[0] + "B_ms2.msalign")

    ba_spec_list = read_msalign.read_spec_file(args[0] + "BA_ms2.msalign")

    result = pd.read_csv(args[0] + "Result_resolved.tsv", delimiter="\t", index_col=0)

    a_result = pd.read_csv(args[0] + "A_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)

    ab_result = pd.read_csv(args[0] + "AB_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)

    b_result = pd.read_csv(args[0] + "B_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)

    ba_result = pd.read_csv(args[0] + "BA_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)

    pairs = []
    for spec in a_spec_list:
        pairs.append(frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)))

    accessA = {}
    for pair in pairs:
        for spec in a_spec_list:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessA[pair] = spec.header.spec_scan

    accessB = {}
    for pair in pairs:
        for spec in b_spec_list:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessB[pair] = spec.header.spec_scan
    spec_dict_a = {}
    for spec in a_spec_list:
        spec_dict_a[str(spec.header.spec_scan)] = spec

    spec_dict_ab = {}
    for spec in ab_spec_list:
        spec_dict_ab[str(spec.header.spec_scan)] = spec

    spec_dict_b = {}
    for spec in b_spec_list:
        spec_dict_b[str(spec.header.spec_scan)] = spec

    spec_dict_ba = {}
    for spec in ba_spec_list:
        spec_dict_ba[str(spec.header.spec_scan)] = spec

    output_list_1 = []
    output_list_2 = []
    output = []
    for pair in pairs:
        # If no prsm might fail, FIXME for later
        if (result[result["Pair"] == str(pair)].iloc[0]["choice"] == "A"):
            main_spec = copy.deepcopy((spec_dict_a[str(accessA[pair])]))
            main_matchedList, nonMatchedList = getMatchedPeaks(a_result[a_result["Scan(s)"] == int(accessA[pair])].iloc[0]["Prsm ID"], a_dir, main_spec)
            
            sub_spec = copy.deepcopy((spec_dict_ab[str(accessA[pair])]))
            sub_matchedList, noiseList = getMatchedPeaks(ab_result[ab_result["Scan(s)"] == int(accessA[pair])].iloc[0]["Prsm ID"], ab_dir, sub_spec)
        else:
            main_spec = copy.deepcopy((spec_dict_b[str(accessB[pair])]))
            main_matchedList, nonMatchedList = getMatchedPeaks(b_result[b_result["Scan(s)"] == int(accessB[pair])].iloc[0]["Prsm ID"], b_dir, main_spec)
            
            sub_spec = copy.deepcopy((spec_dict_ba[str(accessB[pair])]))
            sub_matchedList, noiseList = getMatchedPeaks(ba_result[ba_result["Scan(s)"] == int(accessB[pair])].iloc[0]["Prsm ID"], ba_dir, sub_spec)

        output.append(copy.deepcopy(main_spec))
        
        ratio = math.floor(len(main_matchedList) / (len(main_matchedList) + len(sub_matchedList)) * len(noiseList))
        random.shuffle(noiseList)
        main_spec.peak_list = main_matchedList + noiseList[:ratio]
        sub_spec.peak_list = sub_matchedList + noiseList[ratio:]
        output_list_1.append(main_spec)
        output_list_2.append(sub_spec)

    sorted1 = read_msalign.sortScans(output_list_1)

    sorted2 = read_msalign.sortScans(output_list_2)

    read_msalign.write_spec_file(args[0] + "temp_ms2.msalign", output)

    read_msalign.write_spec_file(args[0] + "resolved1_ms2.msalign", sorted1)

    read_msalign.write_spec_file(args[0] + "resolved2_ms2.msalign", sorted2)


if __name__ == "__main__":
    main()