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


def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the directory")
    
    inputdf = pd.read_csv(args[0] + "Result.tsv", delimiter="\t", index_col=0)

    a_dir = args[0] + "A_html/toppic_prsm_cutoff/data_js/prsms/"

    ab_dir = args[0] + "AB_html/toppic_prsm_cutoff/data_js/prsms/"

    b_dir = args[0] + "B_html/toppic_prsm_cutoff/data_js/prsms/"

    ba_dir = args[0] + "BA_html/toppic_prsm_cutoff/data_js/prsms/"

    a_spec_list = read_msalign.read_spec_file(args[0] + "A_ms2.msalign")

    ab_spec_list = read_msalign.read_spec_file(args[0] + "AB_ms2.msalign")

    b_spec_list = read_msalign.read_spec_file(args[0] + "B_ms2.msalign")

    ba_spec_list = read_msalign.read_spec_file(args[0] + "BA_ms2.msalign")

    result_a = pd.read_csv(args[0] + "A_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

    result_ab = pd.read_csv(args[0] + "AB_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

    result_b = pd.read_csv(args[0] + "B_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

    result_ba = pd.read_csv(args[0] + "BA_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=26)

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

    inputdf["F1 Con"] = "-"
    inputdf["F2 Con"] = "-"

    count = 0
    for index, row in inputdf.iterrows():
        scan = row["Scan"]

        if row["A+B_1"] == "-":
            matchedList_a = []
        else:
            main_spec = copy.deepcopy((spec_dict_a[str(scan)]))
            matchedList_a, nonMatchedList = getMatchedPeaks(result_a[result_a["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], a_dir, main_spec)

        if row["A+B_2"] == "-":
            matchedList_ab = []
        else:    
            sub_spec = copy.deepcopy((spec_dict_ab[str(scan)]))
            matchedList_ab, noiseList = getMatchedPeaks(result_ab[result_ab["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], ab_dir, sub_spec)

        if row["B+A_1"] == "-":
            matchedList_b = []
        else:
            main_spec = copy.deepcopy((spec_dict_b[str(scan)]))
            matchedList_b, nonMatchedList = getMatchedPeaks(result_b[result_b["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], b_dir, main_spec)
            
        if row["B+A_2"] == "-":
            matchedList_ba = []
        else:
            sub_spec = copy.deepcopy((spec_dict_ba[str(scan)]))
            matchedList_ba, noiseList = getMatchedPeaks(result_ba[result_ba["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], ba_dir, sub_spec)
    
        seta = set(matchedList_a)
        setab = set(matchedList_ab)
        setb = set(matchedList_b)
        setba = set(matchedList_ba)


        if (len(seta) > 0 and len(setba) > 0 and len(seta.intersection(setba)) / min(float(len(seta)), float(len(setba))) > 0.9):
            inputdf.loc[(inputdf["Scan"] == str(scan)), ["F1 Con"]] = "True"
            count += 1
            continue
        if (len(setb) > 0 and len(setab) > 0 and len(setb.intersection(setab)) / min(float(len(setb)), float(len(setab))) > 0.9):
            inputdf.loc[(inputdf["Scan"] == str(scan)), ["F2 Con"]] = "True"
            count += 1
            continue
    
    print(count)

    inputdf["choice"] = "-"

    mask = (inputdf['F1 Con'] == "True") & (inputdf['F2 Con'] == "-")
    inputdf.loc[mask, "choice"] = "A" 

    mask = (inputdf['F1 Con'] == "-") & (inputdf['F2 Con'] == "True")
    inputdf.loc[mask, "choice"] = "B"

    
    mask = (((inputdf['F1 Con'] == "True") & (inputdf['F2 Con'] == "True")) | ((inputdf['F1 Con'] == "-") & (inputdf['F2 Con'] == "-"))) & ((inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] > inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) | ((inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] == inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) & (inputdf["A+B_1 E-value"] <= inputdf["B+A_1 E-value"])))
    inputdf.loc[mask, "choice"] = "A"

    mask = (((inputdf['F1 Con'] == "True") & (inputdf['F2 Con'] == "True")) | ((inputdf['F1 Con'] == "-") & (inputdf['F2 Con'] == "-"))) & ((inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] < inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) | ((inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] == inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) & (inputdf["A+B_1 E-value"] > inputdf["B+A_1 E-value"])))
    inputdf.loc[mask, "choice"] = "B"

    inputdf.to_csv(args[0] + "Result_final.tsv", sep="\t")

    output_list_1 = []
    output_list_2 = []

    for index, row in inputdf.iterrows():
        scan = row["Scan"]
        if (row["choice"] == "A"):
            if (row["A+B_2"] == "-"):
                if (row["A+B_1"] == "-"):
                    continue
                output_list_1.append(spec_dict_a[str(scan)])
                continue
            elif (row["A+B_1"] == "-"):
                if (row["A+B_2"] == "-"):
                    continue
                output_list_2.append(spec_dict_ab[str(scan)])
                continue
            
            main_spec = copy.deepcopy((spec_dict_a[str(scan)]))
            # main_matchedList, nonMatchedList = getMatchedPeaks(result_a[result_a["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], a_dir, main_spec)
            
            sub_spec = copy.deepcopy((spec_dict_ab[str(scan)]))
            # sub_matchedList, noiseList = getMatchedPeaks(result_ab[result_ab["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], ab_dir, sub_spec)
        elif (row["choice"] == "B"):
            if (row["B+A_2"] == "-"):
                if (row["B+A_1"] == "-"):
                    continue
                output_list_1.append(spec_dict_b[str(scan)])
                continue
            elif (row["B+A_1"] == "-"):
                if (row["B+A_2"] == "-"):
                    continue
                output_list_2.append(spec_dict_ba[str(scan)])
                continue

            main_spec = copy.deepcopy((spec_dict_b[str(scan)]))
            # main_matchedList, nonMatchedList = getMatchedPeaks(result_b[result_b["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], b_dir, main_spec)
            
            sub_spec = copy.deepcopy((spec_dict_ba[str(scan)]))
            # sub_matchedList, noiseList = getMatchedPeaks(result_ba[result_ba["Scan(s)"] == int(scan)].iloc[0]["Prsm ID"], ba_dir, sub_spec)

        # output.append(copy.deepcopy(main_spec))

        # main_spec.peak_list = main_matchedList + noiseList
        # sub_spec.peak_list = sub_matchedList + noiseList
        output_list_1.append(main_spec)
        output_list_2.append(sub_spec)

    # print(len(output_list_1), len(output_list_2))

    sorted1 = read_msalign.sortScans(output_list_1)
    sorted2 = read_msalign.sortScans(output_list_2)

    # for spec in sorted2:
    #     spec.header.spec_scan += 10000000
    #     spec.header.spec_id += 10000000

    # read_msalign.write_spec_file(args[0] + "temp_ms2.msalign", output)

    read_msalign.write_spec_file(args[0] + "resolved1_ms2.msalign", sorted1)
    read_msalign.write_spec_file(args[0] + "resolved2_ms2.msalign", sorted2)

        

if __name__ == "__main__":
    main()