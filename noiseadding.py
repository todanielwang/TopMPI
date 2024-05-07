import read_msalign
import pandas as pd
import sys
import re
import random
import math
import copy


def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts

""" def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file and the prsm_single.tsv")
    
    random.seed(0)
    
    raw_spec_list = read_msalign.read_spec_file(args[0])

    prsm_df = pd.read_csv(args[1], delimiter="\t")

    scan_list = prsm_df["Scan(s)"].tolist()

    spec_list = []

    for spec in raw_spec_list:
        if int(spec.header.spec_scan) in scan_list:
            spec_list.append(spec)

    for spec in spec_list:
        for counter in range(0, math.ceil(len(spec.peak_list)*1.0)):
            peaklist = spec_list[random.randint(0, len(spec_list) - 1)].peak_list
            peak = peaklist[random.randint(0, len(peaklist) - 1)]
            while (peak in spec.peak_list):
                peaklist = spec_list[random.randint(0, len(spec_list) - 1)].peak_list
                peak = peaklist[random.randint(0, len(peaklist) - 1)]
            spec.peak_list.append(peak)

    read_msalign.write_spec_file(args[0], spec_list) """


def main():
    args = sys.argv[1:]
    if (len(args) != 3):
        raise Exception(
            "Please pass in the ms.align file, the prsm_single.tsv, and the amount of noise in decimal")
    
    random.seed(0)
    
    raw_spec_list = read_msalign.read_spec_file(args[0])

    prsm_df = pd.read_csv(args[1], delimiter="\t", skiprows=29)

    scan_list = prsm_df["Scan(s)"].tolist()

    prsms = []

    for spec in raw_spec_list:
        if int(spec.header.spec_scan) in scan_list:
            prsms.append(spec)

    nonmulti_dict = {}
    toplength = 0
    topprotein = ""
    topscan = 0
    for spec in prsms:
        if (len(spec.header.pre_inte_list) == 1 or float(spec.header.pre_inte_list[0]) > 0.8 * sum(map(float, spec.header.pre_inte_list))):
            nonmulti_dict[spec.header.spec_scan] = spec
            if len(spec.peak_list) > toplength:
                toplength = len(spec.peak_list)
                topscan = int(spec.header.spec_scan)
                topprotein = prsm_df[prsm_df["Scan(s)"] == int(spec.header.spec_scan)].iloc[0]["Protein accession"]

    print(len(nonmulti_dict), toplength, topscan)

    candidate_dict = {}
    topnoise = 1
    for spec in nonmulti_dict.values():
        if (len(spec.peak_list) * topnoise < toplength) and (prsm_df[prsm_df["Scan(s)"] == int(spec.header.spec_scan)].iloc[0]["Protein accession"] != topprotein): #To get rid of the longest one itself
            candidate_dict[spec.header.spec_scan] = spec
    print("We have total of " + str(len(candidate_dict)) + " candidates") 

    outputList = []

    for scan in candidate_dict.keys():
        spec = copy.deepcopy(candidate_dict[scan])
        protein = prsm_df[prsm_df["Scan(s)"] == int(scan)].iloc[0]["Protein accession"]
        length = len(spec.peak_list)
        numpeaksToAdd = math.ceil(length * float(args[2]))
        diff = 1.5
        possibleList = [x for x in nonmulti_dict.keys() if (x != scan) and (float(nonmulti_dict[x].header.pre_mz_list[0]) - float(spec.header.pre_mz_list[0]) < diff) and (len(nonmulti_dict[x].peak_list) >= numpeaksToAdd) and (prsm_df[prsm_df["Scan(s)"] == int(x)].iloc[0]["Protein accession"] != protein)]
        # (float(x.header.pre_mass_list[0]) - float(spec.header.pre_mass_list[0])) < 100 and 
        while (len(possibleList) == 0):
            diff *= 2
            possibleList = [x for x in nonmulti_dict.keys() if (x != scan) and (float(nonmulti_dict[x].header.pre_mz_list[0]) - float(spec.header.pre_mz_list[0]) < diff) and (len(nonmulti_dict[x].peak_list) >= numpeaksToAdd) and (prsm_df[prsm_df["Scan(s)"] == int(x)].iloc[0]["Protein accession"] != protein)]
        
        random_spec = nonmulti_dict[possibleList[random.randint(0, len(possibleList) - 1)]]

        other_scan = random_spec.header.spec_scan

        spec.header.title = str(other_scan)
        
        noiselist = random.sample(random_spec.peak_list, numpeaksToAdd)
        
        spec.peak_list = spec.peak_list + noiselist
        
        outputList.append(spec)


    read_msalign.write_spec_file(args[0], outputList)
    

if __name__ == "__main__":
    main()