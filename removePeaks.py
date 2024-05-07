import read_msalign
import pandas as pd
import sys
import re
import random
import math
import copy
import os
import json


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
            "Please pass in the ms.align file, the prsm dir, and the decimal % amount of peaks to remove")
    
    random.seed(0)
    
    raw_spec_list = read_msalign.read_spec_file(args[0])

    for index, filename in enumerate(sorted(os.listdir(args[1]), key=numericalSort)):
        try:
            with open(args[1] + filename) as file:
                file.readline()
                toppic = json.loads(file.read())
        except FileNotFoundError:
            print("Wrong directory, file not found!")

        matched = []
        non_matched = []
        
        peak_list = toppic["prsm"]["ms"]["peaks"]["peak"]
        for idx in range(0, len(peak_list)):
            if "matched_ions" in peak_list[idx]:
                matched.append(idx)
            else:
                non_matched.append(idx)

        random_idx = sorted(random.sample(matched, int(len(matched) * (1 - float(args[2])))) + random.sample(non_matched, int(len(non_matched) * (1 - float(args[2])))))

        new_peak_list = []
        old_peak_list = copy.deepcopy(raw_spec_list[index].peak_list)
        for idx in random_idx:
            new_peak_list.append(old_peak_list[idx])

        raw_spec_list[index].peak_list = new_peak_list

    read_msalign.write_spec_file(args[0], raw_spec_list)


if __name__ == "__main__":
    main()