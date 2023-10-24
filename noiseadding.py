import read_msalign
import pandas as pd
import sys
import re
import random
import math


def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts

def main():
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
        for counter in range(0, math.ceil(len(spec.peak_list)*2.0)):
            peaklist = spec_list[random.randint(0, len(spec_list) - 1)].peak_list
            peak = peaklist[random.randint(0, len(peaklist) - 1)]
            while (peak in spec.peak_list):
                peaklist = spec_list[random.randint(0, len(spec_list) - 1)].peak_list
                peak = peaklist[random.randint(0, len(peaklist) - 1)]
            spec.peak_list.append(peak)

    read_msalign.write_spec_file(args[0], spec_list)

if __name__ == "__main__":
    main()