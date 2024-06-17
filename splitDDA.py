import read_msalign
import sys
import os
import json
import re

def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts

def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the ms.align file")
    
    jsfolder = args[0].rsplit("_", 1)[0] + "_html/toppic_proteoform_cutoff/data_js/prsms/"
    
    spec_list = read_msalign.read_spec_file(args[0])

    spec_list = [spec for spec in spec_list if len(spec.header.pre_mz_list) > 1]

    for spec in spec_list:
        spec.header.pre_mz_list = spec.header.pre_mz_list[1:]
        spec.header.pre_charge_list = spec.header.pre_charge_list[1:]
        spec.header.pre_mass_list = spec.header.pre_mass_list[1:]
        spec.header.pre_inte_list = spec.header.pre_inte_list[1:]
        spec.header.pre_id_list = spec.header.pre_id_list[1:]

    curr_spec = 0
    count = 0
    for x in sorted(os.listdir(jsfolder), key=numericalSort):
        if x.endswith(".js"):
            with open(jsfolder + x) as file:
                file.readline()
                toppic = json.loads(file.read())
                deleted = False
                while(curr_spec < len(spec_list)):
                    filescan = spec_list[curr_spec].header.spec_scan
                    toppicscan = int(toppic["prsm"]["ms"]["ms_header"]["scans"])
                    if (filescan < toppicscan):
                        curr_spec += 1
                    else:
                        if (filescan > toppicscan):
                            deleted = True
                        break
                if (curr_spec >= len(spec_list)):
                    break
                if (deleted):
                    continue
                peak_list = toppic["prsm"]["ms"]["peaks"]["peak"]
                for idx in range(len(peak_list) -1, -1, -1):
                    if "matched_ions" in peak_list[idx]:
                        del spec_list[curr_spec].peak_list[idx]
                count += 1
    print("Number of scans with peaks removed is {}".format(count))


    read_msalign.write_spec_file(args[0], spec_list)

if __name__ == "__main__":
    main()