# There are 2 inputs to this function, one the directory to the proteoform.js files, usually under topsearchfolder_html/toppic_proteoform_cutoff/data_js/proteoform/
#The second is the exact protein accession that you are looking for
#example: I want to look at ecoli.mzml and protein 123, the command I use would be
#python3 ubq.py DIRECTORY_TO_THE_NEXT_FOLDER/ecoli_html/toppic_proteoform_cutoff/data_js/proteoforms/ "123"


import sys
import pandas as pd
import json
import re
import os

def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts


def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass the directory and the protein accession, you can read the example by opening up the .py file and read the top command line")
    
    outputlist = []
    for filename in sorted(os.listdir(args[0]), key=numericalSort):
        try:
            with open(args[0] + filename) as file:
                file.readline()
                toppic = json.loads(file.read())
        except FileNotFoundError:
            print("Wrong directory, file not found!")
        
        if not toppic["compatible_proteoform"]["sequence_name"] == str(args[1]):
            continue

        dictlist = []
        prsm = toppic["compatible_proteoform"]["prsm"]
        if (not int(toppic["compatible_proteoform"]["prsm_number"]) == 1):
            prsm = prsm[0]
        peak_list = prsm["ms"]["peaks"]["peak"]
        for peak in peak_list:
            if "matched_ions" in peak:
                dictlist.append(peak)

        for dict in dictlist:
            matched_ion_num = int(dict["matched_ions_num"])
            if matched_ion_num > 1:
                for idx in range(0, matched_ion_num):
                    newdict = {}
                    newdict["Proteoform ID"] = filename
                    for key in dict.keys():
                        if not key == "matched_ions":
                            newdict[key] = dict[key]
                        else:
                            for smallkey in dict[key]["matched_ion"][idx].keys():
                                newdict[smallkey] = dict[key]["matched_ion"][idx][smallkey]
                    outputlist.append(newdict)
            else:
                newdict = {}
                newdict["Proteoform ID"] = filename
                for key in dict.keys():
                    if not key == "matched_ions":
                        newdict[key] = dict[key]
                    else:
                        for smallkey in dict[key]["matched_ion"].keys():
                            newdict[smallkey] = dict[key]["matched_ion"][smallkey]
                outputlist.append(newdict)
    
    df = pd.DataFrame(outputlist)

    df.to_csv("ion list.tsv", sep="\t")

if __name__ == "__main__":
    main()