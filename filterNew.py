#Map the scans with second feature

import read_msalign
import sys
import re
import pandas as pd
import os
import math
import json
 
def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts

def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))


def main():
    args = sys.argv[1:]
    if (len(args) != 4):
        raise Exception(
            "Please pass in the ms.align file, the feature file, the start ratio, and the js file directory")
    
    spec_list = read_msalign.read_spec_file(args[0])

    mzwindow = 1

    startRatio = int(args[2])

    feature_file = pd.read_csv(args[1], sep=",", index_col=False, converters={'XIC': lambda string: list(map(float, string[:].split(';')))})

    for idx in range(len(spec_list) - 1, -1, -1):
        ms_one_id = int((spec_list[idx].header.spec_scan - 1) / 21)
        isolationWindow = (spec_list[idx].header.spec_scan - 1) % 21
        lowerbound = startRatio + (isolationWindow - 1) * mzwindow 
        upperbound = startRatio + (isolationWindow) * mzwindow
        rt = spec_list[idx].header.retention_time / 60

        query = feature_file[(feature_file.XIC.apply(lambda x: not math.isnan(x[ms_one_id]) and x[ms_one_id] > 1)) & (feature_file['MonoMz'] > lowerbound) & (feature_file['MonoMz'] < upperbound) \
                             & (feature_file['rtLo'] < rt) & (feature_file['rtHi'] > rt)]
    
        if (len(query) <= 1):
            del spec_list[idx]
            continue
        
        
        query = query.drop([query['XIC'].apply(lambda x : x[ms_one_id]).idxmax()])

        secMax = query.loc[query['XIC'].apply(lambda x : x[ms_one_id]).idxmax()]

        condition = False
        while (overlap(lowerbound, upperbound, secMax["mzLo"], secMax["mzHi"]) / (upperbound - lowerbound) < 0.2):
            query = query.drop([query['XIC'].apply(lambda x : x[ms_one_id]).idxmax()])
            if (query.empty):
                condition = True
                break
            secMax = query.loc[query['XIC'].apply(lambda x : x[ms_one_id]).idxmax()]

        if (condition):
            continue

        spec_list[idx].header.ms_one_id = ms_one_id
        spec_list[idx].header.ms_one_scan = int(secMax['ID']) 
        spec_list[idx].header.mono_mz = secMax['MonoMz']
        spec_list[idx].header.charge = int(secMax['Charge'])
        spec_list[idx].header.mono_mass = secMax['Mass']
        spec_list[idx].header.inte = secMax['XIC'][ms_one_id]



    #peak removal
    
    curr_spec = 0
    count = 0
    for x in sorted(os.listdir(args[3]), key=numericalSort):
        if x.endswith(".js"):
            with open(args[3] + x) as file:
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