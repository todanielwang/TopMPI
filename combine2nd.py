import read_msalign
import sys
import re
import pandas as pd
import os
import json
 
def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts



def main():
    args = sys.argv[1:]
    if (len(args) != 4):
        raise Exception(
            "Please pass in the ms.align file, the feature file, the start ratio, and the js file directory")
    
    spec_list = read_msalign.read_spec_file(args[0])

    mzwindow = 4

    startRatio = int(args[2])

    feature_file = pd.read_csv(args[1], sep="\t")

    for idx in range(len(spec_list) - 1, -1, -1):
        isolationWindow = (spec_list[idx].header.spec_scan - 1) % 21
        lowerbound = startRatio + (isolationWindow - 1) * mzwindow 
        upperbound = startRatio + (isolationWindow) * mzwindow
        rt = spec_list[idx].header.retention_time / 60

        query = feature_file[(feature_file['RepMz'] > lowerbound) & (feature_file['RepMz'] < upperbound) \
                             & (feature_file['MinElutionTime'] < rt) & (feature_file['MaxElutionTime'] > rt)]
    
        if (len(query) <= 1):
            del spec_list[idx]
            continue
        
        temp = query.drop([query['Abundance'].idxmax()])
        secMax = temp.loc[temp['Abundance'].idxmax()]

        spec_list[idx].header.ms_one_id = int(secMax['FeatureID'])
        spec_list[idx].header.ms_one_scan = int((spec_list[idx].header.spec_scan - 1) / 21)
        spec_list[idx].header.mono_mz = secMax['RepMz']
        spec_list[idx].header.charge = int(secMax['MinCharge'])
        spec_list[idx].header.mono_mass = secMax['MonoMass']
        spec_list[idx].header.inte = secMax['Abundance']

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

    
    print(count)
    read_msalign.write_spec_file(args[0], spec_list)   
    
    

if __name__ == "__main__":
    main()