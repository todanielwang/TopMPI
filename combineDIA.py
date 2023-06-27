#Map the scans with first feature

import read_msalign
import sys
import pandas as pd
import math

def main():
    args = sys.argv[1:]
    if (len(args) != 3):
        raise Exception(
            "Please pass in the ms.align file, the feature file, and the start ratio")
    
    
    
    spec_list = read_msalign.read_spec_file(args[0])

    mzwindow = 1

    startRatio = int(args[2])

    feature_file = pd.read_csv(args[1], sep=",", index_col=False, converters={'XIC': lambda string: list(map(float, string[:].split(';')))})

    for idx in range(len(spec_list) - 1, -1, -1):
        ms1round = int((spec_list[idx].header.spec_scan - 1) / 21)
        isolationWindow = (spec_list[idx].header.spec_scan - 1) % 21
        lowerbound = startRatio + (isolationWindow - 1) * mzwindow 
        upperbound = startRatio + (isolationWindow) * mzwindow
        rt = spec_list[idx].header.retention_time / 60

        query = feature_file[(feature_file.XIC.apply(lambda x: not math.isnan(x[ms1round]) and x[ms1round] > 1)) & (feature_file['MonoMz'] > lowerbound) & (feature_file['MonoMz'] < upperbound) \
                             & (feature_file['rtLo'] < rt) & (feature_file['rtHi'] > rt)]
        
        if (query.empty):
            del spec_list[idx]
            continue

        maxMatch = query.loc[query['XIC'].apply(lambda x : x[ms1round]).idxmax()]
        spec_list[idx].header.ms_one_id = int(maxMatch['ID'])
        spec_list[idx].header.mono_mz = maxMatch['MonoMz']
        spec_list[idx].header.charge = int(maxMatch['Charge'])
        spec_list[idx].header.mono_mass = maxMatch['Mass']
        spec_list[idx].header.inte = maxMatch['XIC'][ms1round]


    read_msalign.write_spec_file(args[0], spec_list)

if __name__ == "__main__":
    main()