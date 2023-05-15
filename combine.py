#Map the scans with first feature

import read_msalign
import sys
import pandas as pd

def main():
    args = sys.argv[1:]
    if (len(args) != 3):
        raise Exception(
            "Please pass in the ms.align file, the feature file, and the start ratio")
    
    
    
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
        
        if (query.empty):
            del spec_list[idx]
            continue
        
        maxMatch = query.loc[query['Abundance'].idxmax()]
        spec_list[idx].header.ms_one_id = int((spec_list[idx].header.spec_scan - 1) / 21)
        spec_list[idx].header.ms_one_scan = int(maxMatch['FeatureID']) 
        spec_list[idx].header.mono_mz = maxMatch['RepMz']
        spec_list[idx].header.charge = int(maxMatch['MinCharge'])
        spec_list[idx].header.mono_mass = maxMatch['MonoMass']
        spec_list[idx].header.inte = maxMatch['Abundance']


    read_msalign.write_spec_file(args[0], spec_list)

if __name__ == "__main__":
    main()