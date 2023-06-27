#Map the scans with first feature

import read_msalign
import sys
import pandas as pd
import math
import pymzml

def main():
    args = sys.argv[1:]
    if (len(args) != 3):
        raise Exception(
            "Please pass in the ms.align file, the feature file, and the mzml file")
    
    mzml = pymzml.run.Reader(
        args[2] ,
        MS1_Precision = 5e-6 ,
        MSn_Precision = 20e-6
    )

    mzdict = {}
    for spec in mzml:
        if spec.ms_level == 2:
            id = spec.ID
            targetmz = spec.__getitem__("MS:1000827")
            loweroffset = spec.__getitem__("MS:1000828")
            upperoffset = spec.__getitem__("MS:1000829")
            mzdict[id] = (targetmz, loweroffset, upperoffset)
    
    spec_list = read_msalign.read_spec_file(args[0])

    feature_file = pd.read_csv(args[1], sep=",", index_col=False, converters={'XIC': lambda string: list(map(float, string[:].split(';')))})

    for idx in range(len(spec_list) - 1, -1, -1):
        scan = spec_list[idx].header.spec_scan
        ms1round = int((spec_list[idx].header.ms_one_scan - 1) / 7)
        lowerbound = mzdict[scan][0] - mzdict[scan][1]
        upperbound = mzdict[scan][0] + mzdict[scan][2]
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