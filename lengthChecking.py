import read_msalign
import sys
import numpy as np
import pandas as pd

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file and the prsm file directory")
    
    prsm_df = pd.read_csv(args[1], delimiter="\t")

    scan_list = prsm_df["Scan(s)"].tolist()

    spec_list = read_msalign.read_spec_file(args[0])

    fragmentlengths = []
    fragmentmasses = []
    for spec in spec_list:
        if spec.header.spec_scan in scan_list:
            fragmentlengths.append(len(spec.peak_list))
            for peak in spec.peak_list:
                fragmentmasses.append(peak.mass)

    print(np.average(fragmentlengths), np.average(fragmentmasses))


if __name__ == "__main__":
    main()