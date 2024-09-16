import sys
import pandas as pd
import read_msalign

def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the directory")
    
    result_1 = pd.read_csv(args[0] + "resolved1_ms2_toppic_proteoform_single.tsv", delimiter="\t", skiprows=29)
    spec_list1 = read_msalign.read_spec_file(args[0] + "resolved1_ms2.msalign")
    spec_list2 = read_msalign.read_spec_file(args[0] + "resolved2_ms2.msalign")

    spectrumlist = []
    spec1_dict = {}
    for spec in spec_list1:
        spectrumlist.append(int(spec.header.spec_scan))
        spec1_dict[int(spec.header.spec_scan)] = spec

    prsmlist = result_1["Scan(s)"].tolist()
    
    count = 0
    outputlist = []
    for spec in spec_list2:
        scan = int(spec.header.spec_scan)
        if (scan not in spectrumlist) or (scan in prsmlist):
            count += 1
            outputlist.append(spec)

    print("We have " + str(count) + " spectra in second round after final filtering")

    read_msalign.write_spec_file(args[0] + "resolved2_ms2.msalign", outputlist)

if __name__ == "__main__":
    main()