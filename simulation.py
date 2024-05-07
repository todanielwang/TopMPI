import sys
import csv
import read_msalign
import copy

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the prsm file and then the msalign file")

    with open(args[0]) as f:
        prsms = [{k: v for k, v in row.items()}
            for row in csv.DictReader(f, skipinitialspace=True, delimiter="\t")]

    pairs = set()

    for idx, prsm in list(enumerate(prsms))[:-1]:
        for idxother, prsmother in list(enumerate(prsms))[idx:]:
            if abs(float(prsmother["m/z"]) - float(prsm["m/z"])) < 1.5:
                if (((int(prsmother["Charge"]) != int(prsm["Charge"]))) & (prsmother["Protein accession"] != prsm["Protein accession"])):
                    if (int(prsm["#matched peaks"]) >= int(prsmother["#matched peaks"])):
                        pairs.add(tuple([prsm["Scan(s)"], prsmother["Scan(s)"]]))
                    else:
                        pairs.add(tuple([prsmother["Scan(s)"], prsm["Scan(s)"]]))

    spec_list = read_msalign.read_spec_file(args[1])

    spec_dict = {}
    for spec in spec_list:
        if (len(spec.header.pre_inte_list) == 1 or float(spec.header.pre_inte_list[0]) > 0.8 * sum(map(float, spec.header.pre_inte_list))):
            spec_dict[str(spec.header.spec_scan)] = spec
    temp = []
    import pandas as pd
    df = pd.DataFrame(prsms)
    for spec in spec_dict.keys():
        if spec in df["Scan(s)"].tolist():
            temp.append(1)
    print(len(temp))
    
    output_list = []
    for pair in pairs:
        pair = list(pair)
        if (pair[0] not in spec_dict or pair[1] not in spec_dict):
            continue
        mainSpec = copy.deepcopy(spec_dict[pair[0]])
        sideSpec = spec_dict[pair[1]]
        mainSpec.header.title = str(sideSpec.header.spec_scan)
        mainSpec.header.pre_mz_list = [mainSpec.header.pre_mz_list[0], sideSpec.header.pre_mz_list[0]]
        mainSpec.header.pre_charge_list = [mainSpec.header.pre_charge_list[0], sideSpec.header.pre_charge_list[0]]
        mainSpec.header.pre_mass_list = [mainSpec.header.pre_mass_list[0], sideSpec.header.pre_mass_list[0]]
        mainSpec.header.pre_inte_list = [mainSpec.header.pre_inte_list[0], sideSpec.header.pre_inte_list[0]]
        mainSpec.header.pre_id_list = [mainSpec.header.pre_id_list[0], sideSpec.header.pre_id_list[0]]
        mainSpec.peak_list = mainSpec.peak_list + sideSpec.peak_list
        output_list.append(mainSpec)

    sortedlist = read_msalign.sortScans(output_list)

    print(len(sortedlist))

    read_msalign.write_spec_file("A_ms2.msalign", sortedlist)

    for spec in sortedlist:
        title = int(spec.header.title) % 100000
        spec.header.title = str(int(spec.header.spec_scan) % 100000)
        spec.header.spec_scan = title

        pre_mz_list = spec.header.pre_mz_list
        pre_charge_list = spec.header.pre_charge_list
        pre_mass_list = spec.header.pre_mass_list
        pre_inte_list = spec.header.pre_inte_list
        pre_id_list = spec.header.pre_id_list

        pre_mz_list[0], pre_mz_list[1] = pre_mz_list[1], pre_mz_list[0]
        pre_charge_list[0], pre_charge_list[1] = pre_charge_list[1], pre_charge_list[0]
        pre_mass_list[0], pre_mass_list[1] = pre_mass_list[1], pre_mass_list[0]
        pre_inte_list[0], pre_inte_list[1] = pre_inte_list[1], pre_inte_list[0]
        pre_id_list[0], pre_id_list[1] = pre_id_list[1], pre_id_list[0]

        spec.header.pre_mz_list = pre_mz_list
        spec.header.pre_charge_list = pre_charge_list
        spec.header.pre_mass_list = pre_mass_list
        spec.header.pre_inte_list = pre_inte_list
        spec.header.pre_id_list = pre_id_list

    othersortedlist = read_msalign.sortScans(sortedlist)

    read_msalign.write_spec_file("B_ms2.msalign", othersortedlist)


        

if __name__ == "__main__":
    main()