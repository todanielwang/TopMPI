import sys
import pandas as pd
import read_msalign

def main():
    args = sys.argv[1:]

    AB = pd.read_csv(args[0], delimiter="\t")

    BA = pd.read_csv(args[1], delimiter="\t")

    a_list = AB[AB["Verified_ProteinB_r1"] & AB["Verified_ProteinB"]]["Scan(s)"].to_list()

    b_list = BA[BA["Verified_ProteinA_r1"] & BA["Verified_ProteinA"]]["Scan(s)"].to_list()

    print(len(a_list), len(b_list))

    spec_list_A = read_msalign.read_spec_file(args[2])

    spec_list_B = read_msalign.read_spec_file(args[3])

    spec_dict_A = {}
    for spec in spec_list_A:
        spec_dict_A[str(spec.header.spec_scan)] = spec

    spec_dict_B = {}
    for spec in spec_list_B:
        spec_dict_B[str(spec.header.spec_scan)] = spec

    A = pd.read_csv(args[4], delimiter="\t")
    B = pd.read_csv(args[5], delimiter="\t")

    scan_list_A = A[A["Verified_ProteinB"]]["Scan(s)"].to_list()
    scan_list_B = B[B["Verified_ProteinA"]]["Scan(s)"].to_list()

    accessA = {}
    for spec in spec_list_A:
        accessA[frozenset((int(spec.header.spec_scan) % 100000, int(spec.header.title) % 100000))] = int(spec.header.spec_scan)

    accessB = {}
    for spec in spec_list_B:
        accessB[frozenset((int(spec.header.spec_scan) % 100000, int(spec.header.title) % 100000))] = int(spec.header.spec_scan)

    output_a = []
    for scan in scan_list_A:
        if B.loc[B["Scan(s)"] == accessB[frozenset((int(scan) % 100000, int(spec_dict_A[str(scan)].header.title)))]]["Verified_ProteinB"].bool():
            output_a.append(scan)

    output_b = []
    for scan in scan_list_B:
        print(A.loc[A["Scan(s)"] == accessA[frozenset((int(scan) % 100000, int(spec_dict_B[str(scan)].header.title)))]]["Verified_ProteinA"], scan)
            # output_b.append(scan)    

    def Union(lst1, lst2):
        final_list = list(set(lst1) | set(lst2))
        return final_list
    
    flipalist = Union(a_list, output_a)
    flipblist = Union(b_list, output_b)

    for spec in spec_list_A:
        if spec.header.spec_scan in flipalist:
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

    for spec in spec_list_B:
        if spec.header.spec_scan in flipblist:
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

    outputAlist = read_msalign.sortScans(spec_list_A)
    outputBlist = read_msalign.sortScans(spec_list_B)

    read_msalign.write_spec_file("A.msalign", spec_list_A)
    read_msalign.write_spec_file("B.msalign", spec_list_B)


if __name__ == "__main__":
    main()