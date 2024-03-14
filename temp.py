import read_msalign
import sys
import numpy as np
import pandas as pd

def main():
    args = sys.argv[1:]
    
    a_spec_list = read_msalign.read_spec_file(args[0])
    b_spec_list = read_msalign.read_spec_file(args[1])
    ab_spec_list = read_msalign.read_spec_file(args[2])
    ba_spec_list = read_msalign.read_spec_file(args[3])

    A = pd.read_csv(args[4], delimiter="\t")
    B = pd.read_csv(args[5], delimiter="\t")
    AB = pd.read_csv(args[6], delimiter="\t")
    BA = pd.read_csv(args[7], delimiter="\t")

    
    pairs = []
    for spec in a_spec_list:
        pairs.append(frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)))

    accessA = {}
    for pair in pairs:
        for spec in a_spec_list:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessA[pair] = spec.header.spec_scan

    accessB = {}
    for pair in pairs:
        for spec in b_spec_list:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessB[pair] = spec.header.spec_scan

    accessAB = {}
    for pair in pairs:
        for spec in ab_spec_list:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessAB[pair] = spec.header.spec_scan

    accessBA = {}
    for pair in pairs:
        for spec in ba_spec_list:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessBA[pair] = spec.header.spec_scan

    output_dict = {"Pair": [], 
                   "A+B_1": [], 
                   "A+B_1 peaks": [], 
                   "A+B_1 E-value": [], 
                   "A+B_2": [], 
                   "A+B_2 peaks": [], 
                   "A+B_2 E-value": [], 
                   "B+A_1": [], 
                   "B+A_1 peaks": [],
                   "B+A_1 E-value": [], 
                   "B+A_2": [], 
                   "B+A_2 peaks": [], 
                   "B+A_2 E-value": []}
    
    def checkCondition(prsm):
        if prsm.empty:
            return "-"
        elif prsm.iloc[0].at["Verified_ProteinA"]:
            return "A"
        elif prsm.iloc[0].at["Verified_ProteinB"]:
            return "B"
        else:
            return "C"
        
    def getValue(prsm, key):
        if prsm.empty:
            if key == "#matched peaks":
                return 0
            else:
                return 1
        return prsm.iloc[0].at[key]
    

    for pair in pairs:
        output_dict["Pair"].append(pair)
        prsm_A = A.loc[A["Scan(s)"] == accessA[pair]]
        output_dict["A+B_1"].append(checkCondition(prsm_A))
        output_dict["A+B_1 peaks"].append(getValue(prsm_A, "#matched peaks"))
        output_dict["A+B_1 E-value"].append(getValue(prsm_A, "E-value"))
        prsm_AB = AB.loc[AB["Scan(s)"] == accessAB[pair]]
        output_dict["A+B_2"].append(checkCondition(prsm_AB))
        output_dict["A+B_2 peaks"].append(getValue(prsm_AB, "#matched peaks"))
        output_dict["A+B_2 E-value"].append(getValue(prsm_AB, "E-value"))
        prsm_B = B.loc[B["Scan(s)"] == accessB[pair]]
        output_dict["B+A_1"].append(checkCondition(prsm_B))
        output_dict["B+A_1 peaks"].append(getValue(prsm_B, "#matched peaks"))
        output_dict["B+A_1 E-value"].append(getValue(prsm_B, "E-value"))
        prsm_BA = BA.loc[BA["Scan(s)"] == accessBA[pair]]
        output_dict["B+A_2"].append(checkCondition(prsm_BA))
        output_dict["B+A_2 peaks"].append(getValue(prsm_BA, "#matched peaks"))
        output_dict["B+A_2 E-value"].append(getValue(prsm_BA, "E-value"))

    outputdf = pd.DataFrame(output_dict)

    outputdf.to_csv("Result.tsv", sep="\t")


    

    # a2b = set()
    # b2a = set()
    # condition = "#matched peaks"
    # for pair in pairs:
    #     a_scan = A.loc[A["Scan(s)"] == accessA[pair]]
    #     b_scan = B.loc[B["Scan(s)"] == accessB[pair]]
    #     ab_scan = AB.loc[AB["Scan(s)"] == accessAB[pair]]
    #     ba_scan = BA.loc[BA["Scan(s)"] == accessBA[pair]]
    #     if (len(a_scan) > 0 and len(b_scan) > 0 and a_scan.iloc[0]["Protein accession"] == b_scan.iloc[0]["Protein accession"]):
    #         if a_scan.iloc[0][condition] > b_scan.iloc[0][condition]:
    #             b2a.add(pair)
    #         elif a_scan.iloc[0][condition] < b_scan.iloc[0][condition]:
    #             a2b.add(pair)
    #     if (len(a_scan) > 0 and len(ab_scan) > 0 and a_scan.iloc[0]["Protein accession"] == ab_scan.iloc[0]["Protein accession"]):
    #         if a_scan.iloc[0][condition] > ab_scan.iloc[0][condition]:
    #             b2a.add(pair)
    #         elif a_scan.iloc[0][condition] < ab_scan.iloc[0][condition]:
    #             a2b.add(pair)
    #     if (len(b_scan) > 0 and len(ba_scan) > 0 and b_scan.iloc[0]["Protein accession"] == ba_scan.iloc[0]["Protein accession"]):
    #         if b_scan.iloc[0][condition] > ba_scan.iloc[0][condition]:
    #             a2b.add(pair)
    #         elif b_scan.iloc[0][condition] < ba_scan.iloc[0][condition]:
    #             b2a.add(pair)

    # print(len(a2b), len(b2a))
    # a2blist = []
    # b2alist = []
    # for pair in a2b:
    #     a2blist.append(accessA[pair])
    # for pair in b2a:
    #     b2alist.append(accessB[pair])

    # a2ba = A[A["Scan(s)"].isin(a2blist) & A["Verified_ProteinA"]]
    # a2bb = A[A["Scan(s)"].isin(a2blist) & A["Verified_ProteinB"]]
    # print(len(a2ba) / len(a2b), len(a2bb) / len(a2b))
    # b2aa = B[B["Scan(s)"].isin(b2alist) & B["Verified_ProteinA"]]
    # b2ab = B[B["Scan(s)"].isin(b2alist) & B["Verified_ProteinB"]]
    # print(len(b2aa) / len(b2a), len(b2ab) / len(b2a))







    # spec_dict_a = {}
    # for spec in a_spec_list:
    #         spec_dict_a[spec.header.spec_scan] = spec

    # spec_dict_b = {}
    # for spec in b_spec_list:
    #         spec_dict_b[spec.header.spec_scan] = spec

    # for pair in a2b:
    #     spec = spec_dict_a[accessA[pair]]
    #     title = int(spec.header.title) % 100000
    #     spec.header.title = str(spec.header.spec_scan % 100000)
    #     spec.header.spec_scan = title
    #     spec.header.spec_id = title
    #     pre_mz_list = spec.header.pre_mz_list
    #     pre_charge_list = spec.header.pre_charge_list
    #     pre_mass_list = spec.header.pre_mass_list
    #     pre_inte_list = spec.header.pre_inte_list
    #     pre_id_list = spec.header.pre_id_list

    #     pre_mz_list[0], pre_mz_list[1] = pre_mz_list[1], pre_mz_list[0]
    #     pre_charge_list[0], pre_charge_list[1] = pre_charge_list[1], pre_charge_list[0]
    #     pre_mass_list[0], pre_mass_list[1] = pre_mass_list[1], pre_mass_list[0]
    #     pre_inte_list[0], pre_inte_list[1] = pre_inte_list[1], pre_inte_list[0]
    #     pre_id_list[0], pre_id_list[1] = pre_id_list[1], pre_id_list[0]

    #     spec.header.pre_mz_list = pre_mz_list
    #     spec.header.pre_charge_list = pre_charge_list
    #     spec.header.pre_mass_list = pre_mass_list
    #     spec.header.pre_inte_list = pre_inte_list
    #     spec.header.pre_id_list = pre_id_list

    # sortedA = read_msalign.sortScans(a_spec_list)

    # read_msalign.write_spec_file(args[0], sortedA)

    # for pair in b2a:
    #     spec = spec_dict_b[accessB[pair]]
    #     title = int(spec.header.title) % 100000
    #     spec.header.title = str(spec.header.spec_scan % 100000)
    #     spec.header.spec_scan = title
    #     spec.header.spec_id = title
    #     pre_mz_list = spec.header.pre_mz_list
    #     pre_charge_list = spec.header.pre_charge_list
    #     pre_mass_list = spec.header.pre_mass_list
    #     pre_inte_list = spec.header.pre_inte_list
    #     pre_id_list = spec.header.pre_id_list

    #     pre_mz_list[0], pre_mz_list[1] = pre_mz_list[1], pre_mz_list[0]
    #     pre_charge_list[0], pre_charge_list[1] = pre_charge_list[1], pre_charge_list[0]
    #     pre_mass_list[0], pre_mass_list[1] = pre_mass_list[1], pre_mass_list[0]
    #     pre_inte_list[0], pre_inte_list[1] = pre_inte_list[1], pre_inte_list[0]
    #     pre_id_list[0], pre_id_list[1] = pre_id_list[1], pre_id_list[0]

    #     spec.header.pre_mz_list = pre_mz_list
    #     spec.header.pre_charge_list = pre_charge_list
    #     spec.header.pre_mass_list = pre_mass_list
    #     spec.header.pre_inte_list = pre_inte_list
    #     spec.header.pre_id_list = pre_id_list

    # sortedB = read_msalign.sortScans(b_spec_list)

    # read_msalign.write_spec_file(args[1], sortedB)
    

if __name__ == "__main__":
    main()