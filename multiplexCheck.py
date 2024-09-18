import sys
import pandas as pd
import read_msalign

def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the directory")
    
    result_a = pd.read_csv(args[0] + "A_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)
    result_ab = pd.read_csv(args[0] + "AB_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)
    result_b = pd.read_csv(args[0] + "B_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)
    result_ba = pd.read_csv(args[0] + "BA_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)

    a_scans = set(result_a["Scan(s)"].tolist())
    ab_scans = set(result_ab["Scan(s)"].tolist())
    b_scans = set(result_b["Scan(s)"].tolist())
    ba_scans = set(result_ba["Scan(s)"].tolist())

    scanlist = sorted(list(a_scans.union(ab_scans).union(b_scans).union(ba_scans)))

    # spec_list_a = read_msalign.read_spec_file(args[0] + "A_ms2.msalign")
    # spec_list_ab = read_msalign.read_spec_file(args[0] + "AB_ms2.msalign")
    # spec_list_b = read_msalign.read_spec_file(args[0] + "B_ms2.msalign")
    # spec_list_ba = read_msalign.read_spec_file(args[0] + "BA_ms2.msalign")


    # spec_dict_a = {}
    # for spec in spec_list_a:
    #     spec_dict_a[str(spec.header.spec_scan)] = spec
    
    # spec_dict_ab = {}
    # for spec in spec_list_ab:
    #     spec_dict_ab[str(spec.header.spec_scan)] = spec
    
    # spec_dict_b = {}
    # for spec in spec_list_b:
    #     spec_dict_b[str(spec.header.spec_scan)] = spec

    # spec_dict_ba = {}
    # for spec in spec_list_ba:
    #     spec_dict_ba[str(spec.header.spec_scan)] = spec

    output_dict = {"Scan": [], 
                    "A+B_1": [], 
                    "A+B_1_mass": [], 
                    "A+B_1 peaks": [], 
                    "A+B_1 E-value": [], 
                    "A+B_2": [], 
                    "A+B_2_mass": [], 
                    "A+B_2 peaks": [], 
                    "A+B_2 E-value": [], 
                    "B+A_1": [], 
                    "B+A_1_mass": [], 
                    "B+A_1 peaks": [], 
                    "B+A_1 E-value": [], 
                    "B+A_2": [], 
                    "B+A_2_mass": [], 
                    "B+A_2 peaks": [], 
                    "B+A_2 E-value": []}
    
            
    def getValue(prsm, key):
        if prsm.empty:
            if key == "Protein accession":
                return "-"
            elif key == "#matched peaks" or key == "Proteoform mass":
                return 0
            elif key == "E-value":
                return 1
        return prsm.iloc[0].at[key]
        

    for scan in scanlist:
        output_dict["Scan"].append(scan)
        prsm_A = result_a.loc[result_a["Scan(s)"] == int(scan)]
        output_dict["A+B_1"].append(getValue(prsm_A, "Protein accession"))
        output_dict["A+B_1_mass"].append(getValue(prsm_A, "Proteoform mass"))
        output_dict["A+B_1 peaks"].append(getValue(prsm_A, "#matched peaks"))
        output_dict["A+B_1 E-value"].append(getValue(prsm_A, "E-value"))
        prsm_AB = result_ab.loc[result_ab["Scan(s)"] == int(scan)]
        output_dict["A+B_2"].append(getValue(prsm_AB, "Protein accession"))
        output_dict["A+B_2_mass"].append(getValue(prsm_AB, "Proteoform mass"))
        output_dict["A+B_2 peaks"].append(getValue(prsm_AB, "#matched peaks"))
        output_dict["A+B_2 E-value"].append(getValue(prsm_AB, "E-value"))
        prsm_B = result_b.loc[result_b["Scan(s)"] == int(scan)]
        output_dict["B+A_1"].append(getValue(prsm_B, "Protein accession"))
        output_dict["B+A_1_mass"].append(getValue(prsm_B, "Proteoform mass"))
        output_dict["B+A_1 peaks"].append(getValue(prsm_B, "#matched peaks"))
        output_dict["B+A_1 E-value"].append(getValue(prsm_B, "E-value"))
        prsm_BA = result_ba.loc[result_ba["Scan(s)"] == int(scan)]
        output_dict["B+A_2"].append(getValue(prsm_BA, "Protein accession"))
        output_dict["B+A_2_mass"].append(getValue(prsm_BA, "Proteoform mass"))
        output_dict["B+A_2 peaks"].append(getValue(prsm_BA, "#matched peaks"))
        output_dict["B+A_2 E-value"].append(getValue(prsm_BA, "E-value"))

    outputdf = pd.DataFrame(output_dict)

    outputdf.to_csv(args[0] + "Evaluecutoff10000/Result.tsv", sep="\t")

if __name__ == "__main__":
    main()