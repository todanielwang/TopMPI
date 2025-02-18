import sys
import pandas as pd
import read_msalign

def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the directory")
    
    verifier = pd.read_csv("/home/daniel/Desktop/datafiles/RealData/ecoli/DDA/20231215_ecoli_400ng/20231215_ecoli_400ng_daniel_1_ms2_toppic_prsm_single.tsv", delimiter="\t")

    result_a = pd.read_csv(args[0] + "A_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)
    result_ab = pd.read_csv(args[0] + "AB_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)
    result_b = pd.read_csv(args[0] + "B_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)
    result_ba = pd.read_csv(args[0] + "BA_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)

    spec_list_a = read_msalign.read_spec_file(args[0] + "A_ms2.msalign")
    spec_list_ab = read_msalign.read_spec_file(args[0] + "AB_ms2.msalign")
    spec_list_b = read_msalign.read_spec_file(args[0] + "B_ms2.msalign")
    spec_list_ba = read_msalign.read_spec_file(args[0] + "BA_ms2.msalign")


    spec_dict_a = {}
    for spec in spec_list_a:
        spec_dict_a[str(spec.header.spec_scan)] = spec
    
    spec_dict_ab = {}
    for spec in spec_list_ab:
        spec_dict_ab[str(spec.header.spec_scan)] = spec
    
    spec_dict_b = {}
    for spec in spec_list_b:
        spec_dict_b[str(spec.header.spec_scan)] = spec

    spec_dict_ba = {}
    for spec in spec_list_ba:
        spec_dict_ba[str(spec.header.spec_scan)] = spec


    for index, row in result_a.iterrows():
        result_a.loc[index, "Scan(other)"] = int(spec_dict_a[str(row["Scan(s)"])].header.title) % 100000

    for index, row in result_ab.iterrows():
        result_ab.loc[index, "Scan(other)"] = int(spec_dict_ab[str(row["Scan(s)"])].header.title) % 100000

    for index, row in result_b.iterrows():
        result_b.loc[index, "Scan(other)"] = int(spec_dict_b[str(row["Scan(s)"])].header.title) % 100000

    for index, row in result_ba.iterrows():
        result_ba.loc[index, "Scan(other)"] = int(spec_dict_ba[str(row["Scan(s)"])].header.title) % 100000

    result_a["Scan_true"] = result_a["Scan(s)"] % 100000
    result_ab["Scan_true"] = result_ab["Scan(s)"] % 100000
    result_b["Scan_true"] = result_b["Scan(s)"] % 100000
    result_ba["Scan_true"] = result_ba["Scan(s)"] % 100000

    result_a = result_a.sort_values(by=["Scan_true", "Scan(other)"])
    result_ab = result_ab.sort_values(by=["Scan_true", "Scan(other)"])
    result_b = result_b.sort_values(by=["Scan_true", "Scan(other)"])
    result_ba = result_ba.sort_values(by=["Scan_true", "Scan(other)"])


    verifiedResult_a = pd.merge(result_a,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan_true', right_on="Scan(s)", how='left', suffixes=("", "_TrueA"))
    finalResult_a = pd.merge(verifiedResult_a,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan(other)', right_on="Scan(s)", how='left', suffixes=("", "_TrueB"))

    verifiedResult_ab = pd.merge(result_ab,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan_true', right_on="Scan(s)", how='left', suffixes=("", "_TrueA"))
    finalResult_ab = pd.merge(verifiedResult_ab,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan(other)', right_on="Scan(s)", how='left', suffixes=("", "_TrueB"))

    verifiedResult_b = pd.merge(result_b,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan_true', right_on="Scan(s)", how='left', suffixes=("", "_TrueB"))
    finalResult_b = pd.merge(verifiedResult_b,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan(other)', right_on="Scan(s)", how='left', suffixes=("", "_TrueA"))

    verifiedResult_ba = pd.merge(result_ba,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan_true', right_on="Scan(s)", how='left', suffixes=("", "_TrueB"))
    finalResult_ba = pd.merge(verifiedResult_ba,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan(other)', right_on="Scan(s)", how='left', suffixes=("", "_TrueA"))

    # finalResult["Verified_ProteoformA"] = (finalResult["Protein accession"] == finalResult["Protein accession_TrueA"]) & (abs(finalResult["Proteoform mass"] - finalResult["Proteoform mass_TrueA"]) < 1.5)
    # finalResult["Verified_ProteoformB"] = (finalResult["Protein accession"] == finalResult["Protein accession_TrueB"]) & (abs(finalResult["Proteoform mass"] - finalResult["Proteoform mass_TrueB"]) < 1.5)

    finalResult_a["Verified_ProteinA"] = (finalResult_a["Protein accession"] == finalResult_a["Protein accession_TrueA"])
    finalResult_a["Verified_ProteinB"] = (finalResult_a["Protein accession"] == finalResult_a["Protein accession_TrueB"])

    finalResult_ab["Verified_ProteinA"] = (finalResult_ab["Protein accession"] == finalResult_ab["Protein accession_TrueA"])
    finalResult_ab["Verified_ProteinB"] = (finalResult_ab["Protein accession"] == finalResult_ab["Protein accession_TrueB"])

    finalResult_b["Verified_ProteinA"] = (finalResult_b["Protein accession"] == finalResult_b["Protein accession_TrueA"])
    finalResult_b["Verified_ProteinB"] = (finalResult_b["Protein accession"] == finalResult_b["Protein accession_TrueB"])

    finalResult_ba["Verified_ProteinA"] = (finalResult_ba["Protein accession"] == finalResult_ba["Protein accession_TrueA"])
    finalResult_ba["Verified_ProteinB"] = (finalResult_ba["Protein accession"] == finalResult_ba["Protein accession_TrueB"])

    finalResult_a.drop(columns=["Scan_true"], inplace=True)
    finalResult_ab.drop(columns=["Scan_true"], inplace=True)
    finalResult_b.drop(columns=["Scan_true"], inplace=True)
    finalResult_ba.drop(columns=["Scan_true"], inplace=True)


    finalResult_a.to_csv(args[0] + "A_SimulationResult.csv", sep="\t")
    finalResult_ab.to_csv(args[0] + "AB_SimulationResult.csv", sep="\t")
    finalResult_b.to_csv(args[0] + "B_SimulationResult.csv", sep="\t")
    finalResult_ba.to_csv(args[0] + "BA_SimulationResult.csv", sep="\t")
        
    pairs = []
    for spec in spec_list_a:
        pairs.append(frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)))

    accessA = {}
    for pair in pairs:
        for spec in spec_list_a:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessA[pair] = spec.header.spec_scan

    accessAB = {}
    for pair in pairs:
        for spec in spec_list_ab:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessAB[pair] = spec.header.spec_scan

    accessB = {}
    for pair in pairs:
        for spec in spec_list_b:
            if frozenset((spec.header.spec_scan % 100000, int(spec.header.title) % 100000)) == pair:
                accessB[pair] = spec.header.spec_scan

    accessBA = {}
    for pair in pairs:
        for spec in spec_list_ba:
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
        elif "DECOY" == str(prsm.iloc[0].at["Protein accession"]).split("_")[0]:
            return "DECOY"
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
        # pair = eval(pair)
        output_dict["Pair"].append(pair)
        prsm_A = finalResult_a.loc[finalResult_a["Scan(s)"] == accessA[pair]]
        output_dict["A+B_1"].append(checkCondition(prsm_A))
        output_dict["A+B_1 peaks"].append(getValue(prsm_A, "#matched peaks"))
        output_dict["A+B_1 E-value"].append(getValue(prsm_A, "E-value"))
        prsm_AB = finalResult_ab.loc[finalResult_ab["Scan(s)"] == accessAB[pair]]
        output_dict["A+B_2"].append(checkCondition(prsm_AB))
        output_dict["A+B_2 peaks"].append(getValue(prsm_AB, "#matched peaks"))
        output_dict["A+B_2 E-value"].append(getValue(prsm_AB, "E-value"))
        prsm_B = finalResult_b.loc[finalResult_b["Scan(s)"] == accessB[pair]]
        output_dict["B+A_1"].append(checkCondition(prsm_B))
        output_dict["B+A_1 peaks"].append(getValue(prsm_B, "#matched peaks"))
        output_dict["B+A_1 E-value"].append(getValue(prsm_B, "E-value"))
        prsm_BA = finalResult_ba.loc[finalResult_ba["Scan(s)"] == accessBA[pair]]
        output_dict["B+A_2"].append(checkCondition(prsm_BA))
        output_dict["B+A_2 peaks"].append(getValue(prsm_BA, "#matched peaks"))
        output_dict["B+A_2 E-value"].append(getValue(prsm_BA, "E-value"))

    outputdf = pd.DataFrame(output_dict)

    outputdf.to_csv(args[0] + "Result.tsv", sep="\t")

if __name__ == "__main__":
    main()