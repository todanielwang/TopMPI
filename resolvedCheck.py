import sys
import pandas as pd
import read_msalign

def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the directory")
    
    verifier = pd.read_csv("/home/daniel/Desktop/datafiles/RealData/ecoli/DDA/20231215_ecoli_400ng/20231215_ecoli_400ng_daniel_1_ms2_toppic_prsm_single.tsv", delimiter="\t")

    result_a = pd.read_csv(args[0] + "resolved1_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)
    result_b = pd.read_csv(args[0] + "resolved2_ms2_toppic_prsm_single.tsv", delimiter="\t", skiprows=29)

    spec_list_a = read_msalign.read_spec_file(args[0] + "resolved1_ms2.msalign")
    spec_list_b = read_msalign.read_spec_file(args[0] + "resolved2_ms2.msalign")


    spec_dict_a = {}
    for spec in spec_list_a:
        spec_dict_a[str(spec.header.spec_scan)] = spec
    
    
    spec_dict_b = {}
    for spec in spec_list_b:
        spec_dict_b[str(spec.header.spec_scan)] = spec


    for index, row in result_a.iterrows():
        result_a.loc[index, "Scan(other)"] = int(spec_dict_a[str(row["Scan(s)"])].header.title) % 100000

    for index, row in result_b.iterrows():
        result_b.loc[index, "Scan(other)"] = int(spec_dict_b[str(row["Scan(s)"])].header.title) % 100000

    result_a["Scan_true"] = result_a["Scan(s)"] % 100000
    result_b["Scan_true"] = result_b["Scan(s)"] % 100000

    result_a = result_a.sort_values(by=["Scan_true", "Scan(other)"])
    result_b = result_b.sort_values(by=["Scan_true", "Scan(other)"])


    verifiedResult_a = pd.merge(result_a,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan_true', right_on="Scan(s)", how='left', suffixes=("", "_TrueCorrect"))
    finalResult_a = pd.merge(verifiedResult_a,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan(other)', right_on="Scan(s)", how='left', suffixes=("", "_TrueWrong"))

    verifiedResult_b = pd.merge(result_b,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan_true', right_on="Scan(s)", how='left', suffixes=("", "_TrueCorrect"))
    finalResult_b = pd.merge(verifiedResult_b,verifier[['Scan(s)','Proteoform', "Protein accession", "Proteoform mass", "#matched peaks", "#unexpected modifications"]], left_on='Scan(other)', right_on="Scan(s)", how='left', suffixes=("", "_TrueWrong"))

    # finalResult["Verified_ProteoformA"] = (finalResult["Protein accession"] == finalResult["Protein accession_TrueA"]) & (abs(finalResult["Proteoform mass"] - finalResult["Proteoform mass_TrueA"]) < 1.5)
    # finalResult["Verified_ProteoformB"] = (finalResult["Protein accession"] == finalResult["Protein accession_TrueB"]) & (abs(finalResult["Proteoform mass"] - finalResult["Proteoform mass_TrueB"]) < 1.5)

    finalResult_a["Verified_ProteinCorrect"] = (finalResult_a["Protein accession"] == finalResult_a["Protein accession_TrueCorrect"])
    finalResult_a["Verified_ProteinWrong"] = (finalResult_a["Protein accession"] == finalResult_a["Protein accession_TrueWrong"])

    finalResult_b["Verified_ProteinCorrect"] = (finalResult_b["Protein accession"] == finalResult_b["Protein accession_TrueCorrect"])
    finalResult_b["Verified_ProteinWrong"] = (finalResult_b["Protein accession"] == finalResult_b["Protein accession_TrueWrong"])

    finalResult_a.drop(columns=["Scan_true"], inplace=True)
    finalResult_b.drop(columns=["Scan_true"], inplace=True)


    finalResult_a.to_csv(args[0] + "Resolved1_Result.csv", sep="\t")
    finalResult_b.to_csv(args[0] + "Resolved2_Result.csv", sep="\t")

if __name__ == "__main__":
    main()