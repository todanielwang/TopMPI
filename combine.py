import pandas as pd
import util
import os
def combine(combined_name, input_files, spectrumcutofftype, spectrumcutoffvalue, proteoformcutofftype, proteoformcutoffvalue):
    outputdir = combined_name + "_TopMPI"
    
    os.makedirs(outputdir, exist_ok=True)

    totalprsm_list = []
    
    print("Combining files")
    for index, input_file in enumerate(input_files):
        base_dir = os.path.dirname(input_file)
        common_prefix = os.path.basename(input_file).rsplit('_', 1)[0]
        new_sub_dir = os.path.join(base_dir, f"{common_prefix}_TopMPI")

        primary = util.read_tsv(os.path.join(new_sub_dir, "Primary_ms2_temp_prsm.tsv"))
        primary["Spectrum ID"] += index * 10000000
        primary["Feature ID"] += index * 10000000

        totalprsm_list.append(primary)

        secondary = util.read_tsv(os.path.join(new_sub_dir, "Secondary_ms2_toppic_prsm_single.tsv"))

        secondary["Spectrum ID"] += index * 10000000
        secondary["Feature ID"] += index * 10000000
        totalprsm_list.append(secondary)

    totalprsm = pd.concat(totalprsm_list, ignore_index=True)

    totalproteoform = util.getProteoforms(totalprsm)

    if spectrumcutofftype == "EVALUE":
        combined_prsm = totalprsm[~totalprsm['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
        combined_prsm_output = combined_prsm[combined_prsm["E-value"] < spectrumcutoffvalue]
        combined_prsm_output.to_csv(os.path.join(outputdir, "TotalPrSM.tsv"), sep="\t", index=False)
    else:
        totalprsm["Spectrum-level Q-value"] = util.calculate_q_values(totalprsm)
        combined_prsm = totalprsm[~totalprsm['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
        combined_prsm_output = combined_prsm[combined_prsm["Spectrum-level Q-value"] < spectrumcutoffvalue]
        combined_prsm_output.to_csv(os.path.join(outputdir, "TotalPrSM.tsv"), sep="\t", index=False)

    if proteoformcutofftype == "EVALUE":
        combined_proteoform = totalproteoform[~totalproteoform['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
        combined_proteoform_output = combined_proteoform[combined_proteoform["E-value"] < proteoformcutoffvalue]
        combined_proteoform_output.to_csv(os.path.join(outputdir, "TotalProteoform.tsv"), sep="\t", index=False)
    else:
        totalproteoform["Proteoform-level Q-value"] = util.calculate_q_values(totalproteoform)
        combined_proteoform = totalproteoform[~totalproteoform['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
        combined_proteoform_output = combined_proteoform[combined_proteoform["Proteoform-level Q-value"] < proteoformcutoffvalue]
        combined_proteoform_output.to_csv(os.path.join(outputdir, "TotalProteoform.tsv"), sep="\t", index=False)

    print("Under spectrum-level cutoff type " + spectrumcutofftype + " and cutoff value " + str(spectrumcutoffvalue))
    print("Outputing " + str(combined_prsm_output.shape[0]) + " combined total PrSMs")

    print("Under proteoform-level cutoff type " + proteoformcutofftype + " and cutoff value " + str(proteoformcutoffvalue))
    print("Outputing " + str(combined_proteoform_output.shape[0]) + " combined total proteoforms")