import pandas as pd
import argparse
import util
import os
import util

def main(args_list=None):
    # Create the argument parser
    parser = argparse.ArgumentParser(description="The merging the rounds and some postprocessing of TopMPI, users should not see this")

    # Add the positional argument
    parser.add_argument("directory", help="The directory to the TopMPI folder")
    parser.add_argument("filterbyFeature", type=util.str_to_bool, help="Proteoform filter by feature or not")
    parser.add_argument("keepDecoys", type=util.str_to_bool, help="Keep decoy ids or not")

    # Add the optional flags
    parser.add_argument("-t", "--spectrum-cutoff-type", choices=["EVALUE", "FDR"], default="EVALUE", help="Spectrum-level cutoff type.")
    parser.add_argument("-v", "--spectrum-cutoff-value", type=float, default=0.01, help="Spectrum-level cutoff value.")
    parser.add_argument("-T", "--proteoform-cutoff-type", choices=["EVALUE", "FDR"], default="EVALUE", help="Proteoform-level cutoff type.")
    parser.add_argument("-V", "--proteoform-cutoff-value", type=float, default=0.01, help="Proteoform-level cutoff value.")
    parser.add_argument("-p", "--proteoform-error-tolerance", type=float, default=1.2, help="Error tolerance for PrSM clusters in Dalton.")

    # Parse the arguments
    args = parser.parse_args(args_list)

    primary_prsm_full = util.read_tsv(os.path.join(args.directory, "Primary_ms2_temp_prsm.tsv"))

    primary_proteoform_full = util.getProteoforms(primary_prsm_full, threshold=args.proteoform_error_tolerance, filterbyfeature=args.filterbyFeature)

    secondary_prsm_full = util.read_tsv(os.path.join(args.directory, "Secondary_ms2_toppic_prsm_single.tsv"))

    secondary_prsm_full.to_csv(os.path.join(args.directory, "Secondary_ms2_temp_prsm.tsv"), sep="\t", index=False)

    secondary_proteoform_full = util.getProteoforms(secondary_prsm_full, threshold=args.proteoform_error_tolerance, filterbyfeature=args.filterbyFeature)

    combined_prsm_full = pd.concat([primary_prsm_full, secondary_prsm_full]).reset_index(drop=True)

    combined_proteoform_full = util.getProteoforms(combined_prsm_full, threshold=args.proteoform_error_tolerance, filterbyfeature=args.filterbyFeature)

    #Removes all files in the specified directory that contain both 'Secondary' and 'toppic' in the filename.
    for file in os.listdir(args.directory):
        file_path = os.path.join(args.directory, file)
        if os.path.isfile(file_path) and "Secondary" in file and "toppic" in file:
            # print(f"Removing: {file_path}")
            os.remove(file_path)
    
    if args.spectrum_cutoff_type == "EVALUE":
        if not args.keepDecoys:
            primary_prsm_full = primary_prsm_full[~primary_prsm_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            secondary_prsm_full = secondary_prsm_full[~secondary_prsm_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            combined_prsm_full = combined_prsm_full[~combined_prsm_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
        
        primary_prsm_output = primary_prsm_full[primary_prsm_full["E-value"] < args.spectrum_cutoff_value].reset_index(drop=True)
        secondary_prsm_output = secondary_prsm_full[secondary_prsm_full["E-value"] < args.spectrum_cutoff_value]
        combined_prsm_output = combined_prsm_full[combined_prsm_full["E-value"] < args.spectrum_cutoff_value]

        primary_prsm_output.to_csv(os.path.join(args.directory, "Primary_ms2_toppic_prsm_single.tsv"), sep="\t", index=False)
        secondary_prsm_output.to_csv(os.path.join(args.directory, "Secondary_ms2_toppic_prsm_single.tsv"), sep="\t", index=False)
        combined_prsm_output.to_csv(os.path.join(args.directory, "TotalPrSM.tsv"), sep="\t", index=False)

    else:
        primary_prsm_full["Spectrum-level Q-value"] = util.calculate_q_values(primary_prsm_full)
        secondary_prsm_full["Spectrum-level Q-value"] = util.calculate_q_values(secondary_prsm_full)
        # combined_prsm_full["Spectrum-level Q-value"] = util.calculate_q_values(combined_prsm_full)

        if not args.keepDecoys:
            primary_prsm_full = primary_prsm_full[~primary_prsm_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            secondary_prsm_full = secondary_prsm_full[~secondary_prsm_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            # combined_prsm_full = combined_prsm_full[~combined_prsm_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)

        primary_prsm_output = primary_prsm_full[primary_prsm_full["Spectrum-level Q-value"] < args.spectrum_cutoff_value]
        secondary_prsm_output = secondary_prsm_full[secondary_prsm_full["Spectrum-level Q-value"] < args.spectrum_cutoff_value]
        # combined_prsm_output = combined_prsm_full[combined_prsm_full["Spectrum-level Q-value"] < args.spectrum_cutoff_value]
        combined_prsm_output = pd.concat([primary_prsm_output, secondary_prsm_output]).reset_index(drop=True)

        primary_prsm_output.to_csv(os.path.join(args.directory, "Primary_ms2_toppic_prsm_single.tsv"), sep="\t", index=False)
        secondary_prsm_output.to_csv(os.path.join(args.directory, "Secondary_ms2_toppic_prsm_single.tsv"), sep="\t", index=False)
        combined_prsm_output.to_csv(os.path.join(args.directory, "TotalPrSM.tsv"), sep="\t", index=False)

    if args.proteoform_cutoff_type == "EVALUE":
        if not args.keepDecoys:
            primary_proteoform_full = primary_proteoform_full[~primary_proteoform_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            secondary_proteoform_full = secondary_proteoform_full[~secondary_proteoform_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            combined_proteoform_full = combined_proteoform_full[~combined_proteoform_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)

        primary_proteoform_output = primary_proteoform_full[primary_proteoform_full["E-value"] < args.proteoform_cutoff_value]
        secondary_proteoform_output = secondary_proteoform_full[secondary_proteoform_full["E-value"] < args.proteoform_cutoff_value]
        combined_proteoform_output = combined_proteoform_full[combined_proteoform_full["E-value"] < args.proteoform_cutoff_value]

        primary_proteoform_output.to_csv(os.path.join(args.directory, "Primary_ms2_toppic_proteoform_single.tsv"), sep="\t", index=False)
        secondary_proteoform_output.to_csv(os.path.join(args.directory, "Secondary_ms2_toppic_proteoform_single.tsv"), sep="\t", index=False)
        combined_proteoform_output.to_csv(os.path.join(args.directory, "TotalProteoform.tsv"), sep="\t", index=False)

    else:
        primary_proteoform_full["Proteoform-level Q-value"] = util.calculate_q_values(primary_proteoform_full)
        secondary_proteoform_full["Proteoform-level Q-value"] = util.calculate_q_values(secondary_proteoform_full)
        combined_proteoform_full["Proteoform-level Q-value"] = util.calculate_q_values(combined_proteoform_full)

        if not args.keepDecoys:
            primary_proteoform_full = primary_proteoform_full[~primary_proteoform_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            secondary_proteoform_full = secondary_proteoform_full[~secondary_proteoform_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)
            combined_proteoform_full = combined_proteoform_full[~combined_proteoform_full['Protein accession'].str.contains('DECOY')].reset_index(drop=True)

        primary_proteoform_output = primary_proteoform_full[primary_proteoform_full["Proteoform-level Q-value"] < args.proteoform_cutoff_value]
        secondary_proteoform_output = secondary_proteoform_full[secondary_proteoform_full["Proteoform-level Q-value"] < args.proteoform_cutoff_value]
        combined_proteoform_output = combined_proteoform_full[combined_proteoform_full["Proteoform-level Q-value"] < args.proteoform_cutoff_value]

        primary_proteoform_output.to_csv(os.path.join(args.directory, "Primary_ms2_toppic_proteoform_single.tsv"), sep="\t", index=False)
        secondary_proteoform_output.to_csv(os.path.join(args.directory, "Secondary_ms2_toppic_proteoform_single.tsv"), sep="\t", index=False)
        combined_proteoform_output.to_csv(os.path.join(args.directory, "TotalProteoform.tsv"), sep="\t", index=False)


    print("Under spectrum-level cutoff type " + args.spectrum_cutoff_type + " and cutoff value " + str(args.spectrum_cutoff_value))
    print("Outputing " + str(primary_prsm_output.shape[0]) + " primary PrSMs")
    print("Outputing " + str(secondary_prsm_output.shape[0]) + " secondary PrSMs")
    print("Outputing " + str(combined_prsm_output.shape[0]) + " total PrSMs")

    print("Under proteoform-level cutoff type " + args.proteoform_cutoff_type + " and cutoff value " + str(args.proteoform_cutoff_value))
    print("Outputing " + str(primary_proteoform_output.shape[0]) + " primary proteoforms")
    print("Outputing " + str(secondary_proteoform_output.shape[0]) + " secondary proteoforms")
    print("Outputing " + str(combined_proteoform_output.shape[0]) + " total proteoforms")



if __name__ == "__main__":
    main()