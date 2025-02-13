import os
import shutil
import subprocess
import argparse
import shlex
import glob
import switchPrecursor
import checkAndRemovePeaks
import merge
import combine
import removeTemp

def copy_and_rename(src, dst):
    """Copy and rename files if they exist."""
    if os.path.isfile(src):
        shutil.copy(src, dst)
        print(f"Copied and renamed: {os.path.basename(dst)} from {os.path.basename(src)}")
    else:
        print(f"File not found: {src}")


import argparse

def main():
    parser = argparse.ArgumentParser(description="Python driver for TopMPI.")
    
    parser.add_argument("TopPIC", help="A TopPIC executable")
    parser.add_argument("database", help="A protein database file in the FASTA format")
    parser.add_argument("input_files", nargs="+", help="A mass spectrum data file in the msalign format (wildcard '*' can be used).")
    parser.add_argument("-a", "--activation", choices=["CID", "HCD", "ETD", "UVPD", "FILE"], default="FILE", help="Set the fragmentation method(s) of MS/MS spectra.")
    parser.add_argument("-f", "--fixed-mod", help="Set fixed modifications: C57, C58, or a file containing modifications.")
    parser.add_argument("-n", "--n-terminal-form", default="NONE,M_ACETYLATION,NME,NME_ACETYLATION", help="Set N-terminal forms of proteins.")
    parser.add_argument("-s", "--num-shift", type=int, choices=[0, 1, 2], default=1, help="Maximum number of unexpected mass shifts.")
    parser.add_argument("-m", "--min-shift", type=float, default=-500, help="Minimum value for unexpected mass shifts in Dalton.")
    parser.add_argument("-M", "--max-shift", type=float, default=500, help="Maximum value for unexpected mass shifts in Dalton.")
    parser.add_argument("-S", "--variable-ptm-num", type=int, default=3, help="Maximum number of variable PTM sites.")
    parser.add_argument("-b", "--variable-ptm-file-name", help="Variable PTM file.")
    parser.add_argument("-d", "--decoy", action="store_true", help="Use a shuffled decoy protein database for FDR estimation.")
    parser.add_argument("-e", "--mass-error-tolerance", type=int, default=10, help="Set error tolerance for precursor and fragment masses in ppm.")
    parser.add_argument("-p", "--proteoform-error-tolerance", type=float, default=1.2, help="Error tolerance for PrSM clusters in Dalton.")
    parser.add_argument("-t", "--spectrum-cutoff-type", choices=["EVALUE", "FDR"], default="EVALUE", help="Spectrum-level cutoff type.")
    parser.add_argument("-v", "--spectrum-cutoff-value", type=float, default=0.01, help="Spectrum-level cutoff value.")
    parser.add_argument("-T", "--proteoform-cutoff-type", choices=["EVALUE", "FDR"], default="EVALUE", help="Proteoform-level cutoff type.")
    parser.add_argument("-V", "--proteoform-cutoff-value", type=float, default=0.01, help="Proteoform-level cutoff value.")
    parser.add_argument("-A", "--approximate-spectra", action="store_true", help="Use approximate spectra for protein filtering.")
    parser.add_argument("-l", "--lookup-table", action="store_true", help="Use a lookup table for computing E-values.")
    parser.add_argument("-B", "--local-ptm-file-name", help="File containing common PTMs for proteoform characterization.")
    parser.add_argument("-H", "--miscore-threshold", type=float, default=0.15, help="Set the MIScore threshold for PTM characterization.")
    parser.add_argument("-u", "--thread-number", type=int, default=1, help="Number of threads for computation.")
    parser.add_argument("-r", "--num-combined-spectra", type=int, default=1, help="Number of combined spectra.")
    parser.add_argument("-c", "--combined-file-name", help="Output file for combined identifications.")
    parser.add_argument("-x", "--no-topfd-feature", action="store_true", help="Specify that there are no TopFD feature files.")
    parser.add_argument("-k", "--keep-temp-files", action="store_true", help="Keep intermediate files.")
    parser.add_argument("-K", "--keep-decoy-ids", action="store_true", help="Keep decoy identifications.")
    parser.add_argument("-g", "--skip-html-folder", action="store_true", help="Skip HTML folder generation.")
    parser.add_argument("--alpha", type=float, default=0.2, help="The intensity ratio between the first and second precursor required for a spectrum to be treated as multiplexed.")
    parser.add_argument("--beta", type=float, default=0.9, help="The percentage of shared matched experimental peaks required for identifications of two precursors to be treated as inconsistent.")
    parser.add_argument("--delta", type=int, default=5, help="The offset to calculate the number of normalized matched fragment masses (NMFMs) based on number of unknown mass shifts.")
    parser.add_argument("--gamma", type=int, default=4, help="The number of normalized matched fragment masses (NMFMs) difference required to switch to the second precursor.")
    
    args = parser.parse_args()

    if not args.decoy and args.spectrum_cutoff_type == "FDR":
        print("Spectrum-level cutoff type FDR error! FDR cutoff cannot be used when no decoy database is used! Please add argument '-d' in the command.")
        exit(1)
    elif not args.decoy and args.proteoform_cutoff_type == "FDR":
        print("Proteoform-level cutoff type FDR error! FDR cutoff cannot be used when no decoy database is used! Please add argument '-d' in the command.")
        exit(1)

    input_files = []
    for file_pattern in args.input_files:
        matched_files = glob.glob(file_pattern)  # Expands '*' to matching files
        if matched_files:
            input_files.extend(matched_files)
        else:
            input_files.append(file_pattern)  # Keep original if no match found

    for input_file in input_files:
        base_dir = os.path.dirname(input_file)
        common_prefix = os.path.basename(input_file).rsplit('_', 1)[0]
        new_sub_dir = os.path.join(base_dir, f"{common_prefix}_TopMPI")
        
        os.makedirs(new_sub_dir, exist_ok=True)

        copy_and_rename(os.path.join(base_dir, f"{common_prefix}_ms2.msalign"), os.path.join(new_sub_dir, "First_ms2.msalign"))

        if not args.no_topfd_feature:
            for ext in ["feature.xml", "ms1.feature", "ms1.msalign", "ms2.feature"]:
                copy_and_rename(os.path.join(base_dir, f"{common_prefix}_{ext}"), os.path.join(new_sub_dir, f"First_{ext}"))

        switchPrecursor.main([os.path.join(new_sub_dir, "First_ms2.msalign")])

        os.rename(os.path.join(new_sub_dir, "First_ms2_modified.msalign"), os.path.join(new_sub_dir, "Second_ms2.msalign"))

        if not args.no_topfd_feature:
            for ext in ["feature.xml", "ms1.feature", "ms1.msalign", "ms2.feature"]:
                copy_and_rename(os.path.join(new_sub_dir, f"First_{ext}"), os.path.join(new_sub_dir, f"Second_{ext}"))

        #split flags
        TopPIC_args = {k: v for k, v in vars(args).items() if k not in ["alpha", "beta", "delta", "gamma", "spectrum_cutoff_type", "spectrum_cutoff_value", "proteoform_cutoff_type", "proteoform_cutoff_value", "TopPIC", "database", "input_files", "combined_file_name"]}
        # extra_args = {k: v for k, v in vars(args).items() if k in ["alpha", "beta", "delta", "gamma", "spectrum_cutoff_type", "spectrum_cutoff_value", "proteoform_cutoff_type", "proteoform_cutoff_value"]}
        
        TopPIC_flags = []
        for k, v in TopPIC_args.items():
            formatted_key = k.replace("_", "-")
            if isinstance(v, bool) and v:
                TopPIC_flags.append(f"--{formatted_key}")
            elif v not in [None, False]:
                TopPIC_flags.extend(shlex.split(f"--{formatted_key} {v}"))

        if args.decoy and "--keep-decoy-ids" not in TopPIC_flags:
            TopPIC_flags.append("--keep-decoy-ids")

        # extra_flags = [f"--{k}" if isinstance(v, bool) and v else f"--{k} {v}" for k, v in extra_args.items() if v is not None]

        # print(TopPIC_flags)
        subprocess.run([args.TopPIC, args.database, os.path.join(new_sub_dir, "First_ms2.msalign"), "-v", "100000", "-V", "100000"] + TopPIC_flags, check=True)
        subprocess.run([args.TopPIC, args.database, os.path.join(new_sub_dir, "Second_ms2.msalign"), "-v", "100000", "-V", "100000"] + TopPIC_flags, check=True)

        
        # print("\nExecuting extra parameters program with:")
        # print("extra.exe", " ".join(extra_flags))
        # subprocess.run(["extra.exe"] + extra_flags, check=True)

        checkAndRemovePeaks.main([new_sub_dir, "-a", str(args.alpha), "-b", str(args.beta), "-d", str(args.delta), "-g", str(args.gamma)])

        os.rename(os.path.join(new_sub_dir, "Primary_ms2_modified.msalign"), os.path.join(new_sub_dir, "Primary_ms2.msalign"))
        os.rename(os.path.join(new_sub_dir, "Secondary_ms2_modified.msalign"), os.path.join(new_sub_dir, "Secondary_ms2.msalign"))

        if not args.no_topfd_feature:
            for ext in ["feature.xml", "ms1.feature", "ms1.msalign", "ms2.feature"]:
                copy_and_rename(os.path.join(new_sub_dir, f"First_{ext}"), os.path.join(new_sub_dir, f"Secondary_{ext}"))

        subprocess.run([args.TopPIC, args.database, os.path.join(new_sub_dir, "Secondary_ms2.msalign"), "-v", "100000", "-V", "100000"] + TopPIC_flags, check=True)

        # Merge results
        filterbyFeature = str(not args.no_topfd_feature)

        merge.main([new_sub_dir, filterbyFeature, "-t", args.spectrum_cutoff_type, "-v", str(args.spectrum_cutoff_value), "-T", args.proteoform_cutoff_type, "-V", str(args.proteoform_cutoff_value)])


    if args.combined_file_name:
        combine.main([args.combined_file_name, input_files, "-t", args.spectrum_cutoff_type, "-v", str(args.spectrum_cutoff_value), "-T", args.proteoform_cutoff_type, "-V", str(args.proteoform_cutoff_value)])

    for input_file in input_files:
        if not args.keep_temp_files:
            base_dir = os.path.dirname(input_file)
            common_prefix = os.path.basename(input_file).rsplit('_', 1)[0]
            new_sub_dir = os.path.join(base_dir, f"{common_prefix}_TopMPI")
            removeTemp.main([new_sub_dir])


if __name__ == "__main__":
    main()