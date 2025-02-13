import pandas as pd
import argparse
import util
import os

def main(args_list=None):
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Removing the temp files, users should not see this")

    # Add the positional argument
    parser.add_argument("directory", help="The directory to the TopMPI folder")

    # Parse the arguments
    args = parser.parse_args(args_list)

    os.remove(os.path.join(args.directory, "Primary_ms2_temp_prsm.tsv"))
    os.remove(os.path.join(args.directory, "Secondary_ms2_temp_prsm.tsv"))



if __name__ == "__main__":
    main()