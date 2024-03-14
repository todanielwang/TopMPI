import sys
import pandas as pd

def main():
    args = sys.argv[1:]

    r2 = pd.read_csv(args[0], delimiter="\t")

    r1 = pd.read_csv(args[1], delimiter="\t")

    scan_list = r1["Scan(s)"].to_list()

    string = str(args[0]).split(".")

    r2[r2["Scan(s)"].isin(scan_list)].to_csv(string[0] + "_r1." + string[1], sep="\t")
    r2[~r2["Scan(s)"].isin(scan_list)].to_csv(string[0] + "_r2." + string[1], sep="\t")

        

if __name__ == "__main__":
    main()