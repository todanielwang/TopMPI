import sys
import pandas as pd
import numpy as np


def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the result.csv file")

    inputdf = pd.read_csv(args[0], delimiter="\t", index_col=0)

    mask = (inputdf['A+B_1'] != "-") & (inputdf['A+B_1'] == inputdf["A+B_2"]) & (inputdf["A+B_1 peaks"] < inputdf["B+A_1 peaks"])
    inputdf.loc[mask, "A+B_1"] = "-"
    inputdf.loc[mask, ["A+B_1 peaks", "A+B_1 E-value"]] = 0, 1

    mask = (inputdf['A+B_1'] != "-") & (inputdf['A+B_1'] == inputdf["A+B_2"]) & (inputdf["A+B_1 peaks"] >= inputdf["B+A_1 peaks"])
    inputdf.loc[mask, "A+B_2"] = "-"
    inputdf.loc[mask, ["A+B_2 peaks", "A+B_2 E-value"]] = 0, 1

    mask = (inputdf['B+A_1'] != "-") & (inputdf['B+A_1'] == inputdf["B+A_2"]) & (inputdf["A+B_1 peaks"] > inputdf["B+A_1 peaks"])
    inputdf.loc[mask, "B+A_1"] = "-"
    inputdf.loc[mask, ["B+A_1 peaks", "B+A_1 E-value"]] = 0, 1

    mask = (inputdf['B+A_1'] != "-") & (inputdf['B+A_1'] == inputdf["B+A_2"]) & (inputdf["A+B_1 peaks"] <= inputdf["B+A_1 peaks"])
    inputdf.loc[mask, "B+A_2"] = "-"
    inputdf.loc[mask, ["B+A_2 peaks", "B+A_2 E-value"]] = 0, 1

    inputdf["choice"] = ""
    mask = ((inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] > inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) | (inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] == inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) & (inputdf["A+B_1 E-value"] < inputdf["B+A_1 E-value"]))
    inputdf.loc[mask, "choice"] = "A"

    mask = ((inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] < inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) | (inputdf["A+B_1 peaks"] + inputdf["A+B_2 peaks"] == inputdf["B+A_1 peaks"] + inputdf["B+A_2 peaks"]) & (inputdf["A+B_1 E-value"] > inputdf["B+A_1 E-value"]))
    inputdf.loc[mask, "choice"] = "B"

    inputdf.to_csv(args[0].split(".")[0] + "_resolved." + args[0].split(".")[1], sep="\t")
        

if __name__ == "__main__":
    main()