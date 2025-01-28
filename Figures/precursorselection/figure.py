import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
  r1 = pd.read_csv("./20231215_ecoli_400ng_daniel_1_ms2_toppic_prsm_single.tsv", sep="\t")

  matcherror = pd.read_csv("./matcherror_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)

  shifterror = pd.read_csv("./shifterror_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)

  print(matcherror.shape[0], shifterror.shape[0])

  matchmerge = r1.merge(matcherror, on="Scan(s)", how="inner")
  shiftmerge = r1.merge(shifterror, on="Scan(s)", how="inner")

  print(matchmerge.shape[0], shiftmerge.shape[0])

  newmatchmerge = matchmerge[matchmerge["Protein accession_x"] == matchmerge["Protein accession_y"]]
  newshiftmerge = shiftmerge[shiftmerge["Protein accession_x"] == shiftmerge["Protein accession_y"]]

  print(newmatchmerge.shape[0], newshiftmerge.shape[0])

  newmatchmerge["matched peaks diff"] = newmatchmerge["#matched peaks_x"] - newmatchmerge["#matched peaks_y"]
  newmatchmerge["matched ions diff"] = newmatchmerge["#matched fragment ions_x"] - newmatchmerge["#matched fragment ions_y"]
  newmatchmerge["E-value diff"] = np.log10(newmatchmerge["E-value_y"]) - np.log10(newmatchmerge["E-value_x"])
  # print(newmatchmerge["matched peaks diff"].max())
  newshiftmerge["matched peaks diff"] = newshiftmerge["#matched peaks_x"] - newshiftmerge["#matched peaks_y"]
  newshiftmerge["matched ions diff"] = newshiftmerge["#matched fragment ions_x"] - newshiftmerge["#matched fragment ions_y"]
  newshiftmerge["E-value diff"] = np.log10(newshiftmerge["E-value_y"]) - np.log10(newshiftmerge["E-value_x"])

  # print(newshiftmerge["matched ions diff"].max())

  # bins = np.arange(-2, 26, 2)

  # fig, ax1 = plt.subplots(figsize=(10, 6))

  # sns.histplot(r1merge['matched ions diff'], bins=bins, kde=False, edgecolor='black', alpha=0.6, label="Distribution of ECOLI-MATCH vs ECOLI ERROR")
  # sns.histplot(shiftmerge['matched ions diff'], bins=bins, kde=False, edgecolor='black', alpha=0.6, label="Distribution of ECOLI-SHIFT vs ECOLI ERROR")

  # plt.xlabel("Change in # of Matched Theoretical Fragments")
  # plt.ylabel("# of PrSMs")
  # # plt.legend(loc='upper left', fontsize=12)
  plt.rcParams.update({'font.size': 17})

  # Define bins for the histogram
  # bins = np.arange(-2, 28, 2)

  bins_col1 = np.arange(-1.5, 45.5+1, 1)
  bins_col2 = np.arange(-3.5, 25.5+1, 1)
  bins_col3 = np.arange(-4, 16+1, 1)

  fig, ([ax1, ax2, ax3], [ax4, ax5, ax6]) = plt.subplots(2, 3, figsize=(20, 12))

  # Plot the histograms for the first column
  sns.histplot(newmatchmerge['matched peaks diff'], ax=ax1, bins=bins_col1)
  ax1.set_title("MATCH vs MATCH-ERROR")
  ax1.set_ylabel("PrSMs")
  ax1.set_xlabel("")

  sns.histplot(newshiftmerge['matched peaks diff'], ax=ax4, bins=bins_col1)
  ax4.set_title("SHIFT vs SHIFT-ERROR")
  ax4.set_xlabel("Decrease in matched experimental masses")
  ax4.set_ylabel("PrSMs")

  # Plot the histograms for the second column
  sns.histplot(newmatchmerge['matched ions diff'], ax=ax2, bins=bins_col2)
  ax2.set_title("MATCH vs MATCH-ERROR")
  ax2.set_ylabel("")
  ax2.set_xlabel("")

  sns.histplot(newshiftmerge['matched ions diff'], ax=ax5, bins=bins_col2)
  ax5.set_title("SHIFT vs SHIFT-ERROR")
  ax5.set_xlabel("Decrease in matched theoretical masses")
  ax5.set_ylabel("")

  # Plot the histograms for the third column
  sns.histplot(newmatchmerge['E-value diff'], ax=ax3, bins=bins_col3)
  ax3.set_title("MATCH vs MATCH-ERROR")
  ax3.set_ylabel("")
  ax3.set_xlabel("")

  sns.histplot(newshiftmerge['E-value diff'], ax=ax6, bins=bins_col3)
  ax6.set_title("SHIFT vs SHIFT-ERROR")
  ax6.set_xlabel("Increase in log E-value")
  ax6.set_ylabel("")

  plt.tight_layout()

  plt.savefig("PrecursorSelection.svg", format="svg", dpi=800)        
  

if __name__ == "__main__":
    main()
