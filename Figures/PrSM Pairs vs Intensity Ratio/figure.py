import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import read_msalign

def main():
  r1 = pd.read_csv("./A_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)
  r2 = pd.read_csv("./AB_ms2_toppic_prsm_single.tsv", sep="\t", skiprows=29)

  print(r1.shape[0], r2.shape[0])

  pairs = r1.merge(r2, on="Scan(s)", how="inner")

  print(pairs.shape[0])

  scans = pairs["Scan(s)"].tolist()

  df = pd.DataFrame({
      'scan': pd.Series(dtype='int'),
      'multiplexed': pd.Series(dtype='bool'),
      'ratio': pd.Series(dtype='float')
  })

  spectra = read_msalign.read_spec_file("./A_ms2.msalign")

  totalscans = [int(spec.header.spec_scan) for spec in spectra]

  multiplexed = [True if int(spec.header.spec_scan) in scans else False for spec in spectra]

  ratio = [float(spec.header.pre_inte_list[1]) / float(spec.header.pre_inte_list[0]) if not (spec.header.pre_inte_list[0] == '' or float(spec.header.pre_inte_list[0]) == float(0) or len(spec.header.pre_inte_list) == 1 or float(spec.header.pre_inte_list[1]) == float(0)) else 0 for spec in spectra]

  df["scan"] = totalscans
  df["multiplexed"] = multiplexed
  df["ratio"] = ratio

  # Define the bins for the ratio column
  bins = np.linspace(0, 1, 11)  # 10 bins from 0 to 1 (inclusive)
  labels = [f"{round(bins[i], 1)} - {round(bins[i+1], 1)}" for i in range(len(bins) - 1)]

  # Add a 'group' column to categorize each row into the appropriate bin
  df['group'] = pd.cut(df['ratio'], bins=bins, labels=labels, include_lowest=True)

  outputdf = df[df["multiplexed"] == True]

  print(outputdf['ratio'].min())
  # outputdf = pd.read_csv("./data.tsv", sep="\t")

  plt.rcParams.update({'font.size': 17})

  fig, ax1 = plt.subplots(figsize=(10, 6))

  bins = np.arange(0.0, 1.1, 0.1)

  sns.histplot(outputdf['ratio'], bins=bins, kde=False, edgecolor='black', label="Distribution of # of PrSM Pairs")

  ax1.set_xlabel('Intensity ratio')
  ax1.set_ylabel('Scans')

  plt.tight_layout()
  # ax1.set_title('Histogram of Ratio with Proportion of Multiplexed=True')

  # bin_centers = np.arange(0.25, 1.05, 0.1)

  # ax2 = ax1.twinx()
  # ax2.plot(bin_centers, proportions_percentage[2:], color='red', label="Distribution of % of PrSM Pairs")
  # ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1f}%'))
  # ax2.set_ylabel("% of Scans where PrSM Pairs Reported")

  # lines_1, labels_1 = ax1.get_legend_handles_labels()
  # lines_2, labels_2 = ax2.get_legend_handles_labels()
  # ax1.legend(lines_1, labels_1, loc='upper left', fontsize=12)


  plt.savefig("PrSM Pairs vs Intensity Ratio.svg", format="svg", dpi=800)
        
  

if __name__ == "__main__":
    main()
