import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


def main():
  df = pd.read_csv("./data.tsv", sep="\t")

  plt.rcParams.update({'font.size': 25})

  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 10))

  # sns.lineplot(data=df, x="% of Noise Added", y="# of PrSMs g0", ax=ax0, linewidth=3)
  # sns.lineplot(data=df, x="% of Noise Added", y="# of PrSMs g1", ax=ax0, linewidth=3)
  # sns.lineplot(data=df, x="% of Noise Added", y="# of PrSMs g2", ax=ax0, linewidth=3)
  # sns.lineplot(data=df, x="% of Noise Added", y="# of PrSMs g3", ax=ax0, linewidth=3)

  # ax0.set_ylabel('# of PrSMs with Correct Protein Identification')
  # ax0.legend()

  sns.lineplot(data=df, x="% of Noise Added", y="Accuracy g0", ax=ax1, linewidth=3)
  sns.lineplot(data=df, x="% of Noise Added", y="Accuracy g1", ax=ax1, linewidth=3)
  sns.lineplot(data=df, x="% of Noise Added", y="Accuracy g2", ax=ax1, linewidth=3)
  sns.lineplot(data=df, x="% of Noise Added", y="Accuracy g3", ax=ax1, linewidth=3)

  ax1.set_ylabel('Protein Identification Accuracy(%)')
  # ax1.legend()

  ax1.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))

  sns.lineplot(data=df, x="% of Noise Added", y="False g0", ax=ax2, label="Group 1", linewidth=3)
  sns.lineplot(data=df, x="% of Noise Added", y="False g1", ax=ax2, label="Group 2", linewidth=3)
  sns.lineplot(data=df, x="% of Noise Added", y="False g2", ax=ax2, label="Group 3", linewidth=3)
  sns.lineplot(data=df, x="% of Noise Added", y="False g3", ax=ax2, label="Group 4", linewidth=3)

  ax2.set_ylabel('Noise Protein Identification Rate (%)')
  # ax2.legend()

  ax2.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))

  # fig.suptitle('Group 3')
  sns.move_legend(ax2, loc='upper right', fancybox=True, framealpha=1, bbox_to_anchor=(1, 1.125), ncol = 4)
  plt.savefig("noise.svg", dpi=800, format="svg")
        
  

if __name__ == "__main__":
    main()
