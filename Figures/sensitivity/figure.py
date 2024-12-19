import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def main():
  scatter_data = pd.read_csv("./data.tsv", sep="\t")

  sns.set(style='whitegrid')

  # Plotting the scatter plot
  plt.figure(figsize=(12, 8))

  # Print the count of each combination
  print("Combination Counts:\n", scatter_data['combination'].value_counts())

  # Adding transparency for clarity
  plt.figure(figsize=(12, 8))
  sns.scatterplot(data=scatter_data[scatter_data['combination'] == 'Both identifications are reported'], 
                  x='x', y='y', label='Both')

  sns.scatterplot(data=scatter_data[scatter_data['combination'] == 'Only one identification is reported'], 
                  x='x', y='y', label='One')

  sns.scatterplot(data=scatter_data[scatter_data['combination'] == 'No identifications are reported'], 
                  x='x', y='y', label='None')

  # Customize plot
  plt.xlabel('Logged E-value of First PrSM', fontsize=16)
  plt.ylabel('Logged E-value of Second PrSM', fontsize=16)

  plt.legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.06), ncol=3, frameon=False)

  plt.axvline(x=-2, color='black', linestyle='--', linewidth=1)
  plt.axhline(y=-2, color='black', linestyle='--', linewidth=1)

  plt.xticks(fontsize=12)
  plt.yticks(fontsize=12)

  # plt.suptitle('Reported Identifications under 0.01 E-value Cutoff', fontsize=20)

  plt.savefig("sensitivity.svg", dpi=800, format="svg")
        
  

if __name__ == "__main__":
    main()
