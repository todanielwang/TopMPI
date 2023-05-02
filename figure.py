import read_msalign
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mplt
import numpy as np
import random

def main():
  args = sys.argv[1:]
  if (len(args) != 3):
    raise Exception("Please pass in the ms.align file, the feature file, and the start ratio")
    
  spec_list = read_msalign.read_spec_file(args[0])

  mzwindow = 4

  startRatio = int(args[2])

  feature_file = pd.read_csv(args[1], sep="\t")

  mz = []
  intensity = []
  for spec in spec_list:
     if (spec.header.spec_scan == 6808):
        for peak in spec.peak_list:
           mz.append(peak.mass / peak.charge)
           intensity.append(peak.intensity)
  
  intensity = np.array(intensity) / max(intensity)
  Amz = mz[1::2]
  Aintensity = intensity[1::2]
  Bmz = mz[::2]
  Bintensity = intensity[::2]

  fig, (ax1, ax2, ax3) = plt.subplots(3)
  fig.set_figheight(10)
  fig.set_figwidth(10)
  ax1.yaxis.set_major_formatter(mplt.ticker.PercentFormatter(xmax=1.0))
  ax2.yaxis.set_major_formatter(mplt.ticker.PercentFormatter(xmax=1.0))
  ax3.yaxis.set_major_formatter(mplt.ticker.PercentFormatter(xmax=1.0))
  #plt.autoscale()
  ax1.bar(Amz, Aintensity, width=2, color="Red")
  ax1.bar(Bmz, Bintensity, width=2, color="Blue")
  ax1.set(title="Mixture Spectrum", xlabel="m/z", ylabel="Intensity(%)")
  ax1.label_outer()

  ax2.bar(Amz, Aintensity, width=2, color="Red")
  ax2.set(ylim=[0, 1], xlabel="m/z", ylabel="Intensity(%)")
  ax2.label_outer()
  
  ax3.bar(Bmz, Bintensity, width=2, color="Blue")
  ax3.set(xlabel="m/z", ylabel="Intensity(%)")
  ax3.label_outer()

  plt.savefig("temp.png")
  plt.close()

  """
  filtered = feature_file[(feature_file['RepMz'] > 992) & (feature_file['RepMz'] < 996) & (feature_file["MinElutionTime"] > 27) & (feature_file["MaxElutionTime"] < 32)]
  #(feature_file['RepMz'] > 992) & (feature_file['RepMz'] < 996) & (feature_file["MinElutionTime"] > 27) & (feature_file["MaxElutionTime"] < 32)
  featuremz = filtered["RepMz"]
  featureminrt = filtered["MinElutionTime"]
  featuremaxrt = filtered["MaxElutionTime"]
  featurecmap = np.log2(filtered["Abundance"])

  normalizedcmap = (featurecmap - np.min(featurecmap)) / (np.max(featurecmap) - np.min(featurecmap))
  
  f = plt.figure()
  f.set_figwidth(20)
  f.set_figheight(20)
  hsv = mplt.colormaps["cool"]
  colors = [hsv(x) for x in normalizedcmap]
  plt.xticks(fontsize=18)
  plt.yticks(fontsize=18)
  plt.xlabel('m/z Ratio', fontsize=18)
  plt.ylabel('Retention Time', fontsize=18)
  plt.vlines(x = featuremz, ymin = featureminrt, ymax = featuremaxrt, colors = colors)
  plt.autoscale()
  plt.savefig("temp.png")
  plt.close()
  """
        
  

if __name__ == "__main__":
    main()