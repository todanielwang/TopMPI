import matplotlib.pyplot as plt
import numpy as np
import read_msalign
import sys

def main():
    args = sys.argv[1:]
    if (len(args) != 6):
        raise Exception(
            "Please pass in msalign files")
    
    peaks = []

    for i in range(len(args)):
        spec_list = read_msalign.read_spec_file(args[i])

        for idx in range(len(spec_list) - 1, -1, -1):
            peaks.append(len(spec_list[idx].peak_list))
            
    hist, bins = np.histogram(peaks, bins = 50)

    # Plot the resulting histogram
    center = (bins[:-1]+bins[1:])/2
    width = 0.7*(bins[1]-bins[0])
    plt.bar(center, hist, align = 'center', width = width)
    plt.yscale("log")
    plt.xlabel('Fragment Masses Lengths')
    plt.ylabel('Frequency')

    plt.savefig("masslength.png")


if __name__ == "__main__":
    main()