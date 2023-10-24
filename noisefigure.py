import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mtick


def main():
    xpoints = np.array([0, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200])
    ypoints1 = np.array([1891, 1758, 1709, 1671, 1643, 1609, 1566, 1538, 1504, 1457, 1380, 1330, 1259, 1231])
    ypoints2 = np.array([1891, 1750, 1697, 1659, 1633, 1591, 1550, 1524, 1502, 1440, 1366, 1313, 1240, 1229])

    ypoints1 = [x / 1891 * 100 for x in ypoints1]
    ypoints2 = [x / 1891 * 100 for x in ypoints2]

    plt.figure(figsize = (10,5.625))

    plt.title("Percentage of ID v.s. Ratio of Noise") 
    plt.xlabel("Ratio of Noise")
    plt.ylabel("Percentage of PrSMs")
    plt.plot(xpoints, ypoints1)
    plt.plot(xpoints, ypoints2)
    
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=100))
    plt.gca().xaxis.set_major_formatter(mtick.PercentFormatter(xmax=100))


    plt.legend(["Percentage of IDs Identified", "Percentage of IDs correctly identified"])
    plt.savefig("temp.png", dpi=300)
    


if __name__ == "__main__":
    main()