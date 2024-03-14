import sys
import pandas as pd
import numpy as np
import csv

def main():
    args = sys.argv[1:]
    


    r1 = pd.read_csv(args[0], sep="\t")
    r2 = pd.read_csv(args[1], sep="\t")
    
    with open(args[0]) as f:
        left = [{k: v for k, v in row.items()}
            for row in csv.DictReader(f, skipinitialspace=True, delimiter="\t")]
    with open(args[1]) as f:
        right = [{k: v for k, v in row.items()}
            for row in csv.DictReader(f, skipinitialspace=True, delimiter="\t")]
            
    signal = 0
    for rowNew in right:
        for rowOld in left:
            if rowNew["Protein accession"] == rowOld["Protein accession"] and abs(float(rowNew["Proteoform mass"]) - float(rowOld["Proteoform mass"])) < 1:
                signal = 1
                break
        if signal == 0:
            left.append(rowNew)
        else:
            signal = 0

    
    combined = pd.concat([r1, r2], ignore_index=True)
    Protein = combined.drop_duplicates(subset=["Protein accession"])
    with open('Proteoform.csv', 'w', newline='\n') as output_file:
        dict_writer = csv.DictWriter(output_file, left[0].keys(), delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(left)
    Protein.to_csv("Protein.csv", sep="\t")
        


if __name__ == "__main__":
    main()