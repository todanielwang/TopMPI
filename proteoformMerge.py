import sys
import pandas as pd
import numpy as np
import csv


def main():

    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in prsm.tsv files")
   



    with open(args[0]) as f:
        left = [{k: v for k, v in row.items()}
            for row in csv.DictReader(f, skipinitialspace=True, delimiter="\t")]
    with open(args[1]) as f:
        right = [{k: v for k, v in row.items()}
            for row in csv.DictReader(f, skipinitialspace=True, delimiter="\t")]
        
    print(len(left))
    

    for left_row in reversed(left):
        for right_row in right:
            if left_row["Protein accession"] == right_row["Protein accession"] and \
                abs(float(left_row["Proteoform mass"]) - float(right_row["Proteoform mass"])) < 5:
                left.remove(left_row)
                break

    print(len(left))

    with open('left.tsv', 'w', newline='\n') as output_file:
        dict_writer = csv.DictWriter(output_file, left[0].keys(), delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(left)

if __name__ == "__main__":
    main()