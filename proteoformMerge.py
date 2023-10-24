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
        
    left_counter = 0
    right_counter = 0
    while left_counter < len(left):            
        if (left[left_counter]["Scan(s)"]) < right[right_counter]["Scan(s)"]:
            left_counter += 1
        elif (left[left_counter]["Scan(s)"] > right[right_counter]["Scan(s)"]):
            right_counter += 1
        else:
            if (left[left_counter]["Protein accession"] != right[right_counter]["Protein accession"] or \
                abs(float(left[left_counter]["Proteoform mass"]) - float(right[right_counter]["Proteoform mass"])) > 1):
                del left[left_counter]
                left_counter -= 1
            left_counter += 1
            right_counter += 1

    print(len(left))

    with open('prsm_result.tsv', 'w', newline='\n') as output_file:
        dict_writer = csv.DictWriter(output_file, left[0].keys(), delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(left)

if __name__ == "__main__":
    main()