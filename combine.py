import read_msalign
import sys
import csv

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file and the feature file")
    
    spec_list = read_msalign.read_spec_file(args[0])




if __name__ == "__main__":
    main()