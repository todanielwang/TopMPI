import read_msalign
import sys
import numpy as np
import pandas as pd

def main():
    args = sys.argv[1:]
    
    speclist = read_msalign.read_spec_file(args[0])

    sorted = read_msalign.sortScans(speclist)

    read_msalign.write_spec_file(args[0], sorted)

    

if __name__ == "__main__":
    main()