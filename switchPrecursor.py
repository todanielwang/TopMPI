#Map the scans with second feature

import read_msalign
import sys
import re

def main():
    args = sys.argv[1:]
    if (len(args) != 1):
        raise Exception(
            "Please pass in the ms.align file")

    spec_list = read_msalign.read_spec_file(args[0])

    read_msalign.switchPrecursors(spec_list)

    read_msalign.write_spec_file(args[0], spec_list)

    
if __name__ == "__main__":
    main()