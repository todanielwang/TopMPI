import read_msalign
import sys
import os
import json
import re

def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts

def gene_theo_ions(prot_seq):
  left_ions = []
  right_ions = []

  prot_seq = prot_seq.replace('.', '')
  prot_seq = prot_seq.replace('-', '')

  acetylation_weight = 42.0106
  c57_weight = 57.021464
  c58_weight = 58.005479

  weights = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841, 'G': 57.02146, 
    'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203, 'T': 101.04768, 'V': 99.06841,
    'W': 186.07931, 'Y': 163.06333}

  for idx in range(1, len(prot_seq)):
    left_ions.append(left_ions[idx -1] + weights[prot_seq[idx]])
    right_ions.append(right_ions[idx -1] + weights[prot_seq[len(prot_seq)-1-idx]])

  return left_ions, right_ions

def get_modified_fragments(mass_list, shift):
  mod_mass_list = [x + shift for x in mass_list]
  return mod_mass_list

def add_annotation(anno_dict, env_list, mass_list, shift, anno):
    frags = get_modified_fragments(mass_list, shift)

    for idx, env in enumerate(env_list):
      peak_mass = env.mass
      for j in range(len(frags)):
        frag_mass = frags[j]
        tol = (15 * frag_mass) / 1e6
        if (tol < 0.01):
          tol = 0.01
        if (abs(peak_mass - frag_mass) <= tol):
          if not idx in anno_dict:
            anno_dict[idx] = anno
          break

def annotate(env_list, prot_sequence):
    anno_dict = []
    nterm_masses, cterm_masses = gene_theo_ions(prot_sequence)
    Proton = 1.00727647
    H = 1.007825035
    O = 15.99491463
    CO = 12.0000 + O
    NH3 = 14.003074 + H + H + H
    H2O = H + H + O

    # b -ion
    add_annotation(anno_dict, env_list, nterm_masses, 0.0, "B")
    # y -ion
    add_annotation(anno_dict, env_list, cterm_masses, 19.0184-Proton, "Y")
    # c - ion
    add_annotation(anno_dict, env_list, nterm_masses, 18.0344-Proton, "C")
    # z' - ion
    add_annotation(anno_dict, env_list, cterm_masses, 1.9918-Proton+Proton, "Z+1")
    # a -ion
    add_annotation(anno_dict, env_list, nterm_masses, -CO, "A")
    # x -ion
    add_annotation(anno_dict, env_list, cterm_masses, CO, "X")
    # b - waterS
    add_annotation(anno_dict, env_list, nterm_masses, -H2O, "B-Water")
    # y - water
    add_annotation(anno_dict, env_list, cterm_masses, 19.0184-Proton-H2O, "Y-Water")
    # b - ammonia
    add_annotation(anno_dict, env_list, nterm_masses, -NH3, "B-Ammonia")
    # y - ammonia
    add_annotation(anno_dict, env_list, cterm_masses, 19.0184-Proton-NH3, "Y-Ammonia")
    # b - 1
    add_annotation(anno_dict, env_list, nterm_masses, -Proton, "B-1")
    # y - 1
    add_annotation(anno_dict, env_list, cterm_masses, 19.0184-Proton-Proton, "Y-1")
    # b + 1
    add_annotation(anno_dict, env_list, nterm_masses, Proton, "B+1")
    # y + 1
    add_annotation(anno_dict, env_list, cterm_masses, 19.0184, "Y+1")

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file and the js file directory")
    
    spec_list = read_msalign.read_spec_file(args[0])

    spec_list = [spec for spec in spec_list if len(spec.header.pre_mz_list) > 1]

    for spec in spec_list:
        spec.header.pre_mz_list = spec.header.pre_mz_list[1:]
        spec.header.pre_charge_list = spec.header.pre_charge_list[1:]
        spec.header.pre_mass_list = spec.header.pre_mass_list[1:]
        spec.header.pre_inte_list = spec.header.pre_inte_list[1:]
        spec.header.pre_id_list = spec.header.pre_id_list[1:]

    curr_spec = 0
    count = 0
    for x in sorted(os.listdir(args[1]), key=numericalSort):
        if x.endswith(".js"):
            with open(args[1] + x) as file:
                file.readline()
                toppic = json.loads(file.read())
                deleted = False
                while(curr_spec < len(spec_list)):
                    filescan = spec_list[curr_spec].header.spec_scan
                    toppicscan = int(toppic["prsm"]["ms"]["ms_header"]["scans"])
                    if (filescan < toppicscan):
                        curr_spec += 1
                    else:
                        if (filescan > toppicscan):
                            deleted = True
                        break
                if (curr_spec >= len(spec_list)):
                    break
                if (deleted):
                    continue
                peak_list = toppic["prsm"]["ms"]["peaks"]["peak"]
                for idx in range(len(peak_list) -1, -1, -1):
                    if "matched_ions" in peak_list[idx]:
                        del spec_list[curr_spec].peak_list[idx]
                count += 1
    print("Number of scans with peaks removed is {}".format(count))



    read_msalign.write_spec_file(args[0], spec_list)

if __name__ == "__main__":
    main()