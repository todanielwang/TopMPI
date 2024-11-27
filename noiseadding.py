import read_msalign
import pandas as pd
import numpy as np
import sys
import re
import random
import math
import copy

def gene_theo_ions(prot_seq):
    left_ions = []
    right_ions = []

    prot_seq = prot_seq.split(".", 1)[1]
    prot_seq = prot_seq.rsplit(".", 1)[0]
  
    acetylation = False

    if "[Acetyl]-" in prot_seq:
        acetylation = True
        prot_seq = prot_seq.replace('[Acetyl]-', '')

    if "(C)[Carbamidomethylation]" in prot_seq:
        prot_seq = prot_seq.replace("(C)[Carbamidomethylation]", "X")


    acetylation_weight = 42.0106
    # c57_weight = 57.021464

    weights = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841, 'G': 57.02146, 
            'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
            'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203, 'T': 101.04768, 'V': 99.06841,
            'W': 186.07931, 'Y': 163.06333, "X": 160.030654}


    if acetylation: 
        left_ions.append(acetylation_weight) 

    idx = 0
    previousweight = 0
    while idx < len(prot_seq):
        if (prot_seq[idx] == "("):
            startidx = idx
            endidx = prot_seq[idx:].find(")") + idx
            modS = endidx + 1
            modE = prot_seq[modS:].find("]") + modS
            weight = float(prot_seq[modS + 1:modE])
            frags = prot_seq[startidx + 1:endidx]
            for fragIdx in range(0, len(frags)):
                baseweight = previousweight
                left_ions.append(baseweight + weights[frags[fragIdx]] + weight)
                if fragIdx != len(frags) - 1:
                    left_ions.append(baseweight + weights[frags[fragIdx]])
                previousweight = left_ions[-1]
            idx = modE + 1
        else:
            left_ions.append(previousweight + weights[prot_seq[idx]])
            previousweight = left_ions[-1]
            idx += 1

    idx = len(prot_seq) - 1
    previousweight = 0
    while idx >= 0:
        if (prot_seq[idx] == "]"):
            startidx = prot_seq[:idx].rfind("(")
            endidx = prot_seq[:idx].rfind(")")
            modS = endidx + 1
            modE = idx
            weight = float(prot_seq[modS + 1:modE])
            frags = prot_seq[startidx + 1:endidx]
            for fragIdx in range(0, len(frags)):
                baseweight = previousweight
                right_ions.append(baseweight + weights[frags[len(frags) - fragIdx - 1]] + weight)
                if fragIdx != len(frags) - 1:
                    right_ions.append(baseweight + weights[frags[len(frags) - fragIdx - 1]])
                previousweight = right_ions[-1]
            idx = startidx - 1
        else:
            right_ions.append(previousweight + weights[prot_seq[idx]])
            previousweight = right_ions[-1]
            idx -= 1

    return left_ions, right_ions

def get_modified_fragments(mass_list, shift):
    mod_mass_list = [x + shift for x in mass_list]
    return mod_mass_list

def remove(peak_list, mass_list, shift):
    frags = get_modified_fragments(mass_list, shift)
    count = 0
    for idx, peak in reversed(list(enumerate(peak_list))):
        peak_mass = peak.mass
        for j in range(len(frags)):
            frag_mass = frags[j]
            tol = (10 * frag_mass) / 1e6
            if (tol < 0.01):
                tol = 0.01
            if (abs(peak_mass - frag_mass) <= tol):
                del peak_list[idx]
                count += 1
                break
    return count


def removePeaks(peak_list, prot_sequence):
    bions, yions = gene_theo_ions(prot_sequence)
    Proton = 1.00727647
    H = 1.007825035
    O = 15.99491463
    CO = 12.0000 + O
    NH3 = 14.003074 + H + H + H
    H2O = H + H + O
    
    count = 0
    # b -ion
    count += remove(peak_list, bions, 0.0)
    # y -ion
    count += remove(peak_list, yions, 19.0184-Proton)
    # # b - 1
    # count += remove(peak_list, bions, -Proton)
    # # y - 1
    # count += remove(peak_list, yions, 19.0184-Proton-Proton)
    # # b + 1
    # count += remove(peak_list, bions, Proton)
    # # y + 1
    # count += remove(peak_list, yions, 19.0184)

    # if (count > 0):
    #     print(str(count) + " peaks was removed from this spectra")

    return peak_list

def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts

""" def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file and the prsm_single.tsv")
    
    random.seed(0)
    
    raw_spec_list = read_msalign.read_spec_file(args[0])

    prsm_df = pd.read_csv(args[1], delimiter="\t")

    scan_list = prsm_df["Scan(s)"].tolist()

    spec_list = []

    for spec in raw_spec_list:
        if int(spec.header.spec_scan) in scan_list:
            spec_list.append(spec)

    for spec in spec_list:
        for counter in range(0, math.ceil(len(spec.peak_list)*1.0)):
            peaklist = spec_list[random.randint(0, len(spec_list) - 1)].peak_list
            peak = peaklist[random.randint(0, len(peaklist) - 1)]
            while (peak in spec.peak_list):
                peaklist = spec_list[random.randint(0, len(spec_list) - 1)].peak_list
                peak = peaklist[random.randint(0, len(peaklist) - 1)]
            spec.peak_list.append(peak)

    read_msalign.write_spec_file(args[0], spec_list) """


def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file, the prsm_single.tsv")
    
    random.seed(0)
    
    raw_spec_list = read_msalign.read_spec_file(args[0])

    prsm_df = pd.read_csv(args[1], delimiter="\t")

    scan_list = prsm_df["Scan(s)"].tolist()

    prsms = []

    for spec in raw_spec_list:
        if int(spec.header.spec_scan) in scan_list:
            prsms.append(spec)

    nonmulti_dict = {}
    toplength = 0
    for spec in prsms:
        if (len(spec.header.pre_inte_list) == 1 or float(spec.header.pre_inte_list[0]) > 0.85 * sum(map(float, spec.header.pre_inte_list))):
            nonmulti_dict[spec.header.spec_scan] = spec
            if len(spec.peak_list) > toplength:
                toplength = len(spec.peak_list)

    print(len(nonmulti_dict))

    candidate_dict = {}
    topnoise = 2
    for spec in nonmulti_dict.values():
        if (len(spec.peak_list) * topnoise < toplength): #To speed up, since we will need length of 2n, if 2n > N for N being max length, then candidate list would be empty anyway
            candidate_dict[spec.header.spec_scan] = spec

    outputList = []
    for i in range(0, 11):
        outputList.append([])

    max_mass = 0
    count0 = 0
    count1 = 0
    for scan in candidate_dict.keys():
        spec = candidate_dict[scan]
        protein = prsm_df[prsm_df["Scan(s)"] == int(scan)].iloc[0]["Protein accession"]
        proteoform = prsm_df[prsm_df["Scan(s)"] == int(scan)].iloc[0]["Proteoform"]
        length = len(spec.peak_list)
        diff = 20
        possibleList = [x for x in nonmulti_dict.keys() if (x != scan) and (abs(float(nonmulti_dict[x].header.pre_mz_list[0]) - float(spec.header.pre_mz_list[0])) < diff)
                        and (len(nonmulti_dict[x].peak_list) > math.ceil(length * topnoise)) and (prsm_df[prsm_df["Scan(s)"] == int(x)].iloc[0]["Protein accession"] != protein) 
                        and (len(removePeaks(copy.deepcopy(nonmulti_dict[x].peak_list), proteoform)) > length * topnoise)]
        # and (float(nonmulti_dict[x].header.pre_mz_list[0]) - float(spec.header.pre_mz_list[0]) < diff)

        if (len(possibleList) == 0):
            continue
        
        specscan0 = 0
        specscan1 = 0

        spec0list = []
        spec1list = []
        for specscan in possibleList:
            if int(prsm_df[prsm_df["Scan(s)"] == int(specscan)].iloc[0]["#unexpected modifications"]) == 0:
                spec0list.append(specscan)
            elif int(prsm_df[prsm_df["Scan(s)"] == int(specscan)].iloc[0]["#unexpected modifications"]) == 1:
                spec1list.append(specscan)

        if len(spec0list) == 0 or len(spec1list) == 0:
            continue

        if int(prsm_df[prsm_df["Scan(s)"] == int(scan)].iloc[0]["#unexpected modifications"]) == 0:
            count0 += 1
        elif int(prsm_df[prsm_df["Scan(s)"] == int(scan)].iloc[0]["#unexpected modifications"]) == 1:
            count1 += 1

        specscan0 = random.choice(spec0list)
        specscan1 = random.choice(spec1list)
        
        spec0 = nonmulti_dict[specscan0]
        spec1 = nonmulti_dict[specscan1]

        mass_diff = abs(float(spec0.header.pre_mass_list[0]) - float(spec.header.pre_mass_list[0]))

        if (mass_diff > max_mass):
            max_mass = mass_diff

        spec.header.title = str(specscan0)

        noise_peak_list = copy.deepcopy(spec0.peak_list)

        new_noise_peak_list= removePeaks(noise_peak_list, proteoform)

        outputspec = copy.deepcopy(spec)

        outputList[0].append(copy.deepcopy(outputspec))

        step = math.ceil(length * 0.2)
        for idx in range(1, 11):
            random.shuffle(new_noise_peak_list)

            noiselist = new_noise_peak_list[0:step]

            new_noise_peak_list = new_noise_peak_list[step:]
            
            outputspec.peak_list = outputspec.peak_list + noiselist
            
            outputList[idx].append(copy.deepcopy(outputspec))

        mass_diff = abs(float(spec1.header.pre_mass_list[0]) - float(spec.header.pre_mass_list[0]))

        if (mass_diff > max_mass):
            max_mass = mass_diff

        spec.header.title = str(specscan1)

        noise_peak_list = copy.deepcopy(spec1.peak_list)

        new_noise_peak_list= removePeaks(noise_peak_list, proteoform)

        outputspec = copy.deepcopy(spec)

        outputList[0].append(copy.deepcopy(outputspec))

        step = math.ceil(length * 0.2)
        for idx in range(1, 11):
            random.shuffle(new_noise_peak_list)

            noiselist = new_noise_peak_list[0:step]

            new_noise_peak_list = new_noise_peak_list[step:]
            
            outputspec.peak_list = outputspec.peak_list + noiselist
            
            outputList[idx].append(copy.deepcopy(outputspec))

        # random_spec = nonmulti_dict[possibleList[random.randint(0, len(possibleList) - 1)]]

        # other_scan = random_spec.header.spec_scan

        # mass_diff = abs(float(random_spec.header.pre_mass_list[0]) - float(spec.header.pre_mass_list[0]))

        # if (mass_diff > max_mass):
        #     max_mass = mass_diff

        # spec.header.title = str(other_scan)

        # noise_peak_list = copy.deepcopy(random_spec.peak_list)

        # new_noise_peak_list= removePeaks(noise_peak_list, proteoform)

        # outputspec = copy.deepcopy(spec)

        # outputList[0].append(copy.deepcopy(outputspec))

        # step = math.ceil(length * 0.2)
        # for idx in range(1, 11):
        #     random.shuffle(new_noise_peak_list)

        #     noiselist = new_noise_peak_list[0:step]

        #     new_noise_peak_list = new_noise_peak_list[step:]
            
        #     outputspec.peak_list = outputspec.peak_list + noiselist
            
        #     outputList[idx].append(copy.deepcopy(outputspec))

    print("We now have " + str(len(outputList[0])) + " spectra left, and the mass difference is " + str(max_mass))
    print(count0, count1)

    for idx in range(0, 11):
        filenamelist = args[0].rsplit(".", 1)
        read_msalign.write_spec_file(filenamelist[0] + str(idx * 0.2) + filenamelist[1], read_msalign.sortScans(outputList[idx]))
    

if __name__ == "__main__":
    main()