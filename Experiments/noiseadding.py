import read_msalign
import pandas as pd
import sys
import random
import math
import copy
import util

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file, the prsm_single.tsv")
    
    random.seed(0)
    
    raw_spec_list = read_msalign.read_spec_file(args[0])

    prsm_df = util.read_tsv(args[1])

    prsm_df = prsm_df[prsm_df["E-value"] < 0.01]

    prsm_df = prsm_df[~prsm_df['Protein accession'].str.contains('DECOY')].reset_index(drop=True)

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
                        and (len(util.removePeaks(copy.deepcopy(nonmulti_dict[x].peak_list), proteoform)) > length * topnoise)]
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

        new_noise_peak_list= util.removePeaks(noise_peak_list, proteoform)

        outputspec = copy.deepcopy(spec)

        pre_mz_list = outputspec.header.pre_mz_list[:1] + [spec0.header.pre_mz_list[0]]
        pre_charge_list = outputspec.header.pre_charge_list[:1] + [spec0.header.pre_charge_list[0]]
        pre_mass_list = outputspec.header.pre_mass_list[:1] + [spec0.header.pre_mass_list[0]]
        pre_inte_list = outputspec.header.pre_inte_list[:1] + [spec0.header.pre_inte_list[0]]
        pre_id_list = outputspec.header.pre_id_list[:1] + [spec0.header.pre_id_list[0]]


        outputspec.header.pre_mz_list = pre_mz_list
        outputspec.header.pre_charge_list = pre_charge_list
        outputspec.header.pre_mass_list = pre_mass_list
        outputspec.header.pre_inte_list = pre_inte_list
        outputspec.header.pre_id_list = pre_id_list

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

        new_noise_peak_list= util.removePeaks(noise_peak_list, proteoform)

        outputspec = copy.deepcopy(spec)

        pre_mz_list = outputspec.header.pre_mz_list[:1] + [spec1.header.pre_mz_list[0]]
        pre_charge_list = outputspec.header.pre_charge_list[:1] + [spec1.header.pre_charge_list[0]]
        pre_mass_list = outputspec.header.pre_mass_list[:1] + [spec1.header.pre_mass_list[0]]
        pre_inte_list = outputspec.header.pre_inte_list[:1] + [spec1.header.pre_inte_list[0]]
        pre_id_list = outputspec.header.pre_id_list[:1] + [spec1.header.pre_id_list[0]]

        outputspec.header.pre_mz_list = pre_mz_list
        outputspec.header.pre_charge_list = pre_charge_list
        outputspec.header.pre_mass_list = pre_mass_list
        outputspec.header.pre_inte_list = pre_inte_list
        outputspec.header.pre_id_list = pre_id_list

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