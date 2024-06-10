# There are 2 inputs to this function, one the directory to the proteoform.js files, usually under topsearchfolder_html/toppic_proteoform_cutoff/data_js/proteoform/
#The second is the exact protein accession that you are looking for
#example: I want to look at ecoli.mzml and protein 123, the command I use would be
#python3 ubq.py DIRECTORY_TO_THE_NEXT_FOLDER/ecoli_html/toppic_proteoform_cutoff/data_js/proteoforms/ "123"


import pandas as pd
import json
import re
import os
import argparse

class SpecHeader:
  def __init__(self, frac_id, file_name, spec_id, title, spec_scan, retention_time, 
               level, ms_one_id, ms_one_scan, pre_window_begin, pre_window_end, 
               activation, pre_mz_list, pre_charge_list, pre_mass_list, pre_inte_list, pre_id_list):
    self.frac_id = frac_id
    self.file_name = file_name
    self.spec_id = spec_id
    self.title = title
    self.spec_scan = spec_scan
    self.retention_time = retention_time
    self.level = level
    self.ms_one_id = ms_one_id
    self.ms_one_scan = ms_one_scan
    self.pre_window_begin = pre_window_begin
    self.pre_window_end = pre_window_end
    self.activation = activation
    self.pre_mz_list = pre_mz_list
    self.pre_charge_list = pre_charge_list
    self.pre_mass_list = pre_mass_list
    self.pre_inte_list = pre_inte_list
    self.pre_id_list = pre_id_list
  
  @classmethod
  def get_header(cls, frac_id, file_name, spec_id, title, spec_scan, retention_time, 
               level, ms_one_id, ms_one_scan, pre_window_begin, pre_window_end, 
               activation, pre_mz_list, pre_charge_list, pre_mass_list, pre_inte_list, pre_id_list):
    return cls(frac_id, file_name, spec_id, title, spec_scan, retention_time, 
               level, ms_one_id, ms_one_scan, pre_window_begin, pre_window_end, 
               activation, pre_mz_list, pre_charge_list, pre_mass_list, pre_inte_list, pre_id_list)
  
class SpecPeak:
  def __init__(self, mass, intensity, charge, ecscore):
    self.mass = mass
    self.intensity = intensity
    self.charge = charge
    self.ecscore = ecscore

  def __eq__(self, __value: object) -> bool:
    if isinstance(__value, self.__class__):
      if self.mass == __value.mass and self.charge == __value.charge and self.ecscore == __value.ecscore:
          return True
    return False
  
class Spectrum:
  def __init__(self, header, peak_list):
    self.header = header
    self.peak_list = peak_list

  def __del__(self):
    del self.header
    del self.peak_list

  @classmethod
  def get_spec(cls, header, peak_list):
    return cls(header, peak_list)
  
def _read_header(spec_lines):
  frac_id = -1
  file_name = ""
  spec_id = -1
  title = ""
  spec_scan = -1
  retention_time = -1
  level = -1
  ms_one_id = -1
  ms_one_scan = -1
  pre_window_begin = -1
  pre_window_end = -1
  activation = -1
  pre_mz_list = -1
  pre_charge_list = -1
  pre_mass_list = -1
  pre_inte_list = -1
  pre_id_list = -1

  for i in range(len(spec_lines)):
    line = spec_lines[i]
    mono = line.split('=')
    if("FRACTION_ID" in line):
      frac_id = int(mono[1])
    if("FILE_NAME" in line):
      file_name = mono[1]
    if("SPECTRUM_ID" in line):
      spec_id = int(mono[1])
    if("TITLE" in line):
      title = mono[1]
    if("SCANS" in line):
      spec_scan = int(mono[1])
    if("RETENTION_TIME" in line):
      retention_time = float(mono[1])
    if ("LEVEL" in line):
      level = int(mono[1])
    if("MS_ONE_ID" in line):
      ms_one_id = int(mono[1])
    if("MS_ONE_SCAN" in line):
      ms_one_scan = int(mono[1])
    if ("PRECURSOR_WINDOW_BEGIN" in line):
      pre_window_begin = float(mono[1])
    if ("PRECURSOR_WINDOW_END" in line):
      pre_window_end = float(mono[1])
    if("ACTIVATION" in line):
      activation = mono[1]
    if("PRECURSOR_MZ" in line):
      pre_mz_list = mono[1].split(":")
    if("PRECURSOR_CHARGE" in line):
      pre_charge_list = mono[1].split(":")
    if("PRECURSOR_MASS" in line):
      pre_mass_list = mono[1].split(":")
    if("PRECURSOR_INTENSITY" in line):
      pre_inte_list = mono[1].split(":")
    if("PRECURSOR_FEATURE_ID" in line):
      pre_id_list = mono[1].split(":")
  header = SpecHeader.get_header(frac_id, file_name, spec_id, title, spec_scan, retention_time, 
               level, ms_one_id, ms_one_scan, pre_window_begin, pre_window_end, 
               activation, pre_mz_list, pre_charge_list, pre_mass_list, pre_inte_list, pre_id_list)
  return header

def _read_peaks(spec_lines):
  exp_line = "PRECURSOR_FEATURE_ID"
  end_line = "END IONS"
  peak_list = []
  i = 0
  while(i < len(spec_lines) and exp_line not in spec_lines[i]):
    i = i + 1
  i = i + 1
  while(spec_lines[i] != end_line): 
    mono = spec_lines[i].split('\t')
    mass = float(mono[0])
    intensity = float(mono[1])
    charge = int(mono[2])
    ecscore = float(mono[3])
    peak = SpecPeak(mass, intensity, charge, ecscore)
    peak_list.append(peak)
    i = i + 1
  return peak_list

def _parse_spectrum(spec_lines):
  header = _read_header(spec_lines)
  peak_list = _read_peaks(spec_lines)
  spec = Spectrum.get_spec(header, peak_list)
  return spec

def _get_begin_index(all_lines, begin_idx):
  idx = begin_idx
  while (idx < len(all_lines) and "BEGIN IONS" not in all_lines[idx]):
    idx = idx + 1
  return idx

def _get_end_index(all_lines, begin_idx):
  idx = begin_idx
  while (idx < len(all_lines) and "END IONS" not in all_lines[idx]):
    idx = idx + 1
  return idx

def _get_level_one(all_lines, begin_idx, end_idx):
  for idx in range(begin_idx, end_idx): 
    if ("LEVEL=1" in all_lines[idx]):
      return True
  return False

def read_spec_file(filename):
  file = open(filename)
  all_lines = file.readlines()
  all_lines = [x.strip() for x in all_lines] 
  file.close()
  ## Assign file name to header
  spec_list = []
  begin_idx = _get_begin_index(all_lines, 0)
  while (begin_idx < len(all_lines)):
    end_idx = _get_end_index(all_lines, begin_idx)
    spec_lines = all_lines[begin_idx:end_idx +1]
    if (len(spec_lines) < 5):
      break
    if (_get_level_one(all_lines, begin_idx, end_idx)):
      begin_idx = end_idx + 1
      continue
    spec = _parse_spectrum(spec_lines)
    spec_list.append(spec)
    begin_idx = end_idx + 1
  return spec_list


def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1:2] = map(int, parts[1::2])
    return parts


def main():

    # Create the parser
    parser = argparse.ArgumentParser(description="This program will take five positional inputs")

    # Add positional arguments
    parser.add_argument('file', help='The ms2.msalign file')
    parser.add_argument('protein', help='The protein accession as a string, in the form "protein"')
    parser.add_argument('rtbegin', help='The begining of the rt range in minutes')
    parser.add_argument('rtend', help='The end of the rt range in minutes')
    # parser.add_argument('charge', help='The charge state of the precursor')


    # Add optional argument with '-w' flag
    # parser.add_argument('-p', '--precursor', help='An optional flag, use if you would like to limit to sa specific precursor, pass in the MS2 scan number of the specific precursor', required=False)

    # Parse arguments
    args = parser.parse_args()

    try:
        spec_list = read_spec_file(args.file)
    except:
        print("Wrong directory, file not found!")
        exit(1)

    spec_dict = {}
    for spec in spec_list:
        spec_dict[str(spec.header.spec_scan)] = spec

    jsfolder = args.file.rsplit("_", 1)[0] + "_html/toppic_proteoform_cutoff/data_js/proteoforms/"
    outputlist = []
    for filename in sorted(os.listdir(jsfolder), key=numericalSort):
        try:
            with open(str(jsfolder) + filename) as file:
                file.readline()
                toppic = json.loads(file.read())
        except FileNotFoundError:
            print("Wrong directory, file not found!")
            exit(1)
        
        if not toppic["compatible_proteoform"]["sequence_name"] == str(args.protein):
            continue

        dictlist = []
        prsm = toppic["compatible_proteoform"]["prsm"]
        if (int(toppic["compatible_proteoform"]["prsm_number"]) == 1):
            prsm = [prsm]

        # Find the best prsm within the rt range and correct charge state until run out
        rtmin = float(spec_dict[str(prsm[0]["ms"]["ms_header"]["scans"])].header.retention_time) / 60
        # charge = int(spec_dict[str(prsm[0]["ms"]["ms_header"]["scans"])].header.pre_charge_list[0])
        while ((rtmin < float(args.rtbegin) or rtmin > float(args.rtend))):
          prsm = prsm[1:]
          if (len(prsm) == 0):
             break
          rtmin = float(spec_dict[str(prsm[0]["ms"]["ms_header"]["scans"])].header.retention_time) / 60
          # charge = int(spec_dict[str(prsm[0]["ms"]["ms_header"]["scans"])].header.pre_charge_list[0])

        if (len(prsm) == 0):
          continue
        
        prsm = prsm[0]
        
        # tol = (10 * float(args.mass)) / 1e6
        # if (tol < 0.01):
        #   tol = 0.01
        # if not (abs(float(spec_dict[str(prsm["ms"]["ms_header"]["scans"])].header.pre_mass_list[0]) - float(args.mass)) <= tol):
        #    continue
        
        peak_list = prsm["ms"]["peaks"]["peak"]
        for peak in peak_list:
            if "matched_ions" in peak:
                dictlist.append(peak)

        startPos = int(prsm["annotated_protein"]["annotation"]["first_residue_position"])

        for dict in dictlist:
            matched_ion_num = int(dict["matched_ions_num"])
            if matched_ion_num > 1:
                for idx in range(0, matched_ion_num):
                    newdict = {}
                    newdict["Proteoform ID"] = filename
                    for key in dict.keys():
                        if not key == "matched_ions":
                            newdict[key] = dict[key]
                        else:
                            for smallkey in dict[key]["matched_ion"][idx].keys():
                                if smallkey == "ion_position":
                                    newdict[smallkey] = int(dict[key]["matched_ion"][idx][smallkey]) + startPos
                                else:
                                    newdict[smallkey] = dict[key]["matched_ion"][idx][smallkey]
                    outputlist.append(newdict)
            else:
                newdict = {}
                newdict["Proteoform ID"] = filename
                for key in dict.keys():
                    if not key == "matched_ions":
                        newdict[key] = dict[key]
                    else:
                        for smallkey in dict[key]["matched_ion"].keys():
                            if smallkey == "ion_position":
                                newdict[smallkey] = int(dict[key]["matched_ion"][smallkey]) + startPos
                            else:
                                newdict[smallkey] = dict[key]["matched_ion"][smallkey]
                outputlist.append(newdict)
    
    df = pd.DataFrame(outputlist)

    df.to_csv("ion list.tsv", sep="\t")

if __name__ == "__main__":
    main()