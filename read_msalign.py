##Copyright (c) 2014 - 2020, The Trustees of Indiana University.
##
##Licensed under the Apache License, Version 2.0 (the "License");
##you may not use this file except in compliance with the License.
##You may obtain a copy of the License at
##
##    http://www.apache.org/licenses/LICENSE-2.0
##
##Unless required by applicable law or agreed to in writing, software
##distributed under the License is distributed on an "AS IS" BASIS,
##WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##See the License for the specific language governing permissions and
##limitations under the License.

#!/usr/bin/python3

from spec_header import SpecHeader
from spec_peak import SpecPeak
from spectrum import Spectrum
import copy

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

def write_spec_file(filename, spec_list):
  filenamelist = filename.split(".")
  with open(filenamelist[0] + "_modified." + filenamelist[1], "w") as outputfile:
    for spectrum in spec_list:
      outputfile.write("BEGIN IONS\n")
      outputfile.write("FRACTION_ID=" + str(spectrum.header.frac_id) + "\n")
      outputfile.write("FILE_NAME=" + spectrum.header.file_name + "\n")
      outputfile.write("SPECTRUM_ID=" + str(spectrum.header.spec_id) + "\n")
      outputfile.write("TITLE=" + spectrum.header.title + "\n")
      outputfile.write("SCANS=" + str(spectrum.header.spec_scan) + "\n")
      outputfile.write("RETENTION_TIME=" + str(spectrum.header.retention_time) + "\n")
      outputfile.write("LEVEL=" + str(spectrum.header.level) + "\n")
      outputfile.write("MS_ONE_ID=" + str(spectrum.header.ms_one_id) + "\n")
      outputfile.write("MS_ONE_SCAN=" + str(spectrum.header.ms_one_scan) + "\n")
      outputfile.write("PRECURSOR_WINDOW_BEGIN=" + str(spectrum.header.pre_window_begin) + "\n")
      outputfile.write("PRECURSOR_WINDOW_END=" + str(spectrum.header.pre_window_end) + "\n")
      outputfile.write("ACTIVATION=" + str(spectrum.header.activation) + "\n")
      outputfile.write("PRECURSOR_MZ=" + ":".join(spectrum.header.pre_mz_list) + "\n")
      outputfile.write("PRECURSOR_CHARGE=" + ":".join(spectrum.header.pre_charge_list) + "\n")
      outputfile.write("PRECURSOR_MASS=" + ":".join(spectrum.header.pre_mass_list) + "\n")
      outputfile.write("PRECURSOR_INTENSITY=" + ":".join(spectrum.header.pre_inte_list) + "\n")
      outputfile.write("PRECURSOR_FEATURE_ID=" + ":".join(spectrum.header.pre_id_list) + "\n")
      for peak in spectrum.peak_list:
        outputfile.write(str(peak.mass) + "\t" + str(peak.intensity) + "\t" + str(peak.charge) + "\t" + str(peak.ecscore) + "\n")
      outputfile.write("END IONS\n\n")

def switchPrecursors(spec_list):
  for spec in spec_list:
    if len(spec.header.pre_charge_list) < 2:
      continue
    pre_mz_list = spec.header.pre_mz_list
    pre_charge_list = spec.header.pre_charge_list
    pre_mass_list = spec.header.pre_mass_list
    pre_inte_list = spec.header.pre_inte_list
    pre_id_list = spec.header.pre_id_list

    pre_mz_list[0], pre_mz_list[1] = pre_mz_list[1], pre_mz_list[0]
    pre_charge_list[0], pre_charge_list[1] = pre_charge_list[1], pre_charge_list[0]
    pre_mass_list[0], pre_mass_list[1] = pre_mass_list[1], pre_mass_list[0]
    pre_inte_list[0], pre_inte_list[1] = pre_inte_list[1], pre_inte_list[0]
    pre_id_list[0], pre_id_list[1] = pre_id_list[1], pre_id_list[0]

    spec.header.pre_mz_list = pre_mz_list
    spec.header.pre_charge_list = pre_charge_list
    spec.header.pre_mass_list = pre_mass_list
    spec.header.pre_inte_list = pre_inte_list
    spec.header.pre_id_list = pre_id_list

def sortScans(spec_list):
    sort_dict = {}
    for spec in spec_list:
      scan = spec.header.spec_scan % 100000
      temp_spec = copy.deepcopy(spec)
      while (scan in sort_dict.keys()):
        scan += 100000
        temp_spec.header.spec_id = scan
      temp_spec.header.spec_scan = scan
      temp_spec.header.spec_id = scan
      sort_dict[scan] = temp_spec

    return dict(sorted(sort_dict.items())).values()