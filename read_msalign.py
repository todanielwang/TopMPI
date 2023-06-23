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

import os 
from spec_header import SpecHeader
from spec_peak import SpecPeak
from spectrum import Spectrum

def _read_header(spec_lines):
  spec_id = -1
  spec_scan = -1
  retention_time = -1
  activation = -1
  ms_one_id = -1
  ms_one_scan = -1
  mono_mz = -1
  charge = -1
  mono_mass = -1
  inte = -1

  for i in range(len(spec_lines)):
    line = spec_lines[i]
    mono = line.split('=')
    if("ID" == line[0:2]):
      spec_id = int(mono[1])
    if("SCANS" in line):
      spec_scan = int(mono[1])
    if("RETENTION_TIME" in line):
      retention_time = float(mono[1])
    if("ACTIVATION" in line):
      activation = mono[1]
    if("MS_ONE_ID" in line):
      ms_one_id = int(mono[1])
    if("MS_ONE_SCAN" in line):
      ms_one_scan = int(mono[1])
    if("PRECURSOR_MZ" in line):
      mono_mz = float(mono[1])
    if("PRECURSOR_CHARGE" in line):
      charge = int(mono[1])
    if("PRECURSOR_MASS" in line):
      mono_mass = float(mono[1])
    if("PRECURSOR_INTENSITY" in line):
      inte = float(mono[1])
  header = SpecHeader.get_header(spec_id, spec_scan, retention_time, 
               activation, ms_one_id, ms_one_scan, mono_mz, 
               charge, mono_mass, inte)
  return header

def _read_peaks(spec_lines):
  exp_line = "PRECURSOR_INTENSITY"
  end_line = "END IONS"
  peak_list = []
  i = 0
  while(exp_line not in spec_lines[i]):
    i = i + 1
  i = i + 1
  while(spec_lines[i] != end_line): 
    mono = spec_lines[i].split('\t')
    mass = float(mono[0])
    intensity = float(mono[1])
    charge = int(mono[2])
    peak = SpecPeak(mass, intensity, charge)
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
  print("DO NOTE THAT THIS VERSION SKIPS THE MS1 SCANS")
  file = open(filename)
  all_lines = file.readlines()
  all_lines = [x.strip() for x in all_lines] 
  file.close()
  ## Assign file name to header
  filename = os.path.basename(filename)
  spec_list = []
  begin_idx = _get_begin_index(all_lines, 0)
  while (begin_idx < len(all_lines)):
    end_idx = _get_end_index(all_lines, begin_idx)
    spec_lines = all_lines[begin_idx:end_idx +1]
    if (_get_level_one(all_lines, begin_idx, end_idx)):
      begin_idx = end_idx + 1
      continue
    begin_idx = end_idx + 1
    if begin_idx >= len(all_lines):
      break
    
    spec = _parse_spectrum(spec_lines)
    spec.header.file_name = filename
    spec_list.append(spec)
  return spec_list

def write_spec_file(filename, spec_list):
  with open(filename + "_modified", "w") as outputfile:
    for spectrum in spec_list:
      outputfile.write("BEGIN IONS\n")
      outputfile.write("ID=" + str(spectrum.header.spec_id) + "\n")
      outputfile.write("FRACTION_ID=0\nFILE_NAME=/home/daniel/Desktop/datafiles/RealData/yeast/study4/20220718_Yeast_DIA_720-800_75cmCapil_300nL_30kV_01.mzML\n")
      outputfile.write("SCANS=" + str(spectrum.header.spec_scan) + "\n")
      outputfile.write("RETENTION_TIME=" + str(spectrum.header.retention_time) + "\n")
      outputfile.write("LEVEL=2\n")
      outputfile.write("ACTIVATION=" + str(spectrum.header.activation) + "\n")
      outputfile.write("MS_ONE_ID=" + str(spectrum.header.ms_one_id) + "\n")
      outputfile.write("MS_ONE_SCAN=" + str(spectrum.header.ms_one_scan) + "\n")
      outputfile.write("PRECURSOR_MZ=" + str(spectrum.header.mono_mz) + "\n")
      outputfile.write("PRECURSOR_CHARGE=" + str(spectrum.header.charge) + "\n")
      outputfile.write("PRECURSOR_MASS=" + str(spectrum.header.mono_mass) + "\n")
      outputfile.write("PRECURSOR_INTENSITY=" + str(spectrum.header.inte) + "\n")
      for peak in spectrum.peak_list:
        outputfile.write(str(peak.mass) + "\t" + str(peak.intensity) + "\t" + str(peak.charge) + "\n")
      outputfile.write("END IONS\n\n")
