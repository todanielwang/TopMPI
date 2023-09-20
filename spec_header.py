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

class SpecHeader:
  def __init__(self, frac_id, file_name, spec_id, title, spec_scan, retention_time, 
               level, ms_one_id, ms_one_scan, pre_window_begin, pre_window_end, 
               activation, pre_mz, pre_charge, pre_mass, pre_inte, pre_feature_id):
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
    self.pre_mz = pre_mz
    self.pre_charge = pre_charge
    self.pre_mass = pre_mass
    self.pre_inte = pre_inte
    self.pre_feature_id = pre_feature_id
  
  @classmethod
  def get_header(cls, frac_id, file_name, spec_id, title, spec_scan, retention_time, 
               level, ms_one_id, ms_one_scan, pre_window_begin, pre_window_end, 
               activation, pre_mz, pre_charge, pre_mass, pre_inte, pre_feature_id):
    return cls(frac_id, file_name, spec_id, title, spec_scan, retention_time, 
               level, ms_one_id, ms_one_scan, pre_window_begin, pre_window_end, 
               activation, pre_mz, pre_charge, pre_mass, pre_inte, pre_feature_id)
