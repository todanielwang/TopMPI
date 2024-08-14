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

  def __hash__(self):
    return hash(tuple([self.mass, self.intensity, self.charge, self.ecscore]))
  
  def __str__(self):
    return "Test mass:% s charge:% s ecscore:% s\n" % (self.mass, self.charge, self.ecscore)
  
  def __repr__(self):
    return "Test mass:% s charge:% s ecscore:% s\n" % (self.mass, self.charge, self.ecscore)
  
