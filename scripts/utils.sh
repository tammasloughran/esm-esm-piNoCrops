#!/bin/bash
#  Copyright 2021 Scott Wales
#
#  \author  Scott Wales <scott.wales@unimelb.edu.au>
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

get_payu_start_year() {
    python <<EOF
import yaml
with open('config.yaml') as f:
    c = yaml.safe_load(f)
    print(c['calendar']['start']['year'])
EOF
}

set_um_start_year() {
    start_year="$1"
    python <<EOF
import f90nml
nml = f90nml.read('atmosphere/namelists')
nml['NLSTCALL']['MODEL_BASIS_TIME'][0] = $start_year
nml.write('atmosphere/namelists', force=True)
EOF
}

set_cice_start_year() {
    start_year="$1"
    python <<EOF
import f90nml
nml = f90nml.read('ice/input_ice.nml')
nml['coupling']['init_date'] = ${start_year}0101 # Experiment start date
nml['coupling']['inidate']   = ${start_year}0101 # First run date
nml.write('ice/input_ice.nml', force=True)
EOF
}
