#!/bin/bash
#  Copyright 2020 Scott Wales
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

# Sets the start year in each model component

source /etc/profile.d/modules.sh
module use /g/data/hh5/public/modules
module load conda/analysis3
module load nco

set -eu
trap "echo Error in set_restart_year.sh" ERR
export UMDIR=~access/umdir

# Load some helper scripts
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source $SCRIPTDIR/utils.sh

# Starting year
start_year=$1

# Set the restart year in the atmosphere and ice namelists
set_um_start_year $start_year
set_cice_start_year $start_year

# Most recent restart directory
payu_restart=$(ls -d ./archive/restart* | sort -t t -k 3 -n | tail -n 1)

echo "Setting restart year in ${payu_restart} to ${start_year}"

if [ ! -d $payu_restart/ocean ]; then
    echo "No restart directory"
    exit 1
fi

# Update ocean start time - read from a text file
cat > $payu_restart/ocean/ocean_solo.res << EOF
    3
    1 1 1 0 0 0
    $start_year 1 1 0 0 0
EOF

# Update atmos start time - field in the restart file
python scripts/update_um_year.py $start_year $payu_restart/atmosphere/restart_dump.astart 2> /dev/null

cat > $payu_restart/atmosphere/um.res.yaml << EOF
end_date: $(printf %04d $start_year)-01-01 00:00:00
EOF

# Clear ice step count
cat > $payu_restart/ice/cice_in.nml << EOF
&setup_nml
istep0=0,
npt=0,
dt=3600,
/
EOF

# Get the number of seconds since the run start date in ice/input_ice.nml
# init_date is the initial date of the experiment
# inidate is the date of the current run
runtime0=$(python <<EOF
import f90nml
import sys
import cftime
t0 = f90nml.read('ice/input_ice.nml')['coupling']['init_date']
y = t0//10000; m = (t0//100)%100; d = t0%100 
t0 = cftime.datetime(y,m,d,calendar='proleptic_gregorian')
t1 = cftime.datetime($start_year,1,1,calendar='proleptic_gregorian')
diff = t1 - t0
print(f"ice init_date {t0}, runtime0 {diff}", file=sys.stderr)
print(int(diff.total_seconds()))
EOF
)

# Put this in the coupling namelist - used by Payu to generate Oasis namcouple file
cat > $payu_restart/ice/input_ice.nml << EOF
&coupling
runtime0=$runtime0
runtime=0
/
EOF

# Set the date in the cice netcdf file
ncatted -a units,time,o,c,"seconds since ${start_year}-01-01 00:00:00" $payu_restart/ice/mice.nc

# Seconds between init_date and inidate
secs_realyr=0

ice_restart=$(ls $payu_restart/ice/iced.*0101)
mv $ice_restart ${ice_restart}.orig

# Set the date in the cice binary restart file
scripts/cicedumpdatemodify.py -i ${ice_restart}.orig -o $payu_restart/ice/iced.${start_year}0101 --istep0=0 --time=${secs_realyr}. --time_forc=0.

rm ${ice_restart}.orig
