#!/bin/bash

# Pre-run script, runs after payu sets up the work directory and before the model is run

# Sets the appropriate land use for the model start date.
# The model doesn't vary the land use itself, this is instead done as a
# pre-processing step - the old land use values are moved to a new STASH code,
# and new land use values are read in from an external file.

source  /etc/profile.d/modules.sh
module use /g/data/hh5/public/modules
module load conda/analysis3

set -eu

# Input land use file
lu_file=$1

# Get the current year from field t2_year of the restart file
year=$(mule-pumf --component fixed_length_header work/atmosphere/restart_dump.astart | sed -n 's/.*t2_year\s*:\s*//p')

# If that year is in the land use file, save a single timestep to a new netcdf file
year_in_file=$(cdo showdate $lu_file)
if [[ $year_in_file =~ $year ]]; then
    cdo selyear,$(( year )) $lu_file work/atmosphere/land_frac.nc
else
    # Set the year in the file to the current year
    cdo setyear,$year $lu_file work/atmosphere/land_frac.nc
fi

# Back up the original restart file
mv work/atmosphere/restart_dump.astart work/atmosphere/restart_dump.astart.orig

# Use the CSIRO script to set the land use
python scripts/update_cable_vegfrac.py \
        -i work/atmosphere/restart_dump.astart.orig \
        -o work/atmosphere/restart_dump.astart \
        -f work/atmosphere/land_frac.nc \
        -n fraction

# If nothing was done, then copy the restart file back.
if [[ ! -e work/atmosphere/restart_dump.astart ]]; then
    echo "Using the cover fractions from the existing restart file."
    cp work/atmosphere/restart_dump.astart.orig work/atmosphere/restart_dump.astart
fi

# Occasionally, ACCESS-ESM1.5 will crash from a grid cell storm. A grid cell storm is numerical
# instability from a coincidence of conditions in just a few grid cells. The workarround is to
# add some small perturbations to the last reastart file to eventually avoid the storm entirely.

module purge
module use ~access/modules
module load pythonlib/umfile_utils

if [[ $year == 127 ]]; then
    echo "Restart file has been perturbed at year $year to avoid a crash."
    scripts/perturbIC.py work/atmosphere/restart_dump.astart -s $year
fi
