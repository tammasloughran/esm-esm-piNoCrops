#!/bin/bash

# Initialise an ACCESS-ESM Payu run from a CSIRO experiment
# This should be run from the top-level warm-start.sh script, which sets up the $csiro_source etc. environment variables

set -eu
trap "echo Error in warm_start_csiro.sh" ERR

echo "Sourcing restarts from ${csiro_source} / year ${source_year}"

# Load some helper scripts
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source $SCRIPTDIR/utils.sh

# Start year of this run - read from config.yaml
start_year=$(get_payu_start_year)

# =====================================================================

# Setup the restart directory
payu sweep > /dev/null
payu setup --archive > /dev/null
payu_archive=./archive

payu_restart=${payu_archive}/restart000
if [ -d ${payu_restart} ]; then
    echo "ERROR: Restart directory already exists"
    echo "Consider 'payu sweep --hard' to delete all restarts"
    exit 1
fi

mkdir $payu_restart
mkdir $payu_restart/{atmosphere,ocean,ice,coupler}

# Some restart files are marked jan 01, some are dec 31 of the previous year
yearstart="$(printf '%04d' $source_year)0101"
pyearend="$(printf '%04d' $(( source_year - 1 )) )1231"

# Copy in the restart files from the CSIRO run
cp -v $csiro_source/atm/${expname}.astart-${yearstart} $payu_restart/atmosphere/restart_dump.astart

for f in $csiro_source/cpl/*-${pyearend}; do
    cp -v $f $payu_restart/coupler/$(basename ${f%-*})
done

for f in $csiro_source/ocn/*-${pyearend}; do
    cp -v $f $payu_restart/ocean/$(basename ${f%-*})
done

for f in $csiro_source/ice/*-${pyearend}; do
    cp -v $f $payu_restart/ice/$(basename ${f%-*})
done
cp -v $csiro_source/ice/iced.${yearstart} $payu_restart/ice/

# Set the year of each model component to the run start year
$SCRIPTDIR/set_restart_year.sh $start_year

# Cleanup to be ready to run the model
payu sweep
