#!/bin/bash

# Initialise an ACCESS-ESM Payu run from a Payu experiment
# This should be run from the top-level warm-start.sh script, which sets up the $payu_source environment variable

set -eu
trap "echo Error in warm_start_payu.sh" ERR

echo "Sourcing restarts from ${payu_source}"

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

# Copy the source payu restart directory
cp -r "$payu_source" "$payu_restart"

# Set the year of each model component to the run start year
$SCRIPTDIR/set_restart_year.sh $start_year

# Cleanup to be ready to run the model
payu sweep
