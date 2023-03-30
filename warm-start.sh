#!/bin/bash

# Initialise an ACCESS-ESM Payu run from another experiment
#
# This script sets values specific to the experiment, it then calls the common
# 'warm-start' from the scripts directory which performs the copy and sets the
# start dates for all the component models to the date in config.yaml
set -eu

# Either CSIRO or PAYU
source=CSIRO

if [ $source == "CSIRO" ]; then

    # CSIRO job to copy the warm start from
    project=p73
    user=txz599
    export expname=PI-EDC-01            # Source experiment - PI pre-industrial, HI historical
    export source_year=311          # Change this to create different ensemble members
    export csiro_source=/g/data/$project/archive/CMIP6/ACCESS-ESM1-5/$expname/restart

    # Call the main warm-start script
    scripts/warm-start-csiro.sh

else

    # Payu restart directory to copy
    export payu_source=/scratch/w35/saw562/access-esm/archive/esm-historical/restart001

    # Call the main warm-start script
    scripts/warm-start-payu.sh
fi
