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
    project=p66
    user=txz599
    export expname=PI-C2C-1p5r42            # Source experiment - PI pre-industrial, HI historical
    export source_year=700          # Change this to create different ensemble members
    export csiro_source=/g/data/$project/$user/archive/ACCESS-ESM1-5/$expname/restart

    # Call the main warm-start script
    scripts/warm-start-csiro.sh

else

    # Payu restart directory to copy
    export payu_source=/scratch/w35/saw562/access-esm/archive/esm-historical/restart001

    # Call the main warm-start script
    scripts/warm-start-payu.sh
fi
