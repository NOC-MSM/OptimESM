#!/bin/bash

# ================================================================
# run_EC-Earth3-ESM-1_pipeline.sh
#
# Description: Run EC-Earth3-ESM-1 NEMO Pipeline in current process.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-12-15
# ================================================================
set -euo pipefail

# -- Input arguments to NEMO Pipeline -- #
# Define filepaths:
config_file=siarea/config_EC-Earth3-ESM-1_esm-piControl.toml
log_file=EC-Earth3-ESM-1_siarea_pipeline.log

# Run multiple pipelines:
l_multi=false

# Define Experiment IDs [l_multi=true] -> esm-up2p0-gwl
# exp_ids=("esm-up2p0-gwl1p5" "esm-up2p0-gwl2p0" "esm-up2p0-gwl3p0" "esm-up2p0-gwl4p0" "esm-up2p0-gwl5p0" "esm-up2p0-gwl6p0")

# Define Experiment IDs [l_multi=true] -> esm-up2p0-gwl-dn
# exp_ids=("esm-up2p0-gwl1p5-50y-dn2p0" "esm-up2p0-gwl2p0-200y-dn2p0" "esm-up2p0-gwl2p0-50y-dn1p0" "esm-up2p0-gwl2p0-50y-dn2p0" "esm-up2p0-gwl3p0-50y-dn2p0" "esm-up2p0-gwl4p0-200y-dn2p0" "esm-up2p0-gwl4p0-50y-dn1p0" "esm-up2p0-gwl4p0-50y-dn2p0")


# -- Python Environment -- #
# Run this script in the env_optimesm conda virtual environment.

if [ "$l_multi" = false ]; then
    # -- Run NEMO Pipeline CLI -- #
    # nemo_pipeline describe $config_file --log $log_file
    nemo_pipeline run $config_file --log $log_file

else
    # Iterate over all experiment IDs:
    for exp_id in "${exp_ids[@]}"; do
        echo "Running ==> $exp_id"
        # -- Updating Experiment IDs in config.toml -- #
        sed -i "s|esm-[^/_]*|$exp_id|g" $config_file
    
        # -- Run NEMO Pipeline CLI -- #
        # nemo_pipeline describe $config_file --log $log_file
        nemo_pipeline run $config_file --log $log_file
        echo "Completed ==> $exp_id" 
    done
fi
