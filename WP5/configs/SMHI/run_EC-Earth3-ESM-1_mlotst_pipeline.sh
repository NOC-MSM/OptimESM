#!/bin/bash

set -euo pipefail
# ================================================================
# run_EC-Earth3-ESM-1_pipeline.sh
#
# Description: Run EC-Earth3-ESM-1 NEMO Pipeline in current process.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-11-11   
# ================================================================
set -euo pipefail

# -- Input arguments to NEMO Pipeline -- #
# Define filepaths:
config_file=mlotst/config_EC-Earth3-ESM-1_esm-up2p0-gwl-dn.toml
log_file=EC-Earth3-ESM-1_mlotst_pipeline.log

# -- Python Environment -- #
# Run this script in the env_optimesm conda virtual environment.

# IDs -> esm-up2p0-gwl
exp_ids=("esm-up2p0-gwl1p5" "esm-up2p0-gwl2p0" "esm-up2p0-gwl3p0" "esm-up2p0-gwl4p0" "esm-up2p0-gwl5p0" "esm-up2p0-gwl6p0")

# IDs -> esm-up2p0-gwl-dn
exp_ids=("esm-up2p0-gwl1p5-50y-dn2p0" "esm-up2p0-gwl2p0-200y-dn2p0" "esm-up2p0-gwl2p0-50y-dn1p0" "esm-up2p0-gwl2p0-50y-dn2p0" "esm-up2p0-gwl3p0-50y-dn2p0" "esm-up2p0-gwl4p0-200y-dn2p0" "esm-up2p0-gwl4p0-50y-dn1p0" "esm-up2p0-gwl4p0-50y-dn2p0")

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
