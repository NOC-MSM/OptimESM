#!/bin/bash

set -euo pipefail
# ================================================================
# run_IPSL-CM6-ESMCO2_pipeline.sh
#
# Description: Run IPSL-CM6-ESMCO2 NEMO Pipeline in current process.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-11-11   
# ================================================================

# -- Input arguments to NEMO Pipeline -- #
# Define filepaths:
config_file=mlotst/config_IPSL_esm-hist.toml
log_file=IPSL-CM6-ESMCO2_mlotst_pipeline.log

# Run multiple pipelines:
l_multi=false

# Define Experiment IDs [l_multi=true] -> esm-up2p0-gwl
# exp_ids=("esm-up2p0-gwl1p5" "esm-up2p0-gwl2p0" "esm-up2p0-gwl3p0" "esm-up2p0-gwl4p0")

# Define Experiment IDs [l_multi=true] -> esm-up2p0-gwl-dn
# exp_ids=("" "")


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
